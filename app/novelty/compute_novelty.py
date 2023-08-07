import pandas as pd
from datetime import date
import requests
import numpy as np
import traceback

from .known import find_known_results
from .extr_smile_molpro_by_id import mol_to_smile_molpro
from .mol_similarity import find_nearest_neighbors

"""
This script computes the novelty score for a list of results obtained for a 1-H response.
The steps for the ideal workflow are as follows:
1. Distinguish whether the result is a known result or an unknown result.
2. Compute FDA approval status to check if currently under approval or already approved.
3. Compute recency by noting the associated publications of the result.
4. Compute molecular similarity to identify if similar to existing drugs.
 
The end result of this script displays a table with values from different columns and accordingly lists the novelty score as well.
"""


def molecular_sim(known, unknown, message):
    unknown_ids = []
    known_ids = []
    if len(unknown) > 0:
        for drug in unknown:
            s = list(message["results"][drug]["node_bindings"].keys())
            edge_attribute_sn = message["results"][drug]["node_bindings"][s[0]][0]["id"]
            if (
                "PUBCHEM" in edge_attribute_sn
                or "CHEMBL" in edge_attribute_sn
                or "UNII" in edge_attribute_sn
                or "RXNORM" in edge_attribute_sn
                or "UMLS" in edge_attribute_sn
                or not "MONDO" in edge_attribute_sn
            ):
                unknown_ids.append(edge_attribute_sn)
            else:
                unknown_ids.append(
                    message["results"][drug]["node_bindings"][s[1]][0]["id"]
                )

    if len(known) > 0:
        for drug in known:
            s = list(message["results"][drug]["node_bindings"].keys())
            edge_attribute_sn = message["results"][drug]["node_bindings"][s[0]][0]["id"]
            if (
                "PUBCHEM" in edge_attribute_sn
                or "CHEMBL" in edge_attribute_sn
                or "UMLS" in edge_attribute_sn
                or "UNII" in edge_attribute_sn
                or "RXNORM" in edge_attribute_sn
                or not "MONDO" in edge_attribute_sn
            ):
                known_ids.append(edge_attribute_sn)
            else:
                known_ids.append(
                    message["results"][drug]["node_bindings"][s[1]][0]["id"]
                )

    smile_unkown = mol_to_smile_molpro(unknown_ids)
    smile_known = mol_to_smile_molpro(known_ids)

    similarity_map = find_nearest_neighbors(smile_unkown, smile_known, 0, 1)
    return similarity_map


def get_publication_info(pub_id):
    """
    Args: PMID
    Returns: The publication info
    """
    base_url = "https://3md2qwxrrk.us-east-1.awsapprunner.com/publications?pubids="
    request_id = "1df88223-c0f8-47f5-a1f3-661b944c7849"
    full_url = f"{base_url}{pub_id}&request_id={request_id}"
    try:
        response = requests.get(full_url)
        response.raise_for_status()
        response = response.json()
    except Exception:
        response = {
            "_meta": {
                "n_results": 0,
            },
        }
    return response


def sigmoid(x):
    return 1 / (1 + np.exp(x))


def query_id(message):
    for i in message["query_graph"]["nodes"]:
        if "ids" in message["query_graph"]["nodes"][i].keys():
            if message["query_graph"]["nodes"][i]["ids"]:
                known_node = message["query_graph"]["nodes"][i]["categories"][0]
            else:
                unknown_node = message["query_graph"]["nodes"][i]["categories"][0]
        else:
            unknown_node = message["query_graph"]["nodes"][i]["categories"][0]
    if unknown_node in [
        "biolink:ChemicalEntity",
        "biolink:SmallMolecule",
        "biolink:Drug",
    ]:
        chk = 1
    else:
        chk = 0
    return known_node, unknown_node, chk


def recency_function_exp(number_of_publ, age_of_oldest_publ, max_number, max_age):
    """
    Calculates the recency based on number of publications accociated to each drug
    and age of the oldest publication

    Args:
        number_of_publ (float): The current number of publications.
        age_of_oldest_publ (float): The current age of the oldest publication.
        max_number (float): The maximum number of publication: e.g. consider 100 for all drugs.
        max_age (float): The publication with the recent 50 years have been considered.

    Returns:
        float: The recency value of z.
    """
    coef_number = 10
    coef_age = 4
    alter_number_of_publ = sigmoid(coef_number * (number_of_publ / max_number - 0.5))
    alter_age_of_oldest_publ = sigmoid(coef_age * (age_of_oldest_publ / max_age - 0.5))

    if np.isnan(number_of_publ) and not np.isnan(age_of_oldest_publ):
        recency = alter_age_of_oldest_publ
    elif np.isnan(age_of_oldest_publ) and not np.isnan(number_of_publ):
        recency = alter_number_of_publ
    else:
        recency = alter_number_of_publ * alter_age_of_oldest_publ

    return recency


def extracting_drug_fda_publ_date(message, unknown):
    """
    Upon querying, the response is returned as a list containing 10 dictionaries,
    with each dictionary representing the response from an ARA. The function 'extracting_drug_fda_publ_date'
    is designed to extract the drug entity name of each edge. It then checks the EPC attributes of each edge
    to determine the FDA approval status of the drug and the accociated pulications PMID/PMC ID ... .
    And finally extract the publishing date of the publications to get the oldest (Year) date among them.

    Args:
        Dictionary: response for a ARA to a query

    Returns:
        "An DataFrame constructed where each row represents an edge and contains information such as the drug entity
        name, FDA status of the drug, a list of associated publications, the number of associated publications,
        and the oldest publication date (year) linked to each drug."

    """
    attribute_type_id_list_fda = [
        "biolink:FDA_approval_status",
        "biolink:FDA_APPROVAL_STATUS",
    ]
    attribute_type_id_list_pub = [
        "biolink:publications",
        "biolink:Publication",
        "biolink:publication",
    ]
    drug_idx_fda_status = []
    today = date.today()

    res_chk = 1
    query_known, query_unknown, query_chk = query_id(message)
    idi = -1
    for tmp in unknown:
        tmp_res = message["results"][tmp]["analyses"][0]["edge_bindings"]
        for tmp_1 in tmp_res:
            idi += 1
            edge = message["results"][tmp]["analyses"][0]["edge_bindings"][tmp_1][0][
                "id"
            ]
            # edge_list = list(message['knowledge_graph']['edges'].keys())
            # for idx, idi in enumerate(edge_list):
            #     if idx % 20 == 0:
            #         print(f'progressing {idx}')
            # edge = edge_list[idx]
            edge_attribute = message["knowledge_graph"]["edges"][edge]
            # if set(['subject', 'object']).issubset(edge_attribute.keys()):
            if query_chk == 1:
                if (
                    "PUBCHEM" in edge_attribute["subject"]
                    or "CHEMBL" in edge_attribute["subject"]
                    or "UNII" in edge_attribute["subject"]
                    or "RXNORM" in edge_attribute["subject"]
                    or "UMLS" in edge_attribute["subject"]
                    or not "MONDO" in edge_attribute["subject"]
                ):
                    drug_idx = edge_attribute["subject"]
                else:
                    drug_idx = edge_attribute["object"]
                if set(["attributes"]).issubset(edge_attribute.keys()):
                    if len(edge_attribute["attributes"]) > 0:
                        att_type_id = {}
                        fda = []
                        pub = []
                        for i in range(len(edge_attribute["attributes"])):
                            att_type_id[i] = edge_attribute["attributes"][i][
                                "attribute_type_id"
                            ]

                        for key in att_type_id.keys():
                            if att_type_id[key] in attribute_type_id_list_fda:
                                fda.append(key)
                            elif att_type_id[key] in attribute_type_id_list_pub:
                                pub.append(key)

                        if len(fda) > 0:
                            if (
                                edge_attribute["attributes"][fda[0]]["value"]
                                == "FDA Approval"
                            ):
                                fda_status = 0.0
                            else:
                                fda_status = 1.0
                        else:
                            fda_status = None

                        # Publication
                        if len(pub) > 0:
                            publications = edge_attribute["attributes"][pub[0]]["value"]
                            if "|" in publications:
                                publications = publications.split("|")
                            if type(publications) == "str":
                                publications = [publications]

                            # Removal of all publication entries that are links
                            publications = [x for x in publications if "http" not in x]
                            # Removal of all publication entries that are Clinical Trials
                            publications = [
                                x for x in publications if "clinicaltrials" not in x
                            ]
                            number_of_publ = len(publications)

                            if len(publications) > 0:
                                # print(publications)
                                publications_1 = ",".join(publications)
                                try:
                                    response_pub = get_publication_info(publications_1)
                                    if response_pub["_meta"]["n_results"] == 0:
                                        age_oldest = np.nan
                                    else:
                                        publ_year = []
                                        for key in response_pub["results"].keys():
                                            if "not_found" not in key:
                                                publ_year.extend(
                                                    [
                                                        int(
                                                            response_pub["results"][
                                                                key
                                                            ]["pub_year"]
                                                        )
                                                    ]
                                                )
                                        age_oldest = today.year - min(publ_year)
                                except ConnectionError as e:
                                    age_oldest = np.nan
                        else:
                            publications = None
                            number_of_publ = 0.0
                            age_oldest = np.nan
                    drug_idx_fda_status.append(
                        (
                            idi,
                            drug_idx,
                            fda_status,
                            publications,
                            number_of_publ,
                            age_oldest,
                        )
                    )
            else:
                if query_unknown in ["biolink:Gene", "biolink:Protein"]:
                    if (
                        "NCBI" in edge_attribute["subject"]
                        or "GO" in edge_attribute["subject"]
                    ):
                        gene_idx = edge_attribute["subject"]
                    else:
                        gene_idx = edge_attribute["object"]
                    drug_idx_fda_status.append((idi, gene_idx))
                elif query_unknown in ["biolink:Disease", "biolink:Phenotype"]:
                    if "MONDO" in edge_attribute["subject"]:
                        dis_idx = edge_attribute["subject"]
                    else:
                        dis_idx = edge_attribute["object"]
                    drug_idx_fda_status.append((idi, dis_idx))
    if query_chk == 1 and res_chk == 1:
        DF = pd.DataFrame(
            drug_idx_fda_status,
            columns=[
                "edge",
                "drug",
                "fda status",
                "publications",
                "number_of_publ",
                "age_oldest_pub",
            ],
        )
    elif query_chk != 1 and res_chk == 1:
        DF = pd.DataFrame(drug_idx_fda_status, columns=["edge", "result"])
    else:
        DF = pd.DataFrame()
    return DF, query_chk


def extract_results(message, unknown, known):
    results = []
    results_known = []
    kid, ukid = 0, 0
    for idi, i in enumerate(message["results"]):
        if idi in unknown:
            results.append([])
            for idj, j in enumerate(i["analyses"]):
                for idk, k in enumerate(
                    j["edge_bindings"][list(j["edge_bindings"].keys())[0]]
                ):
                    results[ukid].append(k["id"])
            ukid += 1

        elif idi in known:
            results_known.append([])
            for idj, j in enumerate(i["analyses"]):
                for idk, k in enumerate(
                    j["edge_bindings"][list(j["edge_bindings"].keys())[0]]
                ):
                    results_known[kid].append(k["id"])
            kid += 1
    return results, results_known


def result_edge_correlation(results, results_known, df):
    df_res = pd.DataFrame()
    res_known = set()
    res_unknown = set()
    for idi, i in enumerate(results):
        for j in i:
            df_res = pd.concat([df_res, df[df["edge"] == j]])
            res_unknown.add(df.loc[df["edge"] == j, "drug"].iloc[0])

    for idi, i in enumerate(results_known):
        for j in i:
            res_known.add(df.loc[df["edge"] == j, "drug"].iloc[0])
    return df_res, res_unknown, res_known


def novelty_score(fda_status, recency, similarity):
    """
    Calculate the novelty score for each drug entity based on FDA status, recency and similarity of the drug.

    FDA status 0 | 1
        --> 0 to be FDA approved
    0 < recency < 1
        --> 1 to have a high recency where the publications are so new and number of publications is too few.
    0 < similarity < 1
        --> 0 to have a very low molecular structure similarity, so it is novel.

    Args:
        float: fda_status
        float: recency
        float: similarity

    Returns:
        float: novelty_score

    """
    if not np.isnan(recency):
        score = recency
        if not np.isnan(similarity):
            similarity = 1 - similarity
            if similarity > 0.5:
                score = score * (0.73 + similarity)
            if score > 1:
                score = 1
            if fda_status == 0:
                score = score * 0.85
    else:
        if np.isnan(similarity):
            score = 0
        else:
            score = 1 - similarity
    return score


def compute_novelty(message, logger):
    """INPUT: JSON Response with merged annotated results for a 1-H query

    1. load the json file
    2. Give the json to extracting_drug_fda_publ_date(response) function to extract the EPC
    3. Apply the recency function of df, to add a new column as recency to the dataframe
    4. Add a new column to the df as similarity which has random number between 0-1
    5. Now the dataframe df is ready for applying the novelty score on it

    OUTPUT: Pandas DataFrame  with FDA Status, Recency, Similarity and Novelty score per result
    """
    # Step 1
    known, unknown = find_known_results(message)
    #
    # # Step 2

    # start = time.time()
    df, query_chk = extracting_drug_fda_publ_date(message, unknown)
    # print(f"Time to extract fda status and Publication data:{time.time()-start}")
    #         # print(df.head())
    #         # print(query_chk)
    #
    # df.to_excel(f'DATAFRAME.xlsx', header=False, index=False)
    # df = pd.read_excel('DATAFRAME.xlsx', names=['edge', 'drug', 'fda status', 'publications', 'number_of_publ', 'age_oldest_pub'])
    # query_chk = 1

    # #
    # res, res_known = extract_results(mergedAnnotatedOutput, unknown, known)
    #         # print(len(res_known))
    #         # print(len(res))
    #         # # #
    #         df_res, res_unknown, res_known = result_edge_correlation(res, res_known, df)
    # print(len(res_unknown))
    # print(len(res_known))
    # df = df_res
    # print(df.head())
    # print(similarity_map)
    if query_chk == 1:
        # start = time.time()
        try:
            similarity_map = molecular_sim(known, unknown, message)
            df["similarity"] = df.apply(
                lambda row: similarity_map[row["drug"]][0][1]
                if row["drug"] in similarity_map.keys()
                else np.nan,
                axis=1,
            )
        except Exception as e:
            logger.error(traceback.format_exc())
            df = df.assign(similarity=np.nan)

        # print(f"Time to compute Molecular Similarity:{time.time() - start}")
        # Step 3:
        # calculating the recency
        df["recency"] = df.apply(
            lambda row: recency_function_exp(
                row["number_of_publ"], row["age_oldest_pub"], 100, 50
            )
            if not (np.isnan(row["number_of_publ"]) or np.isnan(row["age_oldest_pub"]))
            else np.nan,
            axis=1,
        )
        #
        # # Step 4:
        # # Calculating the Similarity:
        # nearest_neighbours = calculate_nn_distance(res_known, res_unknown, 0, 1)

        # df = df.assign(similarity=np.nan)

        # # Step 5:
        # # Calculating the novelty score:
        df["novelty_score"] = df.apply(
            lambda row: novelty_score(
                row["fda status"], row["recency"], row["similarity"]
            ),
            axis=1,
        )
        # df.to_excel(f'DATAFRAME_result.xlsx', header=False, index=False)

        # # # Step 6
        # # # Just sort them:
        df = df[["drug", "novelty_score"]].sort_values(
            by="novelty_score", ascending=False
        )
    else:
        df = df.assign(novelty_score=0)
    # df.to_excel(f'DATAFRAME_NOVELTY.xlsx', header=False, index=False)
    return df
