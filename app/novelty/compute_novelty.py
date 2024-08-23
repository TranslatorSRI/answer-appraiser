import pandas as pd
from datetime import date
import httpx
import numpy as np
import traceback
import json
from known import find_known_results
from extr_smile_molpro_by_id import mol_to_smile_molpro
from mol_similarity import find_nearest_neighbors
import os, sys
import time
import asyncio

# startup_utils_path = os.path.abspath("app/novelty/gene_nmf_adapter.py")
# startup_utils_dir = os.path.dirname(startup_utils_path)
# sys.path.insert(0, startup_utils_dir)
# print("added to python path the following directory: {}".format(startup_utils_dir))
# import gene_nmf_adapter as adapter
# import dcc.dcc_utils as dutils
# logger = dutils.get_logger(name=__name__)
# logger.info("added to pyton path the following directory: {}".format(startup_utils_dir))
# # print("added to python path the following directory: {}".format(startup_utils_dir))

"""
This script computes the novelty score for a list of results obtained for a 1-H response.
The steps for the ideal workflow are as follows:
1. Distinguish whether the result is a known result or an unknown result.
2. Compute FDA approval status to check if currently under approval or already approved.
3. Compute recency by noting the associated publications of the result.
4. Compute molecular similarity to identify if similar to existing drugs.
 
The end result of this script displays a table with values from different columns and accordingly lists the novelty score as well.
"""


async def molecular_sim(known, unknown, message, query_id):
    unknown_ids = []
    known_ids = []
    if len(unknown) > 0:
        for drug in unknown:
            s = list(message["results"][drug]["node_bindings"].keys())
            if message["results"][drug]['node_bindings'][s[0]][0]['id'] == query_id:
                unknown_ids.append(message["results"][drug]['node_bindings'][s[1]][0]['id'])
            else:
                unknown_ids.append(message["results"][drug]['node_bindings'][s[0]][0]['id'])
            # edge_attribute_sn = message["results"][drug]["node_bindings"][s[0]][0]["id"]
            # if (
            #     "PUBCHEM" in edge_attribute_sn
            #     or "CHEMBL" in edge_attribute_sn
            #     or "UNII" in edge_attribute_sn
            #     or "RXNORM" in edge_attribute_sn
            #     or "UMLS" in edge_attribute_sn
            #     or not "MONDO" in edge_attribute_sn
            # ):
            #     unknown_ids.append(edge_attribute_sn)
            # else:
            #     unknown_ids.append(
            #         message["results"][drug]["node_bindings"][s[1]][0]["id"]
            #     )

    if len(known) > 0:
        for drug in known:
            s = list(message["results"][drug]["node_bindings"].keys())
            if message["results"][drug]['node_bindings'][s[0]][0]['id'] == query_id:
                known_ids.append(message["results"][drug]['node_bindings'][s[1]][0]['id'])
            else:
                known_ids.append(message["results"][drug]['node_bindings'][s[0]][0]['id'])
            # edge_attribute_sn = message["results"][drug]["node_bindings"][s[0]][0]["id"]
            # if (
            #     "PUBCHEM" in edge_attribute_sn
            #     or "CHEMBL" in edge_attribute_sn
            #     or "UMLS" in edge_attribute_sn
            #     or "UNII" in edge_attribute_sn
            #     or "RXNORM" in edge_attribute_sn
            #     or not "MONDO" in edge_attribute_sn
            # ):
            #     known_ids.append(edge_attribute_sn)
            # else:
            #     known_ids.append(
            #         message["results"][drug]["node_bindings"][s[1]][0]["id"]
            #     )
    smile_unkown = await mol_to_smile_molpro(unknown_ids)
    smile_known = await mol_to_smile_molpro(known_ids)
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
        with httpx.AsyncClient(timeout=30) as client:
            response = client.get(full_url)
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
    Calculates the recency based on number of publications associated to each drug
    and age of the oldest publication

    Higher value of recency = 1 - (age of oldest publication * number of publications)

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
       oldest  = alter_age_of_oldest_publ
    elif np.isnan(age_of_oldest_publ) and not np.isnan(number_of_publ):
        oldest = alter_number_of_publ
    else:
        oldest = alter_number_of_publ * alter_age_of_oldest_publ

    recency = 1 - oldest

    return recency

def extracting_publications(message, result):
    """
    Function to extract publications information for every result.
    Four ARAs: ARAX, ARAGORN, BTE and UNSECRET produce results with publications informations
    """
    publications = []
    # print(result)
    for idi, i in enumerate(result['analyses']):
        result_resource = i['resource_id']

        if result_resource=="infores:unsecret-agent":
            for idj, j in enumerate(i['edge_bindings']['t_edge']):
                aux_graph, edges = [], []
                if 'creative' not in j['id']:
                    edges.append(j['id'])
                else:
                    for idl, l in enumerate(message['knowledge_graph']['edges'][j['id']]['attributes']):
                        if l['attribute_type_id'] == "biolink:support_graphs":
                            aux_graph.extend(l['value'])
                            break
                    for idl, l in enumerate(aux_graph):
                        edges.extend(message['auxiliary_graphs'][l]['edges'])

                for e in edges:
                    knowledge_graph_edge = message['knowledge_graph']['edges'][e]
                    for idl, l in enumerate(knowledge_graph_edge['attributes']):
                        if l['attribute_type_id'] == "biolink:publications":
                            publications.extend(l['value'])
                            break

        elif result_resource=="infores:aragorn":
            for idj, j in enumerate(i['edge_bindings']['t_edge']):
                for idl, l in enumerate(message['knowledge_graph']['edges'][j['id']]['attributes']):
                    if l['attribute_type_id'] == "biolink:publications":
                        publications.extend(l['value'])
                        break

        elif result_resource=="infores:arax":
            for idj, j in enumerate(i['support_graphs']):
                for idl, l in enumerate(message['auxiliary_graphs'][j]['edges']):
                    for idm, m in enumerate(message['knowledge_graph']['edges'][l]['attributes']):
                        if m['attribute_type_id'] == "biolink:publications":
                            publications.extend(l['value'])
                            break

        elif result_resource=="infores:biothings-explorer":
            for idj, j in enumerate(i['edge_bindings']['t_edge']):
                aux_graph, edges = [], []
                for idk, k in enumerate(message['knowledge_graph']['edges'][j['id']]['attributes']):
                    if k['attribute_type_id'] == "biolink:support_graphs":
                        aux_graph.extend(k['value'])

                for idk, k in enumerate(aux_graph):
                    edges.extend(message['auxiliary_graphs'][k]['edges'])


                for e in edges:
                    knowledge_graph_edge = message['knowledge_graph']['edges'][e]
                    for idl, l in enumerate(knowledge_graph_edge['attributes']):
                        if l['attribute_type_id'] == "biolink:publications":
                            publications.extend(l['value'])
                            break
    return publications


async def extracting_drug_fda_publ_date(message, unknown):
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

    query_known, query_unknown, query_chk = query_id(message)
    idi = -1
    for tmp in unknown:
        tmp_res = message["results"][tmp]["analyses"][0]["edge_bindings"]
        for tmp_1 in tmp_res:
            idi += 1
            edge_id = message["results"][tmp]["analyses"][0]["edge_bindings"][tmp_1][0][
                "id"
            ]
            edge = message["knowledge_graph"]["edges"][edge_id]
            if query_chk == 1:
                if (
                    "PUBCHEM" in edge["subject"]
                    or "CHEMBL" in edge["subject"]
                    or "UNII" in edge["subject"]
                    or "RXNORM" in edge["subject"]
                    or "UMLS" in edge["subject"]
                    or not "MONDO" in edge["subject"]
                ):
                    drug_idx = edge["subject"]
                else:
                    drug_idx = edge["object"]
                edge_attributes = edge.get("attributes") or []
                if len(edge_attributes) > 0:
                    att_type_id = {}
                    fda = []
                    pub = []
                    for i in range(len(edge_attributes)):
                        att_type_id[i] = edge_attributes[i]["attribute_type_id"]

                    for key in att_type_id.keys():
                        if att_type_id[key] in attribute_type_id_list_fda:
                            fda.append(key)
                        elif att_type_id[key] in attribute_type_id_list_pub:
                            pub.append(key)

                    if len(fda) > 0:
                        if edge_attributes[fda[0]]["value"] == "FDA Approval":
                            fda_status = 0.0
                        else:
                            fda_status = 1.0
                    else:
                        fda_status = None

                    # Publication
                    if len(pub) > 0:
                        publications = edge_attributes[pub[0]]["value"]
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
                                response_pub = await get_publication_info(
                                    publications_1
                                )
                                if response_pub["_meta"]["n_results"] == 0:
                                    age_oldest = np.nan
                                else:
                                    publ_year = []
                                    for key in response_pub["results"].keys():
                                        if "not_found" not in key:
                                            publ_year.extend(
                                                [
                                                    int(
                                                        response_pub["results"][key][
                                                            "pub_year"
                                                        ]
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
                    if "NCBI" in edge["subject"] or "GO" in edge["subject"]:
                        gene_idx = edge["subject"]
                    else:
                        gene_idx = edge["object"]
                    drug_idx_fda_status.append((idi, gene_idx))
                elif query_unknown in ["biolink:Disease", "biolink:Phenotype"]:
                    if "MONDO" in edge["subject"]:
                        dis_idx = edge["subject"]
                    else:
                        dis_idx = edge["object"]
                    drug_idx_fda_status.append((idi, dis_idx))
    if query_unknown == "biolink:SmallMolecule" or query_unknown == "biolink:ChemicalEntity" or query_unknown == "biolink:Drug":
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
    elif query_unknown == "biolink:Gene":
        DF = pd.DataFrame(drug_idx_fda_status, columns=["edge", "gene"])
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
        score = 0

    return score

def clinical_trials_kp(message):
    """
    For every result, check the clinical approval and clinical trials information in the result's nodes attributes
    in association with that disease
    """
    for i in message['query_graph']['nodes']:
        if 'ids' in message['query_graph']['nodes'][i].keys():
            query_id_node = message['query_graph']['nodes'][i]['ids'][0]
            break
    result_clinical_info = []
    for idi, i in enumerate(message['results']):
        node_binding_keys = list(i['node_bindings'].keys())
        if i['node_bindings'][node_binding_keys[0]][0]['id'] == query_id_node:
            result_id_node = i['node_bindings'][node_binding_keys[1]][0]['id']
        else:
            result_id_node = i['node_bindings'][node_binding_keys[0]][0]['id']
        if message['knowledge_graph']['nodes'][result_id_node]['attributes']==[]:
            result_clinical_info.append([idi, result_id_node, "No clinical information"])
        else:
            att_found, found = 0, 0
            for chk in message['knowledge_graph']['nodes'][result_id_node]['attributes']:
                if chk['attribute_type_id'] == "biothings_annotations":
                    att_found = 1
                    break
            if att_found==1:
                attribute_chk = chk['value'][0].keys()
                if 'clinical_approval' in attribute_chk:
                    clinical_approval_list = chk['value'][0]['clinical_approval']
                    for idj, j in enumerate(clinical_approval_list):
                        clinical_app_key = list(j['disease'].keys())[0]
                        if j['disease'][clinical_app_key] == query_id_node:
                            found=1
                            result_clinical_info.append([idi, result_id_node, j['status']])
                            break

                if found==0:
                    if 'clinical_trials' in attribute_chk:
                        clinical_trials_list = chk['value'][0]['clinical_trials']
                        phase_trials = []
                        for idj, j in enumerate(clinical_trials_list):
                            clinical_app_key = list(j[0]['disease'].keys())[0]
                            if j[0]['disease'][clinical_app_key]==query_id_node:
                                for k in j:
                                    phase_trials.append(k['phase'])
                        if len(phase_trials)>0:
                            if len(phase_trials) > 1:
                                phase_trials.sort(reverse=True)
                                result_clinical_info.append([idi, result_id_node, phase_trials[0]])
                            else:
                                result_clinical_info.append([idi, result_id_node, phase_trials])
                            found=1

                if found==0:
                    result_clinical_info.append([idi, result_id_node, "No clinical information"])

            else:
                result_clinical_info.append([idi, result_id_node, "No clinical information"])

    return result_clinical_info

def tdls(message):
    """
        For every result, check the TDLs information in the result's nodes attributes
    """
    for i in message['query_graph']['nodes']:
        if 'ids' in message['query_graph']['nodes'][i].keys():
            query_id_node = message['query_graph']['nodes'][i]['ids'][0]
            break
    result_tdl_info = []
    for idi, i in enumerate(message['results']):
        node_binding_keys = list(i['node_bindings'].keys())
        if i['node_bindings'][node_binding_keys[0]][0]['id'] == query_id_node:
            result_id_node = i['node_bindings'][node_binding_keys[1]][0]['id']
        else:
            result_id_node = i['node_bindings'][node_binding_keys[0]][0]['id']
        if message['knowledge_graph']['nodes'][result_id_node]['attributes'] == []:
            result_tdl_info.append([idi, result_id_node, "No TDLs information"])
        else:
            att_found, found = 0, 0
            for chk in message['knowledge_graph']['nodes'][result_id_node]['attributes']:
                if chk['attribute_type_id'] == "biothings_annotations":
                    att_found = 1
                    break
            if att_found == 1:
                attribute_chk = chk['value'][0].keys()
                if 'pharos' in attribute_chk:
                    tdl = chk['value'][0]['pharos']['tdl']
                    result_tdl_info.append([idi, result_id_node, tdl])
                    found=1
                if found == 0:
                    result_tdl_info.append([idi, result_id_node, "No TDLs information"])

            else:
                result_tdl_info.append([idi, result_id_node, "No TDLs information"])

    return result_tdl_info

async def compute_novelty(message, logger):
    """
    INPUT: JSON Response with merged annotated results for a 1-H query

    OUTPUT: Pandas DataFrame  with FDA Status, Recency, Similarity and Novelty score per result
    """

    today = date.today()

    query_keys = list(message['query_graph']['nodes'].keys())
    if 'ids' in message['query_graph']['nodes'][query_keys[0]].keys():
        query_id_node = message['query_graph']['nodes'][query_keys[0]]['ids'][0]
        result_id_cat = message['query_graph']['nodes'][query_keys[1]]['categories'][0]
    else:
        query_id_node = message['query_graph']['nodes'][query_keys[1]]['ids'][0]
        result_id_cat = message['query_graph']['nodes'][query_keys[0]]['categories'][0]

    df_numpy = []
    result_ids_list = []
    correct_results = []
    known_list = []
    unknown_list = []
    novelty_score = []
    for idi, i in enumerate(message['results']):
        curated=0
        node_binding_keys = list(i['node_bindings'].keys())
        if i['node_bindings'][node_binding_keys[0]][0]['id'] == query_id_node:
            result_id_node = i['node_bindings'][node_binding_keys[1]][0]['id']
        else:
            result_id_node = i['node_bindings'][node_binding_keys[0]][0]['id']
        df_numpy.append([query_id_node, result_id_node])
        result_node_cat = message['knowledge_graph']['nodes'][result_id_node]['categories'][0]
        if (result_node_cat in ["biolink:Protein", "biolink:Gene"] and result_id_cat in ["biolink:Protein", "biolink:Gene"]) or (result_node_cat in ["biolink:SmallMolecule", "biolink:Drug", "biolink:ChemicalEntity"] and result_id_cat in ["biolink:SmallMolecule", "biolink:Drug", "biolink:ChemicalEntity"]):
            result_ids_list.append(result_id_node)
            correct_results.append(idi)
            for idj, j in enumerate(i['analyses']):
                for idk, k in enumerate(j['edge_bindings']['t_edge']):
                    knowledge_graph_edge = message['knowledge_graph']['edges'][k['id']]
                    epc_found = 0
                    for idl, l in enumerate(knowledge_graph_edge['attributes']):
                        if l['attribute_type_id'] == 'biolink:knowledge_level':
                            epc_found = 1
                            if l['value'] != 'prediction':
                                curated = 1
                                df_numpy[idi].extend([l['attribute_type_id'], l['value']])
                                break
                    if curated == 1 and epc_found == 1:
                        break
                    elif curated==0 and epc_found==0:
                        for idl, l in enumerate(knowledge_graph_edge['sources']):
                            if l['resource_role'] == 'primary_knowledge_source':
                                if l['resource_id'] not in ['infores:arax', 'infores:aragorn', 'infores:biothings-explorer','infores:unsecret-agent','infores:improving', 'infores:cqs']:
                                    curated = 1
                                    df_numpy[idi].extend([l['resource_role'], l['resource_id']])
                                    break

                if curated==1:
                    break
            if curated==1:
                df_numpy[idi].extend(["Curated"])
                known_list.append(idi)
            else:
                publications = extracting_publications(message, i)
                df_numpy[idi].extend(["All KLs checked", "All Primary Knowledge Sources Checked"])
                df_numpy[idi].append(f"{publications}")
                unknown_list.append(idi)


        else:
            df_numpy[idi].extend(["Incorrect Category of Result", "Incorrect Category of Result", "Incorrect Category of Result"])


        if df_numpy[idi][4]==[] or df_numpy[idi][4]=='Curated' or df_numpy[idi][4]=='Incorrect Category of Result':
            df_numpy[idi].extend(["No Recency Score"])
            novelty_score.extend([0])
        else:
            number_of_publ = len(publications)
            publications_1 = ",".join(publications)
            try:
                response_pub = get_publication_info(publications_1)
                if response_pub["_meta"]["n_results"] == 0:
                    age_oldest = np.nan
                else:
                    publ_year = []
                    for key in response_pub["results"].keys():
                        if "not_found" not in key:
                            publ_year.extend([int(response_pub["results"][key]["pub_year"])])
                    age_oldest = today.year - min(publ_year)
            except ConnectionError as e:
                age_oldest = np.nan
            if not np.isnan(age_oldest):
                recency_val = recency_function_exp(number_of_publ, age_oldest, 100, 50)
                df_numpy[idi].extend([recency_val])
                novelty_score.extend([recency_val])
            else:
                df_numpy[idi].extend(["No Recency Score"])
                novelty_score.extend([0])

    if result_id_cat in ["biolink:Protein", "biolink:Gene"]:
        column_list = ['Query ID', 'Result ID', 'Location 1', 'Location 2', 'Curated or not ?', 'Recency Score',
                       'Gene Distinctiveness', 'TDLs', 'novelty_score']

        # map_result = adapter.get_gene_nmf_novelty_for_gene_list(list_input_genes=result_ids_list, log=True)
        # print(f"Gene Distinct: time.time() - start")
        # map_result_keys = list(map_result['gene_results'].keys())

        for idi, i in enumerate(message['results']):
            if idi in correct_results:
                node_binding_keys = list(i['node_bindings'].keys())
                if i['node_bindings'][node_binding_keys[0]][0]['id'] == query_id_node:
                    res = i['node_bindings'][node_binding_keys[1]][0]['id']
                else:
                    res = i['node_bindings'][node_binding_keys[0]][0]['id']
                # if res in map_result_keys:
                #     gene_distinct = map_result['gene_results'][res]['novelty_score']
                #     df_numpy[idi].extend([f"{gene_distinct}"])
                #     novelty_score[idi] = novelty_score[idi] + gene_distinct
                # else:
                df_numpy[idi].extend([f"No information"])

                if message['knowledge_graph']['nodes'][res]['attributes'] == []:
                    df_numpy[idi].extend([f"No TDLs Information"])
                else:
                    att_found, found = 0, 0
                    for chk in message['knowledge_graph']['nodes'][res]['attributes']:
                        if chk['attribute_type_id'] == "biothings_annotations":
                            att_found = 1
                            break
                    if att_found == 1:
                        attribute_chk = chk['value'][0].keys()
                        if 'pharos' in attribute_chk:
                            tdl = chk['value'][0]['pharos']['tdl']
                            df_numpy[idi].extend([f"{tdl}"])
                            if tdl == "Tdark":
                                if novelty_score[idi]==0 and df_numpy[idi][5] == "No Recency Score" and df_numpy[idi][4]!= "Curated":
                                    novelty_score[idi] = 0.8
                                elif novelty_score[idi] <= 0.4:
                                    novelty_score[idi] = 0.41
                            elif tdl == "Tclin":
                                if novelty_score[idi]==0 and df_numpy[idi][5] == "No Recency Score" and df_numpy[idi][4]!= "Curated":
                                    novelty_score[idi] = 0.2
                                elif novelty_score[idi] >0.5:
                                    novelty_score[idi] = 0.5
                            elif tdl=="Tbio" or tdl=="Tchem":
                                if novelty_score[idi] == 0 and df_numpy[idi][5] == "No Recency Score" and df_numpy[idi][4]!= "Curated":
                                    novelty_score[idi] = 0.5
                                elif novelty_score[idi] > 0.6:
                                    novelty_score[idi] = 0.6
                            found = 1
                        if found == 0:
                            df_numpy[idi].extend([f"No TDLs Information"])

                    else:
                        df_numpy[idi].extend([f"No TDLs Information"])


            else:
                df_numpy[idi].extend([f"Incorrect Category of Result"])
                df_numpy[idi].extend([f"Incorrect Category of Result"])


    elif result_id_cat in ["biolink:SmallMolecule", "biolink:Drug", "biolink:ChemicalEntity"]:
        column_list = ['Query ID', 'Result ID', 'Location 1', 'Location 2', 'Curated or not ?', 'Recency Score', 'Molecular Similarity', 'novelty_score']

        similarity_map = await molecular_sim(known_list, unknown_list, message, query_id_node)
        similarity_map_keys = list(similarity_map.keys())
        for idi, i in enumerate(message['results']):
            if idi in correct_results:
                node_binding_keys = list(i['node_bindings'].keys())
                if i['node_bindings'][node_binding_keys[0]][0]['id'] == query_id_node:
                    res = i['node_bindings'][node_binding_keys[1]][0]['id']
                else:
                    res = i['node_bindings'][node_binding_keys[0]][0]['id']

                if res in similarity_map_keys:
                    similarity = similarity_map[res][0][1]
                    chem_distinct = 1-similarity
                    df_numpy[idi].extend([f"{chem_distinct}"])
                    if novelty_score[idi] == 0:
                        novelty_score[idi] = chem_distinct
                    else:
                        novelty_score[idi] = novelty_score[idi] * (0.73 + chem_distinct)
                        if novelty_score[idi] > 1:
                            novelty_score[idi] = 1
                else:
                    df_numpy[idi].extend([f"{np.nan}"])
            else:
                df_numpy[idi].extend([f"Incorrect Category of Result"])

    for i in range(len(novelty_score)):
        df_numpy[i].append(novelty_score[i])


    df_numpy = pd.DataFrame(df_numpy, columns= column_list)
    # print(df_numpy.head())
    # df_numpy.to_csv("example_results.csv", index=False)
    return df_numpy
    # print(novelty_score)
    # known, unknown = find_known_results(message)
    #
    # '''
    # STEP 2: Extract FDA Status from knowledge edge
    # '''
    # df, query_chk = await extracting_drug_fda_publ_date(message, unknown)
    #
    # if result_id_cat == "biolink:SmallMolecule" or result_id_cat == "biolink:ChemicalEntity" or result_id_cat == "biolink:Drug":
    #     df["recency"] = df.apply(
    #         lambda row: (
    #             recency_function_exp(
    #                 row["number_of_publ"], row["age_oldest_pub"], 100, 50
    #             )
    #             if not (
    #                     np.isnan(row["number_of_publ"]) or np.isnan(row["age_oldest_pub"])
    #             )
    #             else np.nan
    #         ),
    #         axis=1,
    #     )
    #     try:
    #         similarity_map = await molecular_sim(known, unknown, message)
    #         df["similarity"] = df.apply(
    #             lambda row: (
    #                 similarity_map[row["drug"]][0][1]
    #                 if row["drug"] in similarity_map.keys()
    #                 else np.nan
    #             ),
    #             axis=1,
    #         )
    #     except Exception as e:
    #         logger.error(traceback.format_exc())
    #         df = df.assign(similarity=np.nan)
    #
    # '''
    # STEP 4b: Compute Gene Distinctiveness and TDLs
    #             Pending: Code integration issues for Gene Distinctiveness. Currently being discussed with Marc Duby
    # '''
    # if result_id_cat == "biolink:Gene":
    #     result_tdls_info = tdls(message)
    #     df["tdls"] = np.array(result_tdls_info)[:,2]
    #     result_genes = []
    #     for i in message['results']:
    #         result_keys = list(i['node_bindings'].keys())
    #         if i['node_bindings'][result_keys[0]][0]['id'] == query_id_node:
    #             if 'NCBIGene' in i['node_bindings'][result_keys[1]][0]['id']:
    #                 result_genes.append(i['node_bindings'][result_keys[1]][0]['id'])
    #         else:
    #             if 'NCBIGene' in i['node_bindings'][result_keys[0]][0]['id']:
    #                 result_genes.append(i['node_bindings'][result_keys[0]][0]['id'])
    #     print(df.head())
    #     # map_result = adapter.get_gene_nmf_novelty_for_gene_list(list_input_genes=result_genes, log=True)

        # print(result_genes)
        # print(map_result)

    # Step 1

    # #
    # # # Step 2
    #
    # # start = time.time()

    # # print(f"Time to extract fda status and Publication data:{time.time()-start}")
    # #         # print(df.head())
    # #         # print(query_chk)
    # #
    # # df.to_excel(f'DATAFRAME.xlsx', header=False, index=False)
    # # df = pd.read_excel('DATAFRAME.xlsx', names=['edge', 'drug', 'fda status', 'publications', 'number_of_publ', 'age_oldest_pub'])
    # # query_chk = 1
    #
    # # #
    # # res, res_known = extract_results(mergedAnnotatedOutput, unknown, known)
    # #         # print(len(res_known))
    # #         # print(len(res))
    # #         # # #
    # #         df_res, res_unknown, res_known = result_edge_correlation(res, res_known, df)
    # # print(len(res_unknown))
    # # print(len(res_known))
    # # df = df_res
    # # print(df.head())
    # # print(similarity_map)


        # print(f"Time to compute Molecular Similarity:{time.time() - start}")
        # Step 3:
        # calculating the recency

        #
        # # Step 4:
        # # Calculating the Similarity:
        # nearest_neighbours = calculate_nn_distance(res_known, res_unknown, 0, 1)

        # df = df.assign(similarity=np.nan)

        # # Step 5:
        # # Calculating the novelty score:
    #     df["novelty_score"] = df.apply(
    #         lambda row: novelty_score(
    #             row["fda status"], row["recency"], row["similarity"]
    #         ),
    #         axis=1,
    #     )
    #     # df.to_excel(f'DATAFRAME_result.xlsx', header=False, index=False)
    #
    #     # # # Step 6
    #     # # # Just sort them:
    #     df = df[["drug", "novelty_score"]].sort_values(
    #         by="novelty_score", ascending=False
    #     )
    # else:
    #     df = df.assign(novelty_score=0)
    # # df.to_excel(f'DATAFRAME_NOVELTY.xlsx', header=False, index=False)
    # return df

filename = f'example_GaucherDisease1.json'
message = json.load(open(filename))['fields']['data']['message']

asyncio.run(compute_novelty(message, None))