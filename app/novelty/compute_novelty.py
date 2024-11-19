import pandas as pd
from datetime import date
import httpx
import numpy as np
import traceback
import json
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
import os, sys
import time
import asyncio
import requests
import redis

startup_utils_path = os.path.abspath("app/novelty/gene_nmf_adapter.py")
startup_utils_dir = os.path.dirname(startup_utils_path)
sys.path.insert(0, startup_utils_dir)
print("added to python path the following directory: {}".format(startup_utils_dir))
import gene_nmf_adapter as adapter
import dcc.dcc_utils as dutils

"""
This script computes the novelty score for a list of results obtained for a 1-H response.
The steps for the ideal workflow are as follows:
1. Distinguish whether the result is a known result or an unknown result.
2. Compute FDA approval status to check if currently under approval or already approved.
3. Compute recency by noting the associated publications of the result.
4. Compute molecular similarity to identify if similar to existing drugs.
 
The end result of this script displays a table with values from different columns and accordingly lists the novelty score as well.
"""

def find_nearest_neighbors( unknown_smiles_dict, known_smiles_dict, similarity_cutoff, num_neighbors):
    """

    Returns:
        Dict

    Args:
        unknown_smiles_dict: Dict
        known_smiles_dict: Dict
        similarity_cutoff: float: 0
        num_neighbors: int: 1

    """
    start = time.time()
    unknown_smiles = {
        key: value
        for key, value in unknown_smiles_dict.items()
        if value != "No SMILES could be found"
    }
    known_smiles = {
        key: value
        for key, value in known_smiles_dict.items()
        if value != "No SMILES could be found"
    }

    # Convert input SMILES to a molecule
    known_mols = {}
    for key, value in known_smiles.items():
        known_mol = Chem.MolFromSmiles(value)
        if known_mol is not None:
            known_mols.update({key: known_mol})
        # else:
        #     known_mols.update({key: known_mol})
    nearest_neighbor_mapping = {}
    for unknownkey, value in unknown_smiles.items():
        query_mol = Chem.MolFromSmiles(value)
        if query_mol is None:
            neighbors = []
            neighbors.append((0, 1))
            nearest_neighbor_mapping.update({unknownkey: neighbors})
        else:
            # Calculate fingerprints for the query molecule
            query_fp = AllChem.GetMorganFingerprint(query_mol, 2)

            # Calculate similarity scores between the query molecule and all molecules in the dataset
            similarities = []
            for key, mol in known_mols.items():
                fp = AllChem.GetMorganFingerprint(mol, 2)
                similarity = DataStructs.TanimotoSimilarity(query_fp, fp)
                similarities.append((key, similarity))

            # Sort the similarities in descending order
            similarities.sort(key=lambda x: x[1], reverse=True)

            # Retrieve the nearest neighbors
            neighbors = []
            for i in range(min(num_neighbors, len(similarities))):
                index, similarity = similarities[i]
                if similarity >= similarity_cutoff:
                    neighbors.append((index, similarity))
            nearest_neighbor_mapping.update({unknownkey: neighbors})
    end = time.time()
    print(f"find_nearest_neighbors function time: {end-start}")
    return nearest_neighbor_mapping

async def mol_to_smile_molpro(molecules):
    start=time.time()
    url = "https://molepro-db-transformer.transltr.io/moleprodb/elements/transform"
    smiles = {}
    data_mol = list(set(molecules))
    message = {
        "controls":
            [
            ]
    }
    try:
        for i in data_mol:
            message["controls"].append({
                "name": "id",
                "value": i
            })
        response = requests.post(url, json=message)
        if response.status_code == 200:
            for idi, i in enumerate(response.json()):
                if "smiles" in i['identifiers'].keys():
                    if i['identifiers']['smiles'] != None:
                        smiles[data_mol[idi]] = i['identifiers']['smiles']
                    else:
                        smiles[data_mol[idi]] = "No SMILES could be found"
                else:
                    smiles[data_mol[idi]] = "No SMILES could be found"

        else:
            print(f"Error: {response.status_code} - {response.text}")

    except Exception as e:
        for idx, i in enumerate(data_mol):
            smiles[i] = "No SMILES could be found"
    end = time.time()
    print(f"smiles computation time:{end-start}")
    # print(smiles)
    return smiles

# async def mol_to_smile_molpro(molecules):
#     start=time.time()
#     """
#     Args:
#         List
#
#     Returns:
#         Dict
#
#     """
#
#     url = "https://molepro.transltr.io/molecular_data_provider/compound/by_id"
#     headers = {"accept": "application/json", "Content-Type": "application/json"}
#
#     smiles = {}
#
#     data_mol = list(set(molecules))
#     # async with httpx.AsyncClient(timeout=30) as client:
#     try:
#         while data_mol:
#             data_mol_before = len(data_mol)
#             response = requests.post(url, headers=headers, json=data_mol)
#             if response.status_code == 200:
#                 json_response = response.json()
#                 collec_url = json_response["url"]
#                 temp_collec_response = requests.get(collec_url, timeout=1)
#                 if temp_collec_response.status_code == 200:
#                     collec_response = temp_collec_response.json()
#
#                     for i in range(json_response["size"]):
#                         key_list = ["identifiers"]
#                         if set(key_list).issubset(collec_response["elements"][i].keys()):
#                             identifiers = collec_response["elements"][i]["identifiers"]
#                             smile = identifiers.get(
#                                 "smiles", "No SMILES could be found"
#                             )
#                             smiles[data_mol[i]] = smile
#                         else:
#                             smiles[data_mol[i]] = "No identifiers could be found"
#
#                     # Remove molecules with successfully retrieved smiles from data_mol
#                     data_mol = [mol for mol in data_mol if mol not in smiles]
#                     data_mol_after = len(data_mol)
#                     # print(f'after: {len(data_mol)}')
#
#                     if data_mol_after == data_mol_before:
#                         break
#
#                 else:
#                     print(
#                         f"Error: {temp_collec_response.status_code} - {temp_collec_response.text}"
#                     )
#                     break
#             else:
#                 print(f"Error: {response.status_code} - {response.text}")
#                 break
#     except Exception as e:
#         for idx, i in enumerate(data_mol):
#             smiles[i] = "No identifiers could be found"
#     end = time.time()
#     print(f"smiles computation time:{end-start}")
#     return smiles



async def molecular_sim(known, unknown, message, query_id):
    start = time.time()
    unknown_ids = []
    known_ids = []
    if len(unknown) > 0:
        for drug in unknown:
            s = list(message["results"][drug]["node_bindings"].keys())
            if message["results"][drug]['node_bindings'][s[0]][0]['id'] == query_id:
                unknown_ids.append(message["results"][drug]['node_bindings'][s[1]][0]['id'])
            else:
                unknown_ids.append(message["results"][drug]['node_bindings'][s[0]][0]['id'])

    if len(known) > 0:
        for drug in known:
            s = list(message["results"][drug]["node_bindings"].keys())
            if message["results"][drug]['node_bindings'][s[0]][0]['id'] == query_id:
                known_ids.append(message["results"][drug]['node_bindings'][s[1]][0]['id'])
            else:
                known_ids.append(message["results"][drug]['node_bindings'][s[0]][0]['id'])

    smile_unkown = await mol_to_smile_molpro(unknown_ids)
    smile_known = await mol_to_smile_molpro(known_ids)
    similarity_map = find_nearest_neighbors(smile_unkown, smile_known, 0, 1)
    end = time.time()
    print(f"molecular_sim function time:{end-start}")
    return similarity_map

def get_publication_info(pub_id):
    """
    Args: PMID
    Returns: The publication info
    """

    # Connect to Redis
    r = redis.Redis(host='localhost', port=6379, db=0)
    pmid_years = []
    for key in pub_id:
        val = r.get(key)
        if val:
            pmid_years.append(int(val))
    # print(pmid_years)
    return pmid_years

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
    for idi, i in enumerate(result['analyses']):
        result_resource = i['resource_id']
        edge_keys = list(i['edge_bindings'].keys())
        if result_resource=="infores:unsecret-agent":
            for idj, j in enumerate(i['edge_bindings'][edge_keys[0]]):
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
            for idj, j in enumerate(i['edge_bindings'][edge_keys[0]]):
                for idl, l in enumerate(message['knowledge_graph']['edges'][j['id']]['attributes']):
                    if l['attribute_type_id'] == "biolink:publications":
                        publications.extend(l['value'])
                        break

        elif result_resource=="infores:arax":
            for idj, j in enumerate(i['support_graphs']):
                for idl, l in enumerate(message['auxiliary_graphs'][j]['edges']):
                    for idm, m in enumerate(message['knowledge_graph']['edges'][l]['attributes']):
                        if m['attribute_type_id'] == "biolink:publications":
                            publications.extend(m['value'])
                            break

        elif result_resource=="infores:biothings-explorer":
            for idj, j in enumerate(i['edge_bindings'][edge_keys[0]]):
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

async def compute_novelty(message, logger, wt_rec_tdl = 0.3, wt_gd = 0.7, wt_rec_clin = 0.3, wt_md = 0.7):
    """
    INPUT: JSON Response with merged annotated results for a 1-H query

    OUTPUT: Pandas DataFrame  with FDA Status, Recency, Similarity and Novelty score per result
    """
    start = time.time()
    today = date.today()

    query_keys = list(message['query_graph']['nodes'].keys())
    if 'ids' in message['query_graph']['nodes'][query_keys[0]].keys():
        if message['query_graph']['nodes'][query_keys[0]]['ids'] != None:
            query_id_node = message['query_graph']['nodes'][query_keys[0]]['ids'][0]
            result_id_cat = message['query_graph']['nodes'][query_keys[1]]['categories'][0]
        else:
            query_id_node = message['query_graph']['nodes'][query_keys[1]]['ids'][0]
            result_id_cat = message['query_graph']['nodes'][query_keys[0]]['categories'][0]
    else:
        query_id_node = message['query_graph']['nodes'][query_keys[1]]['ids'][0]
        result_id_cat = message['query_graph']['nodes'][query_keys[0]]['categories'][0]
    try:
        df_numpy = []
        result_ids_list = []
        correct_results = []
        known_list = []
        unknown_list = []
        novelty_score_rec = []
        novelty_score_gd, novelty_score_md = [], []
        novelty_score_rec_tdl, novelty_score_rec_clin = [], []
        for idi, i in enumerate(message['results']):
            curated=0
            node_binding_keys = list(i['node_bindings'].keys())
            if i['node_bindings'][node_binding_keys[0]][0]['id'] == query_id_node:
                result_id_node = i['node_bindings'][node_binding_keys[1]][0]['id']
            else:
                result_id_node = i['node_bindings'][node_binding_keys[0]][0]['id']
            df_numpy.append([query_id_node, result_id_node])
            result_node_cat = message['knowledge_graph']['nodes'][result_id_node]['categories'][0]
            if (result_node_cat in ["biolink:Protein", "biolink:Gene"] and result_id_cat in ["biolink:Protein", "biolink:Gene"]) or (result_node_cat in ["biolink:SmallMolecule", "biolink:Drug", "biolink:ChemicalEntity", "biolink:MolecularMixture", "biolink:MolecularEntity"] and result_id_cat in ["biolink:SmallMolecule", "biolink:Drug", "biolink:ChemicalEntity"]):
                result_ids_list.append(result_id_node)
                correct_results.append(idi)
                for idj, j in enumerate(i['analyses']):
                    edge_keys = list(j['edge_bindings'].keys())
                    for idk, k in enumerate(j['edge_bindings'][edge_keys[0]]):
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
                                    if l['resource_id'] not in ['infores:arax', 'infores:aragorn', 'infores:biothings-explorer','infores:unsecret-agent','infores:improving-agent', 'infores:cqs']:
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
                novelty_score_rec.append(0)
            else:
                number_of_publ = len(publications)
                if number_of_publ!=0:
                    # publications_1 = ",".join(publications)
                    try:
                        response_pub = get_publication_info(
                            publications
                        )
                        if not response_pub:
                            age_oldest = np.nan
                        else:
                            publ_year = response_pub
                            age_oldest = today.year - min(publ_year)
                    except ConnectionError as e:
                        age_oldest = np.nan
                    if not np.isnan(age_oldest):
                        recency_val = recency_function_exp(number_of_publ, age_oldest, 100, 50)
                        df_numpy[idi].extend([recency_val])
                        novelty_score_rec.append(recency_val)
                    else:
                        df_numpy[idi].extend(["No Recency Score"])
                        novelty_score_rec.append(0)
                else:
                    df_numpy[idi].extend(["No Recency Score"])
                    novelty_score_rec.append(0)
        recency_done = time.time()
        print(f"Recency Time: {recency_done-start}")
        if result_id_cat in ["biolink:Protein", "biolink:Gene"]:
            column_list = ['Query ID', 'Result ID', 'Location 1', 'Location 2', 'Curated or not ?', 'Recency Score',
                           'Gene Distinctiveness', 'TDLs', 'novelty_score']

            map_result = adapter.get_gene_nmf_novelty_for_gene_list(list_input_genes=result_ids_list, log=True)
            map_result_keys = list(map_result['gene_results'].keys())
            for idi, i in enumerate(message['results']):
                if idi in unknown_list:
                    node_binding_keys = list(i['node_bindings'].keys())
                    if i['node_bindings'][node_binding_keys[0]][0]['id'] == query_id_node:
                        res = i['node_bindings'][node_binding_keys[1]][0]['id']
                    else:
                        res = i['node_bindings'][node_binding_keys[0]][0]['id']
                    if res in map_result_keys:
                        gene_distinct = 1 - map_result['gene_results'][res]['novelty_score']
                        df_numpy[idi].extend([f"{gene_distinct}"])
                        novelty_score_gd.append(gene_distinct)
                    else:
                        df_numpy[idi].extend([f"No Gene Distinctiveness information"])
                        novelty_score_gd.append(0)

                    if message['knowledge_graph']['nodes'][res]['attributes'] == []:
                        df_numpy[idi].extend([f"No TDLs Information"])
                        novelty_score_rec_tdl.append(novelty_score_rec[idi])
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
                                    if novelty_score_rec[idi]==0 and df_numpy[idi][5] == "No Recency Score" and df_numpy[idi][4]!= "Curated":
                                        novelty_score_rec_tdl.append(0.8)
                                    elif novelty_score_rec[idi] <= 0.4:
                                        novelty_score_rec_tdl.append(0.41)
                                    else:
                                        novelty_score_rec_tdl.append(novelty_score_rec[idi])
                                elif tdl == "Tclin":
                                    if novelty_score_rec[idi]==0 and df_numpy[idi][5] == "No Recency Score" and df_numpy[idi][4]!= "Curated":
                                        novelty_score_rec_tdl.append(0.2)
                                    elif novelty_score_rec[idi] > 0.5:
                                        novelty_score_rec_tdl.append(0.5)
                                    else:
                                        novelty_score_rec_tdl.append(novelty_score_rec[idi])
                                elif tdl=="Tbio" or tdl=="Tchem":
                                    if novelty_score_rec[idi] == 0 and df_numpy[idi][5] == "No Recency Score" and df_numpy[idi][4]!= "Curated":
                                        novelty_score_rec_tdl.append(0.5)
                                    elif novelty_score_rec[idi] > 0.6:
                                        novelty_score_rec_tdl.append(0.6)
                                    else:
                                        novelty_score_rec_tdl.append(novelty_score_rec[idi])
                                found = 1
                            if found == 0:
                                df_numpy[idi].extend([f"No TDLs Information"])
                                novelty_score_rec_tdl.append(novelty_score_rec[idi])

                        else:
                            df_numpy[idi].extend([f"No TDLs Information"])
                            novelty_score_rec_tdl.append(novelty_score_rec[idi])


                else:
                    df_numpy[idi].extend([f"Incorrect Category of Result"])
                    df_numpy[idi].extend([f"Incorrect Category of Result"])
                    novelty_score_gd.append(0)
                    novelty_score_rec_tdl.append(0)

            for c in range(len(novelty_score_rec)):
                novelty_score_final = wt_rec_tdl * novelty_score_rec_tdl[c] + wt_gd * novelty_score_gd[c]
                df_numpy[c].append(novelty_score_final)

        elif result_id_cat in ["biolink:SmallMolecule", "biolink:Drug", "biolink:ChemicalEntity", "biolink:MolecularMixture", "biolink:MolecularEntity"]:
            column_list = ['Query ID', 'Result ID', 'Location 1', 'Location 2', 'Curated or not ?', 'Recency Score', 'Molecular Similarity', 'Clinical Information', 'novelty_score']

            similarity_map = await molecular_sim(known_list, unknown_list, message, query_id_node)
            similarity_map_keys = list(similarity_map.keys())
            for idi, i in enumerate(message['results']):
                if idi in unknown_list:
                    node_binding_keys = list(i['node_bindings'].keys())
                    if i['node_bindings'][node_binding_keys[0]][0]['id'] == query_id_node:
                        res = i['node_bindings'][node_binding_keys[1]][0]['id']
                    else:
                        res = i['node_bindings'][node_binding_keys[0]][0]['id']

                    if res in similarity_map_keys and similarity_map[res]!=[]:
                        similarity = similarity_map[res][0][1]
                        chem_distinct = 1-similarity
                        df_numpy[idi].extend([f"{chem_distinct}"])
                        novelty_score_md.append(chem_distinct)
                    else:
                        df_numpy[idi].extend([f"{np.nan}"])
                        novelty_score_md.append(0)

                    if message['knowledge_graph']['nodes'][res]['attributes'] == []:
                        df_numpy[idi].extend([f"No clinical Information"])
                        novelty_score_rec_clin.append(novelty_score_rec[idi])
                    else:
                        att_found, found = 0, 0
                        for chk in message['knowledge_graph']['nodes'][res]['attributes']:
                            if chk['attribute_type_id'] == "biothings_annotations":
                                att_found = 1
                                break
                        if att_found == 1:
                            attribute_chk = chk['value'][0].keys()
                            if 'clinical_approval' in attribute_chk:
                                clinical_approval_list = chk['value'][0]['clinical_approval']
                                for idj, j in enumerate(clinical_approval_list):
                                    clinical_app_key = list(j['disease'].keys())[0]
                                    if j['disease'][clinical_app_key] == query_id_node:
                                        found = 1
                                        df_numpy[idi].extend([j['status']])
                                        novelty_score_rec_clin.append(0)
                                        novelty_score_md[idi] = 0
                                        break

                            if found == 0:
                                if 'clinical_trials' in attribute_chk:
                                    clinical_trials_list = chk['value'][0]['clinical_trials']
                                    phase_trials = []
                                    for idj, j in enumerate(clinical_trials_list):
                                        clinical_app_key = list(j[0]['disease'].keys())[0]
                                        if j[0]['disease'][clinical_app_key] == query_id_node:
                                            for k in j:
                                                if "not_provided" not in k['phase']:
                                                    phase_trials.append(k['phase'].lstrip("clinical_trial_phase"))
                                    if len(phase_trials) > 0:
                                        phase_trials.sort(reverse=True)
                                        df_numpy[idi].extend([phase_trials[0]])
                                        if "4" in phase_trials[0]:
                                            novelty_score_rec_clin.append(0)
                                            novelty_score_md[idi] = 0
                                        elif "3" in phase_trials[0]:
                                            novelty_score_rec_clin.append(novelty_score_rec[idi]*0.2)
                                        elif "2" in phase_trials[0]:
                                            novelty_score_rec_clin.append(novelty_score_rec[idi]*0.6)
                                        elif "1" in phase_trials[0]:
                                            novelty_score_rec_clin.append(novelty_score_rec[idi]*0.8)
                                        else:
                                            novelty_score_rec_clin.append(novelty_score_rec[idi])
                                        found = 1
                                        # print(idi, phase_trials)
                            if found == 0:
                                df_numpy[idi].extend([f"No clinical information"])
                                novelty_score_rec_clin.append(novelty_score_rec[idi])

                        else:
                            df_numpy[idi].extend([f"No clinical information"])
                            novelty_score_rec_clin.append(novelty_score_rec[idi])

                else:
                    df_numpy[idi].extend([f"Incorrect Category of Result"])
                    df_numpy[idi].extend([f"Incorrect Category of Result"])
                    novelty_score_md.append(0)
                    novelty_score_rec_clin.append(0)
            for c in range(len(novelty_score_rec)):
                novelty_score_final = wt_rec_clin * novelty_score_rec_clin[c] + wt_md * novelty_score_md[c]
                df_numpy[c].append(novelty_score_final)
        print(f"Similarity Computation Time: {time.time()-recency_done}")
        df_numpy = pd.DataFrame(df_numpy, columns= column_list)

    except Exception as e:
        print(f"Error Encountered: {e}")
        if result_id_cat in ["biolink:SmallMolecule", "biolink:Drug", "biolink:ChemicalEntity", "biolink:MolecularMixture", "biolink:MolecularEntity"]:
            column_list = ['Query ID', 'Result ID', 'Location 1', 'Location 2', 'Curated or not ?', 'Recency Score', 'Molecular Similarity', 'Clinical Information', 'novelty_score']
        elif result_id_cat in ["biolink:Protein", "biolink:Gene"]:
            column_list = ['Query ID', 'Result ID', 'Location 1', 'Location 2', 'Curated or not ?', 'Recency Score',
                           'Gene Distinctiveness', 'TDLs', 'novelty_score']
        df_numpy = pd.DataFrame([[0]*len(column_list)]*len(message['results']), columns = column_list)

    return df_numpy

# filename = 'e2ed43f8-b2c7-4f0e-9eb8-169bdca50444.json'
# filename = "7aa4575e-8f98-47fe-931a-85bd6178ed46.json"
# # filename = '2e2cb6f5-10f7-492e-adb8-db37ea99ce31.json'
# message = json.load(open(filename))['fields']['data']['message']
# print(json.load(open(filename))['fields']['status'])
# if json.load(open(filename))['fields']['status'] == "Done":
#    df = asyncio.run(compute_novelty(message, None))
#    df.to_csv(f"{filename.rstrip('.json')}_novelty_results.csv", index=False)
# else:
#    print(json.load(open(filename))['fields']['status'])