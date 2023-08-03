#!/usr/bin/env python
import requests


def mol_to_smile_molpro(molecules):
    """
    
    Returns:
        Dict
    
    Args:
        List    
        Example: Attached sample_mols.json file 
    """
    
    url = "https://translator.broadinstitute.org/molecular_data_provider/compound/by_id"
    headers = {
        "accept": "application/json",
        "Content-Type": "application/json"
    }
    
    data_mol = list(set(molecules))
    smiles = {}

    response = requests.post(url, headers=headers, json=data_mol)

    if response.status_code == 200:
        json_response = response.json()
        #print(json_response)
        collec_url = json_response['url']
        temp_collec_response = requests.get(collec_url)
        if response.status_code == 200:
            collec_response = temp_collec_response.json()
            #print('\n')
            #print(collec_response['elements'][0]['identifiers'])
            if json_response['size'] > 0:
                
                for i in range(json_response['size']):
                    key_list = ['identifiers']
                    if set(key_list).issubset(collec_response['elements'][i].keys()):
                        #if 'smiles' in collec_response['elements'][i]['identifiers'].keys():
                        identifiers = collec_response['elements'][i]['identifiers']
                        smile = identifiers.get('smiles', 'No SMILES could be found')
                        smiles[data_mol[i]] = smile
                    else:
                        smiles[data_mol[i]] = 'No identifiers could be found'
                        
            #Recursion: re-attempt to retrieve the SMILES for those (from the initial data_mol)
            #for which the retrieval was not successful in the first attempt
            if len(list(smiles.keys())) < len(data_mol):
                diff = [mol for mol in data_mol if mol not in list(smiles.keys())]
                diff_smiles = mol_to_smile_molpro(diff)
                if len(diff_smiles) > 0:
                    smiles.update(diff_smiles)           
    else:
        print(f"Error: {response.status_code} - {response.text}")
        
    
    return smiles