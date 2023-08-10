import httpx


def mol_to_smile_molpro(molecules):
    """
    Args:
        List

    Returns:
        Dict

    """

    url = "https://molepro.transltr.io/molecular_data_provider/compound/by_id"
    headers = {"accept": "application/json", "Content-Type": "application/json"}

    smiles = {}

    data_mol = list(set(molecules))
    # print(f'init data: {len(data_mol)}')
    with httpx.Client(timeout=30) as client:
        while data_mol:
            # print(f'before: {len(data_mol)}')
            data_mol_before = len(data_mol)
            response = client.post(url, headers=headers, json=data_mol)

            if response.status_code == 200:
                json_response = response.json()
                collec_url = json_response["url"]
                temp_collec_response = client.get(collec_url)
                if temp_collec_response.status_code == 200:
                    collec_response = temp_collec_response.json()

                    for i in range(json_response["size"]):
                        key_list = ["identifiers"]
                        if set(key_list).issubset(collec_response["elements"][i].keys()):
                            identifiers = collec_response["elements"][i]["identifiers"]
                            smile = identifiers.get("smiles", "No SMILES could be found")
                            smiles[data_mol[i]] = smile
                        else:
                            smiles[data_mol[i]] = "No identifiers could be found"

                    # Remove molecules with successfully retrieved smiles from data_mol
                    data_mol = [mol for mol in data_mol if mol not in smiles]
                    data_mol_after = len(data_mol)
                    # print(f'after: {len(data_mol)}')

                    if data_mol_after == data_mol_before:
                        break

                else:
                    print(
                        f"Error: {temp_collec_response.status_code} - {temp_collec_response.text}"
                    )
                    break
            else:
                print(f"Error: {response.status_code} - {response.text}")
                break

    return smiles
