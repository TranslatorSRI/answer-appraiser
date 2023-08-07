#!/usr/bin/env python

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem


def find_nearest_neighbors(
    unknown_smiles_dict, known_smiles_dict, similarity_cutoff, num_neighbors
):
    """

    Returns:
        Dict

    Args:
        unknown_smiles_dict: Dict
        known_smiles_dict: Dict
        similarity_cutoff: float: 0
        num_neighbors: int: 1

    """
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
        if known_mol is None:
            raise ValueError("Invalid SMILES string for", key)
        else:
            known_mols.update({key: known_mol})
    nearest_neighbor_mapping = {}
    for unknownkey, value in unknown_smiles.items():
        query_mol = Chem.MolFromSmiles(value)
        if query_mol is None:
            raise ValueError("Invalid SMILES string")

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
    return nearest_neighbor_mapping
