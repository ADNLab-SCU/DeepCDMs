# -*- coding: utf-8 -*-

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import os

path = "_"
write_dir = "_"
"""
Path to the sdf file of molecules to be filtered.
write_dir to the save_path of write out sdf file.
"""
filelist = os.listdir(path)

exclude_elements = ['H', 'C', 'O', 'N', 'P', 'S', 'Cl', 'F', 'Br', 'I', 'Si']

for i in filelist:
    """
    load file and read sdf.
    """
    file_path = path + '/' + i
    mols = []
    suppl = Chem.SDMolSupplier(file_path)
    for mol in suppl:
        has_excluded_atoms = any(atom.GetSymbol() not in exclude_elements for atom in mol.GetAtoms())
        if has_excluded_atoms:
            pass
        else:
            mols.append(mol)
    
    for mol in mols:
        exactmass = mol.GetProp("EXACT_MASS")
        p_smiles = mol.GetProp("PRECURSOR_SMILES")
        if float(exactmass) < 1000 and "." not in p_smiles:
            pass
        else:
            mols.remove(mol)

    write_path = write_dir + '/' + i
    write_op =  Chem.SDWriter(write_path)
    for mol in mols:
        write_op.write(mol)
    write_op.close()
    
