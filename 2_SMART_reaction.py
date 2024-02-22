# -*- coding: utf-8 -*-

import os
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def amine_dansylation(mols):
    """amine dansylation"""
    product_amine = []
    inchikey_lst = []
    for mol in mols:
        reactions_1 = rxn_aromatic_amine.RunReactants((dnscl, mol))
        reactions_2 = rxn_aliphatic_amine.RunReactants((dnscl, mol))
        
        for tp in reactions_1:
            inchikey_product = Chem.MolToInchiKey(tp[0])
            if inchikey_product not in inchikey_lst:
                inchikey_lst.append(inchikey_product)
                product_smiles = Chem.MolToSmiles(tp[0])
                product_mol = Chem.MolFromSmiles(product_smiles)
                formula = str(rdMolDescriptors.CalcMolFormula(product_mol))
                exact_mass = str(rdMolDescriptors.CalcExactMolWt(product_mol))
                product_mol.SetProp("INCHIKEY", inchikey_product)
                product_mol.SetProp("FORMULA", formula)
                product_mol.SetProp("DANSYLATED_SMILES", product_smiles)
                product_mol.SetProp("EXACT_MASS", exact_mass)
                product_mol.SetProp("PRECURSOR_FORMULA", mol.GetProp("Mol_Formula"))
                product_mol.SetProp("PRECURSOR_SMILES", Chem.MolToSmiles(mol))
                product_mol.SetProp("PRECURSOR_EXACT_MASS", mol.GetProp("Monoisotopic_Mass"))
                product_mol.SetProp("PRECURSOR_ION_NAME", mol.GetProp("Preferred_name"))
                product_mol.SetProp("PRECURSOR_CASRN", mol.GetProp("CASRN"))
                product_amine.append(product_mol)
                
        for tp in reactions_2:
            inchikey_product = Chem.MolToInchiKey(tp[0])
            if inchikey_product not in inchikey_lst:
                inchikey_lst.append(inchikey_product)
                product_smiles = Chem.MolToSmiles(tp[0])
                product_mol = Chem.MolFromSmiles(product_smiles)
                formula = str(rdMolDescriptors.CalcMolFormula(product_mol))
                exact_mass = str(rdMolDescriptors.CalcExactMolWt(product_mol))
                product_mol.SetProp("INCHIKEY", inchikey_product)
                product_mol.SetProp("FORMULA", formula)
                product_mol.SetProp("DANSYLATED_SMILES", product_smiles)
                product_mol.SetProp("EXACT_MASS", exact_mass)
                product_mol.SetProp("PRECURSOR_FORMULA", mol.GetProp("Mol_Formula"))
                product_mol.SetProp("PRECURSOR_SMILES", Chem.MolToSmiles(mol))
                product_mol.SetProp("PRECURSOR_EXACT_MASS", mol.GetProp("Monoisotopic_Mass"))
                product_mol.SetProp("PRECURSOR_ION_NAME", mol.GetProp("Preferred_name"))
                product_mol.SetProp("PRECURSOR_CASRN", mol.GetProp("CASRN"))
                product_amine.append(product_mol)
    
    write_path = write_dir + '/amine_dansylation.sdf'
    write_op =  Chem.SDWriter(write_path)
    for mol in product_amine:
        write_op.write(mol)
    write_op.close()

def hydroxyl_dansylation(mols):
    """amine dansylation"""
    product_hydroxyl = []
    inchikey_lst = []
    for mol in mols:
        reactions_1 = rxn_aromatic_hydroxyl.RunReactants((dnscl, mol))
        reactions_2 = rxn_aliphatic_hydroxyl.RunReactants((dnscl, mol))
        
        for tp in reactions_1:
            inchikey_product = Chem.MolToInchiKey(tp[0])
            if inchikey_product not in inchikey_lst:
                inchikey_lst.append(inchikey_product)
                product_smiles = Chem.MolToSmiles(tp[0])
                product_mol = Chem.MolFromSmiles(product_smiles)
                formula = str(rdMolDescriptors.CalcMolFormula(product_mol))
                exact_mass = str(rdMolDescriptors.CalcExactMolWt(product_mol))
                product_mol.SetProp("INCHIKEY", inchikey_product)
                product_mol.SetProp("FORMULA", formula)
                product_mol.SetProp("DANSYLATED_SMILES", product_smiles)
                product_mol.SetProp("EXACT_MASS", exact_mass)
                product_mol.SetProp("PRECURSOR_FORMULA", mol.GetProp("Mol_Formula"))
                product_mol.SetProp("PRECURSOR_SMILES", Chem.MolToSmiles(mol))
                product_mol.SetProp("PRECURSOR_EXACT_MASS", mol.GetProp("Monoisotopic_Mass"))
                product_mol.SetProp("PRECURSOR_ION_NAME", mol.GetProp("Preferred_name"))
                product_mol.SetProp("PRECURSOR_CASRN", mol.GetProp("CASRN"))
                product_hydroxyl.append(product_mol)
                
        for tp in reactions_2:
            inchikey_product = Chem.MolToInchiKey(tp[0])
            if inchikey_product not in inchikey_lst:
                inchikey_lst.append(inchikey_product)
                product_smiles = Chem.MolToSmiles(tp[0])
                product_mol = Chem.MolFromSmiles(product_smiles)
                formula = str(rdMolDescriptors.CalcMolFormula(product_mol))
                exact_mass = str(rdMolDescriptors.CalcExactMolWt(product_mol))
                product_mol.SetProp("INCHIKEY", inchikey_product)
                product_mol.SetProp("FORMULA", formula)
                product_mol.SetProp("DANSYLATED_SMILES", product_smiles)
                product_mol.SetProp("EXACT_MASS", exact_mass)
                product_mol.SetProp("PRECURSOR_FORMULA", mol.GetProp("Mol_Formula"))
                product_mol.SetProp("PRECURSOR_SMILES", Chem.MolToSmiles(mol))
                product_mol.SetProp("PRECURSOR_EXACT_MASS", mol.GetProp("Monoisotopic_Mass"))
                product_mol.SetProp("PRECURSOR_ION_NAME", mol.GetProp("Preferred_name"))
                product_mol.SetProp("PRECURSOR_CASRN", mol.GetProp("CASRN"))
                product_hydroxyl.append(product_mol)
    
    write_path = write_dir + '/hydroxyl_dansylation.sdf'
    write_op =  Chem.SDWriter(write_path)
    for mol in product_hydroxyl:
        write_op.write(mol)
    write_op.close()

def carboxyl_dansylation(mols):
    """amine dansylation"""
    product_carboxyl = []
    inchikey_lst = []
    for mol in mols:
        reactions_1 = rxn_aromatic_carboxyl.RunReactants((dnscl, mol))
        reactions_2 = rxn_aliphatic_carboxyl.RunReactants((dnscl, mol))
        
        for tp in reactions_1:
            inchikey_product = Chem.MolToInchiKey(tp[0])
            if inchikey_product not in inchikey_lst:
                inchikey_lst.append(inchikey_product)
                product_smiles = Chem.MolToSmiles(tp[0])
                product_mol = Chem.MolFromSmiles(product_smiles)
                formula = str(rdMolDescriptors.CalcMolFormula(product_mol))
                exact_mass = str(rdMolDescriptors.CalcExactMolWt(product_mol))
                product_mol.SetProp("INCHIKEY", inchikey_product)
                product_mol.SetProp("FORMULA", formula)
                product_mol.SetProp("DANSYLATED_SMILES", product_smiles)
                product_mol.SetProp("EXACT_MASS", exact_mass)
                product_mol.SetProp("PRECURSOR_FORMULA", mol.GetProp("Mol_Formula"))
                product_mol.SetProp("PRECURSOR_SMILES", Chem.MolToSmiles(mol))
                product_mol.SetProp("PRECURSOR_EXACT_MASS", mol.GetProp("Monoisotopic_Mass"))
                product_mol.SetProp("PRECURSOR_ION_NAME", mol.GetProp("Preferred_name"))
                product_mol.SetProp("PRECURSOR_CASRN", mol.GetProp("CASRN"))
                product_carboxyl.append(product_mol)
                
        for tp in reactions_2:
            inchikey_product = Chem.MolToInchiKey(tp[0])
            if inchikey_product not in inchikey_lst:
                inchikey_lst.append(inchikey_product)
                product_smiles = Chem.MolToSmiles(tp[0])
                product_mol = Chem.MolFromSmiles(product_smiles)
                formula = str(rdMolDescriptors.CalcMolFormula(product_mol))
                exact_mass = str(rdMolDescriptors.CalcExactMolWt(product_mol))
                product_mol.SetProp("INCHIKEY", inchikey_product)
                product_mol.SetProp("FORMULA", formula)
                product_mol.SetProp("DANSYLATED_SMILES", product_smiles)
                product_mol.SetProp("EXACT_MASS", exact_mass)
                product_mol.SetProp("PRECURSOR_FORMULA", mol.GetProp("Mol_Formula"))
                product_mol.SetProp("PRECURSOR_SMILES", Chem.MolToSmiles(mol))
                product_mol.SetProp("PRECURSOR_EXACT_MASS", mol.GetProp("Monoisotopic_Mass"))
                product_mol.SetProp("PRECURSOR_ION_NAME", mol.GetProp("Preferred_name"))
                product_mol.SetProp("PRECURSOR_CASRN", mol.GetProp("CASRN"))
                product_carboxyl.append(product_mol)
    
    write_path = write_dir + '/carboxyl_dansylation.sdf'
    write_op =  Chem.SDWriter(write_path)
    for mol in product_carboxyl:
        write_op.write(mol)
    write_op.close()


file_path = "" # path to the dansylation sdf file
write_dir = "" # path to write out

"""load template of COOH & DnsHz reaction."""
template_carboxyl_aliphatic = '[N:1]-[NH2;D1;+0:2].[C:3]-[C;D3:4](=[O:5])-[OH:6]>>[N:1]-[NH;D2;+0:2]-[C;D3:4](-[C:3])=[O:5]'
template_carboxyl_aromatic = '[N:1]-[NH2;D1;+0:2].[c:3]-[C;D3:4](=[O:5])-[OH:6]>>[N:1]-[NH;D2;+0:2]-[C;D3:4](-[c:3])=[O:5]'
rxn_aliphatic_carboxyl = AllChem.ReactionFromSmarts(template_carboxyl_aliphatic)
rxn_aromatic_carboxyl = AllChem.ReactionFromSmarts(template_carboxyl_aromatic)
dnshz = Chem.MolFromSmiles('O=S(C1=CC=CC2=C(N(C)C)C=CC=C21)(NN)=O')

template_hydroxyl_aliphatic = '[Cl:1]-[S:2](=[O:3])(=[O:4])-[c:5].[OH:6][C;!$(C=O),!$(S=O),!$(P=O):7]>>[O:3]=[S:2](=[O:4])(-[c:5])-[O:6][C:7]'
template_hydroxyl_aromatic = '[Cl:1]-[S:2](=[O:3])(=[O:4])-[c:5].[OH:6][c;!$(C=O),!$(S=O),!$(P=O):7]>>[O:3]=[S:2](=[O:4])(-[c:5])-[O:6][c:7]'
rxn_aliphatic_hydroxyl = AllChem.ReactionFromSmarts(template_hydroxyl_aliphatic)
rxn_aromatic_hydroxyl = AllChem.ReactionFromSmarts(template_hydroxyl_aromatic)

template_amine_aromatic = '[Cl:1]-[S:2](=[O:3])(=[O:4])-[c:5].[N!H0;!$([N]C(=O)):6][c:7]>>[O:3]=[S:2](=[O:4])(-[c:5])-[N:6][c:7]'
template_amine_aliphatic = '[Cl:1]-[S:2](=[O:3])(=[O:4])-[c:5].[N!H0;!$([N]C(=O)):6][C:7]>>[O:3]=[S:2](=[O:4])(-[c:5])-[N:6][C:7]'
rxn_aliphatic_amine = AllChem.ReactionFromSmarts(template_amine_aliphatic)
rxn_aromatic_amine = AllChem.ReactionFromSmarts(template_amine_aromatic)
dnscl = Chem.MolFromSmiles('CN(C)C1=C(C=CC=C2S(=O)(Cl)=O)C2=CC=C1')

mols = []
num_failed = 0
max_atoms = 100

suppl = Chem.SDMolSupplier(file_path)
for idx, mol in enumerate(suppl):
    if mol is not None and (mol.GetNumAtoms() <= max_atoms or max_atoms is None):
        mols.append(mol)
    else:
        num_failed += 1

amine_dansylation(mols)
hydroxyl_dansylation(mols)
carboxyl_dansylation(mols)