#!/usr/bin/env python
# coding: utf-8

# MODIFIED BY: Darin Hauner and Neeraj Kumar
# Date: 3/10/2021
#
# added reaction
#     carbonyl

# MODIFIED BY: Darin Hauner and Neeraj Kumar
# Date: 9/22/2020
#
# Added Reactions:
#    alkynes
#    epoxyKetone_alphaSub_reaction
#    epoxyKetone_betaSub_reaction
#    esters
#    haloKetone_SN2_reaction
#
# Added Functionality
#    Reaction type included in resulting pdb file name
#    counter included in pdb file name to prevent duplicate names

# Python/RDkit script to convert SMILES files containing
# electrophilic groups to PDB files containing a covalent
# Cys linkage.# 
# Jerry M. Parks
# 7/7/2020
# 
# The input file should contain two columns: One SMILES string and 
# one molecule name/identifier per line. No header is needed.
# 
# Currently, the script works for the following reactive groups:
# 
# chloroacetamides
# acrylamides
# vinylmethyl esters
# vinyl sulfones
# 
# NOTES:
# 
# RDkit SMARTS reaction examples
# https://github.com/rdkit/rdkit-tutorials/blob/master/notebooks/003_SMARTS_ReactionsExamples.ipynb
# 
# PubChem Sketcher
# https://pubchem.ncbi.nlm.nih.gov/edit3/index.html
# 
# PDB modification:
# https://sourceforge.net/p/rdkit/mailman/message/36404236/
# 
# PDB "flavor":
# https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html

# In[1]:


from rdkit import rdBase
from rdkit import Chem
from rdkit.Chem import AllChem
#from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors

#from rdkit.Chem.Draw import IPythonConsole

import os.path
import sys
import tarfile
#import shutil

# for flattening tuples and lists
from itertools import chain

#IPythonConsole.ipython_useSVG=True


# In[2]:


# helper functions
def to_smiles(mol_tuple):
    return tuple(Chem.MolToSmiles(mol, isomericSmiles=True) for mol in mol_tuple)

def from_smiles(smiles_tuple):
    return tuple(Chem.MolFromSmiles(smiles) for smiles in smiles_tuple)

def mol_with_atom_index(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol

def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))


# Read in a list of SMILES from a file.
# In[3]:
path = os.getcwd()
supplier = Chem.SmilesMolSupplier('Fin_test.smi', titleLine=False)
#supplier = Chem.SmilesMolSupplier('resnick_2019_jacs_with_names.smi', titleLine=False)


# To do:
# De-salt (trifluoroacetate, etc.), remove arylnitro, bromo, amantidine compounds, etc.
# https://www.rdkit.org/docs/source/rdkit.Chem.SaltRemover.html

# In[4]:


#remover = SaltRemover(defnData="[Cl,Br]")


# In[5]:


mols = [x for x in supplier]

# Display the first 16 molecules.
#Draw.MolsToGridImage(mols, legends=[x.GetProp("_Name") for x in mols], molsPerRow = 4, maxMols = 16)
#mol_with_atom_index(suppl[0])


# In[6]:


# Cys side chain
thiol = Chem.MolFromSmiles('CCCS')
thiol_smarts = Chem.MolFromSmarts('[CH3:1][CH2:2][CH2:3][S:4]')

# Define reactions

Reactions = {
   'chloroacetamide': [thiol_smarts, AllChem.ReactionFromSmarts('[CH3:1][CH2:2][CH2:3][S:4].[C:5](-[C:6](=[O:7])-[N:8])-[Cl:9]>>[C:5](-[C:6](=[O:7])-[N:8])-[S:4]-[CH2:3]-[CH2:2]-[CH3:1]')],
   'acrylamide': [thiol_smarts, AllChem.ReactionFromSmarts('[CH3:1][CH2:2][CH2:3][S:4].[C:5](-[C:6](=[O:7])-[N:8])=[C:9]>>[C:5](-[C:6](=[O:7])-[N:8])-[C:9]-[S:4]-[CH2:3]-[CH2:2]-[CH3:1]')],
   'vinylsulfone': [thiol_smarts, AllChem.ReactionFromSmarts('[CH3:1][CH2:2][CH2:3][S:4].[S:5](=[O:6])(=[O:7])-[C:8]=[C:9]>>[S:5](=[O:6])(=[O:7])-[C:8]-[C:9]-[S:4]-[CH2:3]-[CH2:2]-[CH3:1]')],
   'vinylmethylester': [thiol_smarts, AllChem.ReactionFromSmarts('[CH3:1][CH2:2][CH2:3][S:4].[C:5]-[C:6]=[C:7]-[C:8](=[O:9])-[O:10]-[CH3:11]>>[C:5]-[C:6](-[C:7]-[C:8](=[O:9])-[O:10]-[CH3:11])-[S:4]-[CH2:3]-[CH2:2]-[CH3:1]')],
   'alkyne': [thiol_smarts, AllChem.ReactionFromSmarts('[CH3:1][CH2:2][CH2:3][S:4].[#6:5]-[C:6]#[C:7]>>[#6:5]-[C:6](-[S:4]-[CH2:3]-[CH2:2]-[CH3:1])=[CH:7]')],
   'alphaEpoxyKetone': [thiol_smarts, AllChem.ReactionFromSmarts('[CH3:1][CH2:2][CH2:3][S:4].[#6:5]-[#6:6](=[#8:7])-[#6:8]1-[#8:9]-[#6:10]-1-[#6:11]>>[#6:5]-[#6:6](=[#8:7])-[#6:8](-[S:4]-[CH2:3]-[CH2:2]-[CH3:1])-[#6:10](-[#8:9])-[#6:11]')],
   'betaEpoxyKetone': [thiol_smarts, AllChem.ReactionFromSmarts('[CH3:1][CH2:2][CH2:3][S:4].[#6:5]-[#6:6](=[#8:7])-[#6:8]1-[#8:9]-[#6:10]-1-[#6:11]>>[#6:5]-[#6:6](=[#8:7])-[#6:8](-[#8:9])-[#6:10](-[S:4]-[CH2:3]-[CH2:2]-[CH3:1])-[#6:11]')],
   'ester': [thiol_smarts, AllChem.ReactionFromSmarts('[CH3:1][CH2:2][CH2:3][S:4].[#6:5]-[C:6](=[O:7])-[O:8]-[#6:9]>>[#6:5]-[C:6](=[O:7])-[S:4]-[CH2:3]-[CH2:2]-[CH3:1]')],
   'haloKetone': [thiol_smarts, AllChem.ReactionFromSmarts('[CH3:1][CH2:2][CH2:3][S:4].[#6:5]-[C:6](=[O:7])-[#6:8]-[F,Cl,Br,I:9]>>[#6:5]-[C:6](=[O:7])-[#6:8]-[S:4]-[CH2:3]-[CH2:2]-[CH3:1]')],
   'carbonyl': [thiol_smarts, AllChem.ReactionFromSmarts('[CH3:1][CH2:2][CH2:3][S:4].[#6:5]=[O:6]>>[#6:5]([-O:6])-[S:4]-[CH2:3]-[CH2:2]-[CH3:1]')]
}


# Next, iterate over all molecules in mols and convert them to the corresponding thiol addition products.

# In[7]:

newmols = []
for R in Reactions:
   for mol in mols:
      try:
         if mol is None: continue
         D = Reactions[R] # get list holding MolFromSmarts
         RSmart = D[0] # Get Mol Object from Reactions
         RReact = D[1] # Get ChemicalReaction Object from Reactions

         molname = mol.GetProp("_Name")
         all_products_tuples = RReact.RunReactants((RSmart, mol))
         all_products = list(chain.from_iterable(all_products_tuples))
         all_products_smiles = [Chem.MolToSmiles(m, isomericSmiles=True) for m in all_products]
         all_products_unique = [Chem.MolFromSmiles(smiles) for smiles in set(all_products_smiles)]
         #mol.GetProp("_Name")
         #molname = mol.GetProp("_Name")
         [x.SetProp("_Name", molname + "_" + R + "_" + str(all_products_unique.index(x))) for x in all_products_unique]
         newmols.extend(all_products_unique)

      except:
         print("failed to convert" + molname )

#Draw.MolsToGridImage(newmols, legends=[x.GetProp("_Name") for x in newmols], molsPerRow = 4, maxMols = 16)


# Edit the PDB files

# In[8]:


if os.path.exists("pdb"):
    os.system("rm -rf pdb")

os.makedirs('pdb')
    
i = 0
for m in newmols:
    i = i + 1
    print(m.GetProp("_Name"))
    
    try:
        mh = Chem.AddHs(m)
        AllChem.EmbedMolecule(mh, randomSeed=123456)
        blk = Chem.MolToPDBBlock(mh)

        m2 = Chem.MolFromPDBBlock(blk, sanitize=False)

        Chem.SanitizeMol(m2, sanitizeOps=Chem.SANITIZE_ALL^Chem.SANITIZE_SETAROMATICITY)
        AllChem.EmbedMolecule(m2, useRandomCoords=True)
    
   
        AllChem.MMFFOptimizeMolecule(m2)
        
        # Get the indices of the Cys atoms.
        indices = m.GetSubstructMatches(thiol_smarts)
    
        # Convert tuple to list.
        sidechain_indices = list(chain(*indices)) 
        #print("sidechain_indices", sidechain_indices)

        # Edit the PDB file. First, assign all atoms to residue 'LIG', chain 'A' 
        # and residue number 999.
    
        atom_indices = list(range(0, m2.GetNumAtoms()))
    
        [m2.GetAtomWithIdx(j).GetMonomerInfo().SetResidueName('LIG') for j in atom_indices]
        [m2.GetAtomWithIdx(j).GetMonomerInfo().SetChainId('A') for j in atom_indices]
        [m2.GetAtomWithIdx(j).GetMonomerInfo().SetIsHeteroAtom(False) for j in atom_indices]
        [m2.GetAtomWithIdx(j).GetMonomerInfo().SetResidueNumber(999) for j in atom_indices]
    
        # Now fix the Cys atoms. Note that the residues are not in the correct order yet.
    
        [m2.GetAtomWithIdx(k).GetMonomerInfo().SetResidueName('CYS') for k in sidechain_indices]
        [m2.GetAtomWithIdx(k).GetMonomerInfo().SetResidueNumber(145) for k in sidechain_indices]
    
        # Find the index of the S atom in the Cys side chain and use it as an 
        # anchor to rename all the atoms in the side chain.
    
        for l in sidechain_indices:
            atomname = m2.GetAtomWithIdx(l).GetMonomerInfo().GetName()
            if 'S' in atomname:
                is_s = l
        
        if (is_s == sidechain_indices[3]):
            #print('SG is index', is_s)
            sg=sidechain_indices[3]
            cb=sidechain_indices[2]
            ca=sidechain_indices[1]
            c=sidechain_indices[0]
        elif (is_s == sidechain_indices[0]):
            #print('SG is index', is_s)
            sg=sidechain_indices[0]
            cb=sidechain_indices[1]
            ca=sidechain_indices[2]
            c=sidechain_indices[3]
        else:
            print("Problem with side chain indices")
            sys.exit()
        
        m2.GetAtomWithIdx(sg).GetMonomerInfo().SetName(' SG ')
        m2.GetAtomWithIdx(cb).GetMonomerInfo().SetName(' CB ')
        m2.GetAtomWithIdx(ca).GetMonomerInfo().SetName(' CA ')
        m2.GetAtomWithIdx(c).GetMonomerInfo().SetName(' C  ')
    
        pdbname = './pdb/' + m.GetProp("_Name") + '.pdb'

        Chem.MolToPDBFile(m2, pdbname, flavor=2)

    except:
        print("failed to convert", m.GetProp("_Name"))

#make_tarfile("pdb.tgz", "pdb")
#os.system("rm -rf pdb")

# In[9]:

#os.system("./pdb2pdbqt_jmp_0.1.sh")

