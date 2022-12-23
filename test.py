import csv
import os
import Bio.PDB
import statistics
import json
from numpy import mean
import  Bio.PDB.Residue as Residue
path_positive_colab="pdbs/positive/colabfold/"
path_positive_omega="pdbs/positive/omegafold/"
path_positive_alphafold="pdbs/positive/alphafold/"
path_negative_colab="pdbs/negative/colabfold/"
path_negative_omega="pdbs/negative/omegafold/"
path_results="results/"
csv_heaver=['code', 'group','method','N interface residues ','pLDDT', 'IpLDDT']
csv_data=[]
path = '6ZE1_1b825_unrelaxed_rank_1_model_3.pdb' # your file path here
#path= "NEGATIVE_09_fcd4d_unrelaxed_rank_1_model_3.pdb" #negative
#path= "NEGATIVE_05_462ca_unrelaxed_rank_1_model_3.pdb"
#path= '6WX1_4897e_unrelaxed_rank_1_model_1.pdb'
#path='7FDO_5ce41_unrelaxed_rank_1_model_2.pdb'
json_file="6ZE1_1b825_unrelaxed_rank_1_model_3_scores.json"
p = Bio.PDB.PDBParser()
structure = p.get_structure('myStructureName', path)


ids = [a.get_id() for a in structure.get_atoms()]
bfactors= [a.get_bfactor() for a in structure.get_atoms()]
#print(bfactors)
avg=statistics.mean(bfactors)

print("mean of plDDT: ",avg)

f = open(json_file)
  
# returns JSON object as 
# a dictionary
data = json.load(f)
pae_list= [i for i in data['pae']]
# Iterating through the json
# list
#for i in data['pae']:
 #   pae_list.append(i)
# Closing file
f.close()
#print(pae_list)

print(type(pae_list))
print(len(pae_list))
avg=mean(pae_list)
print(avg)

from Bio.PDB import PDBParser

# create parser
parser = PDBParser()

# read structure from file
structure = parser.get_structure('id',path)

model = structure[0]
chain = model['B']
chain2 = model['C']
index_chain1=[]
index_chain2=[]
# this example uses only the first residue of a single chain.
# it is easy to extend this to multiple chains and residues.

"""
for residue1 in chain:
    for residue2 in chain2:
        if residue1 != residue2:
            # compute distance between CA atoms
            try:
                distance = residue1['CA'] - residue2['CA']
            except KeyError:
                ## no CA atom, e.g. for H_NAG
                continue
            if distance < 9:
                if residue1.get_id()[1] not in index_chain1:
                     index_chain1.append(residue1.get_id()[1])
                if residue2.get_id()[1] not in index_chain2:
                    index_chain2.append(residue2.get_id()[1])
                print(residue1, residue2, distance,residue1.get_atoms())
        # stop after first residue
        #break


print(index_chain1)
#print(bfactors)
#print(ids)





for residue in chain:
    for atom in residue.get_atoms():
        new_bfactors.append(atom.get_bfactor())

print(new_bfactors)
avg=mean(new_bfactors)
print(avg)

resid1=Residue_score(chain[0],2) 
print(resid1)


index_chain=0
new_bfactors=[]
residue_score=[]
for score in  bfactors:
    if score not in residue_score:
        residue_score

for residue1 in chain:
    for residue2 in chain2:
        if residue1 != residue2:
            try:
                distance = residue1['CA'] - residue2['CA']
            except KeyError:
                continue
            if distance < 8:
                residue_in_interface=True
                if residue1.get_id()[1] not in index_chain1:
                     index_chain1.append(residue1.get_id()[1])
                     for atom in residue1.get_atoms():
                        new_bfactors.append(atom.get_bfactor()*-1)
                        atom.set_bfactor(atom.get_bfactor()*-1)
                if residue2.get_id()[1] not in index_chain2:
                    index_chain2.append(residue2.get_id()[1])

                for atom in residue1.get_atoms():
                    new_bfactors.append(atom.get_bfactor()*-100)
                print(residue1, residue2, distance)
        # stop after first residue
        #break

bfactors2= [a.get_bfactor() for a in structure.get_atoms()]
"""

def get_plddt_interface(path):
    for file in os.scandir(path):
        if file.is_file():
            name=file.name.split(".")[0]
            if name.split("_")[0].lower()=="negative":
                protein_name=name.split("_")[0]+name.split("_")[1]
                #method="colabaf"
            else:
                protein_name=name.split("_")[0]
                method=path.split("/")[2]
            method=path.split("/")[2]
            protein_name=protein_name.upper()
            path_file=path+name+".pdb"
            
        residues=[]
        parser = PDBParser()
        structure = parser.get_structure('id',path_file)
        print(path_file)
        chains = structure[0]

        if (len(chains)>2):
            return
        chains_model=[]
        for chain in chains:
            chains_model.append(chain)

        chain = chains_model[0]
        chain2 =  chains_model[1]
        index_chain1=[]
        index_chain2=[]
        interface_residue1=[]
        interface_residue2=[]
        scores=[]
        interface_scores=[]
        for residue1 in chain:
                for residue2 in chain2:
                    if residue1 != residue2:
                        try:
                            distance = residue1['CA'] - residue2['CA']
                        except KeyError:
                            continue
                        if distance < 8:
                            if residue1.get_id()[1] not in index_chain1:
                                index_chain1.append(residue1.get_id()[1])
                                interface_residue1.append(residue1)
                            if residue2.get_id()[1] not in index_chain2:
                                index_chain2.append(residue2.get_id()[1])
                                interface_residue2.append(residue2)
                            
                            #print(residue1, residue2, distance)
        
        for residue in chain:
            for atom in residue.get_atoms():
                        scores.append(atom.get_bfactor())
        for residue in chain2:
            for atom in residue.get_atoms():
                        scores.append(atom.get_bfactor())

        if len(scores)<=0:
            pldtt=0
        else:
            pldtt=mean(scores)
        for residue in interface_residue1:
            for atom in residue.get_atoms():
                        interface_scores.append(atom.get_bfactor())
        

        for residue in interface_residue2:
            for atom in residue.get_atoms():
                        interface_scores.append(atom.get_bfactor())        
        if len(interface_scores)<=0:
            iplddt=0
        else:
            iplddt=mean(interface_scores)
        print("length",len(interface_residue1))
        n_residues= len(interface_residue1)+len(interface_residue2)
        

        csv_data.append([protein_name,path.split("/")[1],method,n_residues,pldtt,iplddt])
    
    return interface_residue2    

#print(get_plddt_interface(path_negative_colab))
print(get_plddt_interface(path_positive_colab))

with open('plddt_results_positive_colab.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(csv_heaver)
    writer.writerows(csv_data)