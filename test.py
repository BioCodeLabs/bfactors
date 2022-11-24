import Bio.PDB
import statistics
import json
from numpy import mean
import  Bio.PDB.Residue as Residue
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

"""
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

#print(new_bfactors)
#print(bfactors)
#print(len(bfactors))
#print(len(new_bfactors))
#print(len(index_chain1))
#print(index_chain1)


print("average plDDT before: ",mean(bfactors))

bfactors2= [a.get_bfactor() for a in structure.get_atoms()]
bfactor_test=[]
#print(bfactors2)
for factor in bfactors2:
    if factor<0:
        factor=factor*-0.75
    else:
        factor=factor*0.25
    bfactor_test.append(factor)

print(bfactor_test)
print("average XplDDT after: ",mean(bfactor_test))

print("Interaction?")
if (mean(bfactor_test)>24):
    print("Probably Yes")
else:
    print("Probably No")



