import sys
import pdock
import os
import numpy as np
import csv
path_positive_colab="pdbs/positive/colabaf/"
path_positive_omega="pdbs/positive/omegafold/"
path_positive_alphafold="pdbs/positive/alphafold/"
path_negative="pdbs/negative/"
path_results="results/"
csv_heaver=['code', 'group','method', 'pDockQ', 'ppv','probable interaction?']
csv_data=[]
#chain_coords, chain_plddt = pdock.read_pdb(path_positive)

t=8 
#pdockq, ppv = pdock.calc_pdockq(chain_coords, chain_plddt, t)
def pdockq_path(path):
    for file in os.scandir(path):
        if file.is_file():
            name=file.name.split(".")[0]
            if name.split("_")[0].lower()=="negative":
                protein_name=name.split("_")[0]+name.split("_")[1]
                method="colabaf"
            else:
                protein_name=name.split("_")[0]
                method=path.split("/")[2]
            protein_name=protein_name.upper()
            path_file=path+name+".pdb"
            chain_coords, chain_plddt = pdock.read_pdb(path_file)
            
            t=8 
            if len(chain_coords.keys())<2:
                print('Only one chain in pdbfile', chain_coords)
                sys.exit()
            pdockq, ppv = pdock.calc_pdockq(chain_coords, chain_plddt, t)
            probable_interaction= "yes" if pdockq>=0.23 else "no" 
            csv_data.append([protein_name,path.split("/")[1],method,pdockq,ppv,probable_interaction])

            print('pDockQ =',np.round(pdockq,3),'for',path_file)
            print('This corresponds to a PPV of at least', ppv)

pdockq_path(path_negative)
pdockq_path(path_positive_colab)
pdockq_path(path_positive_omega)
pdockq_path(path_positive_alphafold)
with open('pdockq_results.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(csv_heaver)
    writer.writerows(csv_data)



