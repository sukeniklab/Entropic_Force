import matplotlib.pyplot as plt
import mdtraj as md
import os, sys
import statistics as st
import numpy as np
import pandas as pd
target_dir=r'F:\alpha_yeast'
os.chdir(target_dir)
def read_prediction(filename):
    pdb = md.formats.pdb.pdbstructure.PdbStructure(open(filename))
    ca_list=[]
    predict=[]
    residue=[]
    location=[]
    location_counter=1
    for i in pdb.iter_atoms():
        if i.get_name()=='CA':
            ca_list.append(i)
            predict.append(i.get_temperature_factor())
            residue.append(i.residue_name)
            location.append(location_counter)
            location_counter+=1
    return residue,predict,location
def listsubdirectory(listname, targetdir):
    # Change the directory to the target dir
    os.chdir(targetdir)
    # Generate a full list of folders and files under the directory
    namestemp = os.listdir()
    # Only append folder to the list
    for n in namestemp:
        if (".csv" in n) or (".py" in n) or (".dat" in n) or (".xlsx" in n) or (".json" in n) or (".png" in n) or (
                ".jpg" in n) or ("pycache" in n) or (".svg" in n) or (".txt" in n) or (".pdf" in n) or (".fasta" in n):
            pass
        else:
            listname.append(n)
namestemp=[]
listsubdirectory(namestemp,target_dir)
data_list=len(namestemp)*[0]
print(len(namestemp))
for index,n in enumerate(namestemp):
    print(index)
    if not 'pdb' in n:
        print (n)
        continue
    residue_n,predict_n,location_n=read_prediction(n)
    data_list[index]=pd.DataFrame({'residue':residue_n,'alpha_fold_prediction':predict_n,'location':location_n})
    data_list[index]['protein_name']=n.split('-')[1]
full_df = pd.concat(data_list, ignore_index=True)    
full_df.to_csv('publish_1229.csv',index=False)