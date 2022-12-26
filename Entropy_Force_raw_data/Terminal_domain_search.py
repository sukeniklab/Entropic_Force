# -*- coding: utf-8 -*-
"""
Search the terminal disordered region in the simulation database
"""
import pandas as pd
import json
import os
target_dir='F:\DATA_F\Entropy_Force_raw_data'
os.chdir(target_dir)
df_entry=pd.read_csv("Disprot_DP_simulation_stat.csv")
#%%
import requests
from urllib.parse import urljoin
disprot_url='https://disprot.org/api/'
# Load protein name and seqeunce. Search in Uniprot database for reference.
# If they are disprot entry search Disprot database.
def get_uniprot_id(disprot_name):
    retrive=requests.get(urljoin(disprot_url,disprot_name))
    if retrive.ok:
        return retrive.json()['acc']
    else:
        return False

#%%
uniprot_url='https://rest.uniprot.org/uniprot/'
def get_uniprot_seq(uniprot_id):
    if not uniprot_id:
        return
    retrive=requests.get(urljoin(uniprot_url,uniprot_id))
    if retrive.ok:
        return retrive.json()['sequence']['value']
    else:
        return
#%%
df_entry=pd.read_csv('Terminal_0923_cp.csv')      
for i in range(len(df_entry)):
    if pd.isnull(df_entry.loc[i,'Uniprot_seq']):
        if pd.notnull(df_entry.loc[i, 'Uniprot']):
            df_entry.loc[i,'Uniprot_seq']=get_uniprot_seq(df_entry.loc[i, 'Uniprot'])
#%%
for i in range(len(df_entry)):
    if pd.notnull(df_entry.loc[i, 'Uniprot_seq']):
        print (df_entry.loc[i, 'Uniprot_seq'])
#%%
### Retrive the Uniprot id for Disprot proteins
df_entry['Uniprot']=df_entry['Protein'].map(get_uniprot_id)
df_entry['Uniprot_seq']=df_entry['Uniprot'].map(get_uniprot_seq)
df_entry.to_csv('Terminal_auto.csv',index=False)
### Manully input other Uniport ID for other proteins
### Using https://www.uniprot.org/blast
#%%
df_entry=pd.read_csv('Terminal_all.csv')
for i in range(len(df_entry)):
    if pd.isnull(df_entry.loc[i,'Uniprot_seq']):
        if pd.notnull(df_entry.loc[i, 'Uniprot']):
            df_entry.loc[i,'Uniprot_seq']=get_uniprot_seq(df_entry.loc[i, 'Uniprot'])
df_entry['location']=df_entry.apply(lambda x: x['Uniprot_seq'].find(x['Sequence']) if pd.notnull(x['Uniprot_seq']) else None, axis=1)
df_entry['N_terminal']=(df_entry['location']<2)&(df_entry['location']>=0)
df_entry['C_terminal']=abs((df_entry['location']-(df_entry['Uniprot_seq'].str.len()-df_entry['Sequence'].str.len())))<=1
df_entry['C_terminal_test']=(df_entry['Uniprot_seq'].str.len()-df_entry['Sequence'].str.len())
#%%
from statannot import add_stat_annotation
df_entry['length']=df_entry['Sequence'].str.len()
df_entry['terminal']=df_entry['N_terminal'] | df_entry['C_terminal']
df_entry.to_csv('Terminal.csv',index=False)
#%%
N_terminal_count=df_entry.loc[df_entry['N_terminal']==True].count()
C_terminal_count=df_entry.loc[df_entry['C_terminal']==True].count()
terminal_count=df_entry.loc[df_entry['terminal']==True].count()
#%%
import seaborn as sns
import matplotlib.pyplot as plt
# Pyplot Setup
import matplotlib as mpl
plt.style.use(r'F:\DATA_F\JPCB_fig_publish\publish.mplstyple')
def save_to_svg(filename, root_path='F:\DATA_F\JPCB_fig_publish'):
    plt.savefig(os.path.join(root_path,filename))
#%%
import matplotlib.pyplot as plt
fig, ax = plt.subplots()

region=['terminal', 'non-terminal']
explode=(0.1, 0)
counts = [terminal_count,df_entry.count()-terminal_count]
#ax.set_title('Simulated Regions')
def autopct_format(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round((pct*total/100)))
        return '{p:.2f}% ({v:d})'.format(p=pct,v=val)
    return my_autopct
ax.pie(counts,explode=explode,labels=region,autopct=autopct_format(counts),textprops={'fontsize':40})
#ax.annotate(str(counts[0]),(0.5,0.7),xycoords='axes fraction',size=50)
#ax.annotate(str(counts[1]),(0.4,0.3),xycoords='axes fraction',size=50)
#ax.set_xlabel('IDR')
ax.tick_params(axis='both')
save_to_svg('Fig2A.svg')