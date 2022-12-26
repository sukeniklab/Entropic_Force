# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd

df=pd.read_csv('Disprot_DP_ensemble_data.csv')
df=df[~(df['Protein']=='pumawild_full_repeat')]
df.to_csv('Disprot_DP_ensemble_data.csv', index=False)
df_2=pd.read_csv('Disprot_DP_entropic_force.csv')
df_2=df_2[~(df_2['Protein']=='pumawild_full_repeat')]
df_2.to_csv('Disprot_DP_entropic_force.csv', index=False)
