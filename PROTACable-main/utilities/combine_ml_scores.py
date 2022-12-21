import pandas as pd
import glob, os
import sys


all_filenames = [i for i in glob.glob(f'{str(sys.argv[1])}_*ml_scores.csv')]
list_of_column_names=[]

for i in all_filenames:
    df = pd.read_csv(i) 
    list_of_column_names.append(list(df.columns))
    

df_from_each_file = [pd.read_csv(f) for f in all_filenames]

frame = pd.concat(df_from_each_file, axis=1)

frame.to_csv(f'{str(sys.argv[1])}_all_ml_scores.csv')
