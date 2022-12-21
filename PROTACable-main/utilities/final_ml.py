import pandas as pd
import glob
files=[i for i in glob.glob("./*all_ml_scores.csv")]
df = pd.concat(map(pd.read_csv, files), ignore_index=True)
df.to_csv('all_scores_ready.csv')
exit()
