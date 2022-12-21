
import os

import sys
import glob
import rdkit
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.AllChem import  GetMorganFingerprintAsBitVect, GetErGFingerprint
import matplotlib.pyplot as plt
plt.ion()
import numpy as np
import pandas as pd
import seaborn as sns
import xgboost as xgb
from sklearn.compose import ColumnTransformer, make_column_transformer
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    f1_score,
    make_scorer,
    plot_confusion_matrix,
)
from sklearn.model_selection import (
    GridSearchCV,
    RandomizedSearchCV,
    cross_val_score,
    cross_validate,
    train_test_split,
)
from sklearn.pipeline import Pipeline, make_pipeline
from sklearn.preprocessing import OneHotEncoder, OrdinalEncoder, StandardScaler
from sklearn.svm import SVC
from catboost import CatBoostClassifier
from lightgbm.sklearn import LGBMClassifier
from sklearn.tree import DecisionTreeClassifier
from xgboost import XGBClassifier
from rdkit.Chem import CanonSmiles as canonicalize_smiles
from sklearn.preprocessing import RobustScaler
from sklearn.ensemble import VotingClassifier




df_indexed=pd.read_csv('indexed_final_protacpedia.csv')
df=pd.read_csv('all_experimental_scores.csv')
df_experimental=pd.read_csv('all_scores_ready.csv')
smiles_experimental=pd.read_csv('collected_smiles.smi', sep=' ', names=['mol','identifier'])


smiles_experimental['identifier']=smiles_experimental['identifier'].apply(lambda x: ("_").join(str(x).split(".")[0].split('_')[2:6]))
df_indexed=df_indexed.drop(['identifier'], axis=1)
df_indexed = df_indexed.rename(columns={'pp_plus_linkers_amber_renumbered_mol': 'mol_refined','E3 Binder SMILES': 'e3_smiles', 'Ligand SMILES': 'poi_smiles', 'PROTAC SMILES': 'protac_smiles','new_name':'identifier', 'Active/Inactive':'flag'})
df_indexed=df_indexed.drop(['mol_refined','Copy of mol', 'PROTACDB ID', 'protac_smiles', 'Best PROTAC', 'Cells', 'cLogP', 'Comments', 'Curator', 'Dc50', 'Dmax', 'e3_smiles', 'E3 Ligase', 'Ec50 of Ligand Cells', 'Ec50 of PROTAC Cells', 'exit_vector', 'Hbond acceptors', 'Hbond donors', 'Ic50 of Ligand', 'Ic50 of PROTAC', 'Ligand Name', 'poi_smiles', 'Copy of Ligand SMILES', 'Linker', 'Copy of Linker', 'Linker Type', 'linker_ha', 'linker_no', 'linker_rb', 'MW', 'Off Targets Reported', 'PATENT', 'Ligand PDB', 'Ligand ID', 'Pubmed', 'PROTAC Name', 'Proteomics Data Available', 'Secondary Pubmed', 'Status', 'Target', 'Tested A Non Binding E3 Control', 'Tested Competition With Ligand', 'Tested Engagement In Cells', 'Tested Proteaseome Inhibitor', 'Time', 'TPSA', 'ligand_pdb_copy', 'ligand_pdb'], axis=1)
df_indexed=df_indexed[['identifier','mol','flag']]
#'total_LYS_SASA'
df_indexed['flag'] = df_indexed['flag'].replace(['Inactive'],0)
df_indexed['flag'] = df_indexed['flag'].replace(['Active'],1)
df_indexed["identifier"] = df_indexed["identifier"].str.lower()

df=df.drop(['Unnamed: 0','Unnamed: 0.1','title.1','icm_poi_lig_score.1','icm_e3_lig_score.1','icm_e3_lig_score.1','title.2','icm_poi_lig_score.2','icm_e3_lig_score.2','icm_e3_lig_score.2','title.3','icm_poi_lig_score.3','icm_e3_lig_score.3','icm_e3_lig_score.3','title.4','icm_poi_lig_score.4','icm_e3_lig_score.4','icm_e3_lig_score.4','title.5','icm_poi_lig_score.5','icm_e3_lig_score.5','icm_e3_lig_score.5','title.6','icm_poi_lig_score.6','icm_e3_lig_score.6','icm_e3_lig_score.6','title.7','icm_poi_lig_score.7','icm_e3_lig_score.7','icm_e3_lig_score.7','title.8','icm_poi_lig_score.8','icm_e3_lig_score.8','icm_e3_lig_score.8','title.9','icm_poi_lig_score.9','icm_e3_lig_score.9','icm_e3_lig_score.9','title.10','icm_poi_lig_score.10','icm_e3_lig_score.10','icm_e3_lig_score.10','title.11','icm_poi_lig_score.11','icm_e3_lig_score.11','icm_e3_lig_score.11','title.12','icm_poi_lig_score.12','icm_e3_lig_score.12','icm_e3_lig_score.12','title.13','icm_poi_lig_score.13','icm_e3_lig_score.13','icm_e3_lig_score.13','title.14','icm_poi_lig_score.14','icm_e3_lig_score.14','icm_e3_lig_score.14','title.15','icm_poi_lig_score.15','icm_e3_lig_score.15','icm_e3_lig_score.15','title.16','icm_poi_lig_score.16','icm_e3_lig_score.16','icm_e3_lig_score.16','title.17','icm_poi_lig_score.17','icm_e3_lig_score.17','icm_e3_lig_score.17','title.18','icm_poi_lig_score.18','icm_e3_lig_score.18','icm_e3_lig_score.18','title.19','icm_poi_lig_score.19','icm_e3_lig_score.19','icm_e3_lig_score.19','rmsd_linker.1','rmsd_linker.2','rmsd_linker.3','rmsd_linker.4','rmsd_linker.5','rmsd_linker.6','rmsd_linker.7','rmsd_linker.8','rmsd_linker.9','rmsd_linker.10','rmsd_linker.11','rmsd_linker.12','rmsd_linker.13','rmsd_linker.14','rmsd_linker.15','rmsd_linker.16','rmsd_linker.17','rmsd_linker.18','rmsd_linker.19'], axis=1)


df_experimental=df_experimental.drop(['Unnamed: 0','Unnamed: 0.1','title.1','icm_poi_lig_score.1','icm_e3_lig_score.1','icm_e3_lig_score.1','title.2','icm_poi_lig_score.2','icm_e3_lig_score.2','icm_e3_lig_score.2','title.3','icm_poi_lig_score.3','icm_e3_lig_score.3','icm_e3_lig_score.3','title.4','icm_poi_lig_score.4','icm_e3_lig_score.4','icm_e3_lig_score.4','title.5','icm_poi_lig_score.5','icm_e3_lig_score.5','icm_e3_lig_score.5','title.6','icm_poi_lig_score.6','icm_e3_lig_score.6','icm_e3_lig_score.6','title.7','icm_poi_lig_score.7','icm_e3_lig_score.7','icm_e3_lig_score.7','title.8','icm_poi_lig_score.8','icm_e3_lig_score.8','icm_e3_lig_score.8','title.9','icm_poi_lig_score.9','icm_e3_lig_score.9','icm_e3_lig_score.9','title.10','icm_poi_lig_score.10','icm_e3_lig_score.10','icm_e3_lig_score.10','title.11','icm_poi_lig_score.11','icm_e3_lig_score.11','icm_e3_lig_score.11','title.12','icm_poi_lig_score.12','icm_e3_lig_score.12','icm_e3_lig_score.12','title.13','icm_poi_lig_score.13','icm_e3_lig_score.13','icm_e3_lig_score.13','title.14','icm_poi_lig_score.14','icm_e3_lig_score.14','icm_e3_lig_score.14','title.15','icm_poi_lig_score.15','icm_e3_lig_score.15','icm_e3_lig_score.15','title.16','icm_poi_lig_score.16','icm_e3_lig_score.16','icm_e3_lig_score.16','title.17','icm_poi_lig_score.17','icm_e3_lig_score.17','icm_e3_lig_score.17','title.18','icm_poi_lig_score.18','icm_e3_lig_score.18','icm_e3_lig_score.18','title.19','icm_poi_lig_score.19','icm_e3_lig_score.19','icm_e3_lig_score.19','rmsd_linker.1','rmsd_linker.2','rmsd_linker.3','rmsd_linker.4','rmsd_linker.5','rmsd_linker.6','rmsd_linker.7','rmsd_linker.8','rmsd_linker.9','rmsd_linker.10','rmsd_linker.11','rmsd_linker.12','rmsd_linker.13','rmsd_linker.14','rmsd_linker.15','rmsd_linker.16','rmsd_linker.17','rmsd_linker.18','rmsd_linker.19'], axis=1)

df_experimental['title']=df_experimental['title'].apply(lambda x: ("_").join(str(x).split('_')[0:4]))

df['title']=df['title'].apply(lambda x: ("_").join(str(x).split('_')[0:2]))

df_merged=pd.merge(df,df_indexed,how="left",left_on=["title"], right_on=["identifier"])
df_merged=df_merged.drop(['identifier'],axis=1)
df_merged.columns = df_merged.columns.str.replace("[.]", "_")
df_merged['mol_rdkit']=df_merged['mol'].apply(Chem.MolFromSmiles)
df_merged['Hacc_protac'] = df_merged['mol_rdkit'].apply(Chem.rdMolDescriptors.CalcNumHBA)
df_merged['Hdon_protac'] = df_merged['mol_rdkit'].apply(Chem.rdMolDescriptors.CalcNumHBD)
df_merged['TPSA_protac'] = df_merged['mol_rdkit'].apply(Chem.rdMolDescriptors.CalcTPSA)
df_merged['cLogP_protac'] = df_merged['mol_rdkit'].apply(Crippen.MolLogP)
df_merged['MW_protac'] = df_merged['mol_rdkit'].apply(Chem.rdMolDescriptors.CalcExactMolWt)
df_merged['FSP3'] = df_merged['mol_rdkit'].apply(Chem.rdMolDescriptors.CalcFractionCSP3)

df_merged=df_merged.drop(['mol_rdkit','mol','title'],axis=1)
df_merged=df_merged.fillna(0)


for idx,feat in enumerate(list(df_merged.columns)[1:97]):
    df_merged[f'{feat}_avg']=df_merged[list(df_merged.filter(regex=feat))].mean(axis=1)
    
for idx,feat in enumerate(list(df_merged.columns)[1:97]):
    df_merged[f'{feat}_sum']=df_merged[list(df_merged.filter(regex=feat))].sum(axis=1)
    
df_merged=df_merged.fillna(0)


df_merged2=pd.merge(df_experimental,smiles_experimental,how="left",left_on=["title"], right_on=["identifier"])


df_merged2.columns = df_merged2.columns.str.replace("[.]", "_")
df_merged2['mol'] = df_merged2['mol'].astype(str)
#df_merged2['canon_mol']=df_merged2['mol'].apply(Chem.CanonSmiles)
df_merged2['mol_rdkit']=df_merged2['mol'].apply(Chem.MolFromSmiles)
df_merged2 = df_merged2[df_merged2['mol_rdkit'].notna()]

df_merged2['Hacc_protac'] = df_merged2['mol_rdkit'].apply(Chem.rdMolDescriptors.CalcNumHBA)
df_merged2['Hdon_protac'] = df_merged2['mol_rdkit'].apply(Chem.rdMolDescriptors.CalcNumHBD)
df_merged2['TPSA_protac'] = df_merged2['mol_rdkit'].apply(Chem.rdMolDescriptors.CalcTPSA)
df_merged2['cLogP_protac'] = df_merged2['mol_rdkit'].apply(Crippen.MolLogP)
df_merged2['MW_protac'] = df_merged2['mol_rdkit'].apply(Chem.rdMolDescriptors.CalcExactMolWt)
df_merged2['FSP3'] = df_merged2['mol_rdkit'].apply(Chem.rdMolDescriptors.CalcFractionCSP3)
ref_df=df_merged2.copy()
df_merged2=df_merged2.drop(['identifier'],axis=1)
df_merged2=df_merged2.drop(['mol_rdkit','mol','title'],axis=1)
df_merged2=df_merged2.fillna(0)


for idx,feat in enumerate(list(df_merged2.columns)[1:97]):
    df_merged2[f'{feat}_avg']=df_merged2[list(df_merged2.filter(regex=feat))].mean(axis=1)
    
for idx,feat in enumerate(list(df_merged2.columns)[1:97]):
    df_merged2[f'{feat}_sum']=df_merged2[list(df_merged2.filter(regex=feat))].sum(axis=1)
    
df_merged2=df_merged2.fillna(0)

X_train,X_test, y_train,y_test = train_test_split(df_merged, df_merged['flag'], test_size=0.20, random_state=1,stratify=df_merged['flag'])

X_train.drop("flag",axis=1,inplace=True)
X_test.drop("flag",axis=1,inplace=True)


preprocessor = make_pipeline(RobustScaler())
# best paramters {'logisticregression__C': 0.0009696451837051425, 'logisticregression__penalty': 'l2', 
#                 'logisticregression__solver': 'newton-cg'}

# best paramters {'lgbmclassifier__colsample_bytree': 0.7790752106122303, 'lgbmclassifier__min_child_samples': 102, 
#                 'lgbmclassifier__min_child_weight': 1, 'lgbmclassifier__num_leaves': 45, 'lgbmclassifier__reg_alpha': 1, 
#                 'lgbmclassifier__reg_lambda': 10, 'lgbmclassifier__subsample': 0.5422907926702556}

# best paramters {'xgbclassifier__colsample_bytree': 0.3, 'xgbclassifier__learning_rate': 0.01, 
#                 'xgbclassifier__max_depth': 10, 'xgbclassifier__n_estimators': 100}

# best paramters {'randomforestclassifier__n_estimators': 466, 'randomforestclassifier__min_samples_split': 3, 
#                 'randomforestclassifier__min_samples_leaf': 2, 'randomforestclassifier__max_features': 'auto', 
#                 'randomforestclassifier__max_depth': 55, 'randomforestclassifier__bootstrap': False}

pipe_rf = make_pipeline(preprocessor, RandomForestClassifier(random_state=123, n_estimators= 600, 
                                                             min_samples_split=10, min_samples_leaf=2,max_features= 'sqrt',
                                                             max_depth= 15,bootstrap= False))
pipe_xgb = make_pipeline(
    preprocessor, XGBClassifier(random_state=123,colsample_bytree= 0.6, subsample=0.6,min_child_weight=3,
                                learning_rate= 0.01, max_depth=5, gamma=0.5, reg_lambda=0.8, 
                                n_estimators= 500)
)

pipe_lgbm = make_pipeline(preprocessor, LGBMClassifier(random_state=123,colsample_bytree= 0.7790752106122303, 
                                                       min_child_samples= 102, min_child_weight= 1, 
                                                       num_leaves=45, 
                                                       reg_alpha= 1, reg_lambda= 10, 
                                                       subsample= 0.5422907926702556))

pipe_catboost = make_pipeline(
    preprocessor, CatBoostClassifier(random_state=123, depth=5, iterations=90, learning_rate=0.04)
)


#pipe_lr = grid_lr.best_estimator_
#pipe_rf = grid_rf.best_estimator_
# pipe_xgb = grid_xgb.best_estimator_
# pipe_lgbm = grid_lgbm.best_estimator_
# pipe_catboost = pipe_catboost


classifiers = {
    #"logistic regression": pipe_lr,
    "random forest": pipe_rf,
    "XGBoost": pipe_xgb,
#    "LightGBM": pipe_lgbm,
    "CatBoost": pipe_catboost,
}

averaging_model_soft = VotingClassifier(
    list(classifiers.items()), voting="soft"
)  


averaging_model_hard = VotingClassifier(
    list(classifiers.items()), voting="hard"
)  


trained_soft=averaging_model_soft.fit(X_train, y_train)
#trained_hard=averaging_model_hard.fit(X_train, y_train)
score_soft=trained_soft.score(X_train, y_train)
score_soft


prediction_of_probability=trained_soft.predict_proba(df_merged2)

ref_df['prob_0'] = prediction_of_probability[:,0] 
ref_df['prob_1'] = prediction_of_probability[:,1]

sorted_df=ref_df.sort_values('prob_1', ascending=False)

chosen=list(sorted_df['identifier'].iloc[0:10])

prefix = []
unique_chosen = []
for i in chosen:
    if i.split('_')[0] not in prefix:
        prefix.append(i.split('_')[0])
        unique_chosen.append(i)
print(unique_chosen)   

retained_chosen=unique_chosen[0:3]

textfile = open("chosen.txt", "w")
for element in retained_chosen:
    textfile.write(element + "\n")
textfile.close()


