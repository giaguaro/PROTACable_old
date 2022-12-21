import pandas as pd
import glob, os
import sys

score_columns=["title",
    "icm_poi_lig_score",
    "icm_e3_lig_score",
    "rmsd_linker",
    "propose_score",
    "propose_csiz",
    "propose_vdw",
    "propose_ele",
    "propose_ehb",
    "propose_flw",
    "propose_bpa",
    "propose_bnpa",
    "propose_dint",
    "propose_dcost",
    "propose_erf"
              ]

score_fill=[
    str(sys.argv[1]),
    float(sys.argv[2]),
    float(sys.argv[3]),
    float(sys.argv[4]),
    float(sys.argv[5]),
    float(sys.argv[6]),
    float(sys.argv[7]),
    float(sys.argv[8]),
    float(sys.argv[9]),
    float(sys.argv[10]),
    float(sys.argv[11]),
    float(sys.argv[12]),
    float(sys.argv[13]),
    float(sys.argv[14]),
    float(sys.argv[15])
]
try:
    mmgbsa=pd.read_csv(str(sys.argv[16]), sep=",")
    if len(mmgbsa)==3:
        mmgbsa=mmgbsa.drop("title", axis=1)
        mmgbsa=mmgbsa.drop(index=[0,2],axis=1)
        #mmgbsa.to_csv(f"{str(sys.argv[1]).rsplit( ".", 1 )[ 0 ]}_cleaned.csv",index=False)
    elif len(mmgbsa)==1:
        mmgbsa=mmgbsa.drop("title", axis=1)
except:
    print("no mmgbsa scores found!")
df_scores=pd.DataFrame(dict(zip(score_columns,score_fill)), index=[0])

try:
    mmgbsa.index = df_scores.index
    df=pd.concat([df_scores,mmgbsa], axis=1)
except:
    df=df_scores.copy()
df.to_csv(f"{str(sys.argv[1])}_ml_scores.csv",index=False)

