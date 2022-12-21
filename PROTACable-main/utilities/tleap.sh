#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=tLeap

source ~/.bashrc
conda activate py39


poi=$1
e3=$2

rm xx*
csplit $e3 '/@<TRIPOS>MOLECULE/' '{*}'
rm xx00

for x in xx*; do
	cat $x $poi > temp_model_${x%.*}.mol2
	obabel -imol2 temp_model_${x%.*}.mol2 -omol2 -O model_${x%.*}.mol2 -j
	obabel -imol2 model_${x%.*}.mol2 -opdb -O model_${x%.*}.pdb
	sed -i -e '/UNL/ s/ATOM  /HETATM/' model_${x%.*}.pdb
        grep -n -m 2 -e C99 -e O99 -e N99 model_${x%.*}.pdb >> in.file
        tr -d ' ' < in.file | cut -d':' -f1 | sort -u > out.file
        var1=$(awk 'NR==1' out.file);
        var2=$(awk 'NR==2' out.file);
        var1=$((var1-2));
        var2=$((var2-2));
        sed -i '/END/d' model_${x%.*}.pdb;
	sed -i '/MASTER/d' model_${x%.*}.pdb;
        string="CONECT ${var1} ${var2}"
        echo $string >> model_${x%.*}.pdb;
        echo "END" >> model_${x%.*}.pdb;
        rm in.file out.file
	sed -nr -e '/N99.*UNL/p' -e '/O99.*UNL/p' -e '/C99.*UNL/p' model_${x%.*}.mol2 >> indyc
	column -t indyc > indyc2
	var3=$(cat indyc2 | head -n1 | awk '{print $1;}')
	var4=$(cat indyc2 | head -n2 | awk '{print $1;}')
	var5=$(tail -n 1 model_${x%.*}.mol2| column -t tayl| awk '{print $1}')
	var6=$((var5+1))
	string2="  $var6  $var3  $var4    1"
	echo $string2 >> model_${x%.*}.mol2;
	rm indyc*
	rm temp_model_${x%.*}.mol2 $x
done;

echo "%{j}.out"
collective="model_xx01.pdb"

rm -f *trial* *prep*pdb *ligand* sqm.* ANTECHAMBER* ATOMTYPE* *log *noHET* *sed* *sslink* *txt *diagnosis


cp $collective ${collective%.*}_prep1.pdb

pdb4amber -i ${collective%.*}_prep1.pdb -o ${collective%.*}_r_prep1.pdb
grep "^HETATM" ${collective%.*}_prep1.pdb > ${collective%.*}_ligand.pdb

sed -i '/ACE/d' ${collective%.*}_r_prep1.pdb
sed -i '/NME/d' ${collective%.*}_r_prep1.pdb

sed -i '/HETATM/s/CL/Cl/g' ${collective%.*}_r_prep1.pdb
sed -i '/HETATM/s/BR/Br/g' ${collective%.*}_r_prep1.pdb

############################################################################################################################################################
sed -i 's/CSD/CYS/g' ${collective%.*}_r_prep1.pdb                 #3-SULFINOALANINE
sed -i 's/HYP/PRO/g' ${collective%.*}_r_prep1.pdb                 #4-HYDROXYPROLINE
sed -i 's/BMT/THR/g' ${collective%.*}_r_prep1.pdb                 #4-METHYL-4-[(E)-2-BUTENYL]-4,N-METHYL-THREONINE
sed -i 's/5HP/GLU/g' ${collective%.*}_r_prep1.pdb                 #5-HYDROXYPROLINE
sed -i 's/ABA/ALA/g' ${collective%.*}_r_prep1.pdb                 #ALPHA-AMINOBUTYRIC_ACID
sed -i 's/AIB/ALA/g' ${collective%.*}_r_prep1.pdb                 #ALPHA-AMINOISOBUTYRIC_ACID
sed -i 's/CSW/CYS/g' ${collective%.*}_r_prep1.pdb                 #CYSTEINE-S-DIOXIDE
sed -i 's/OCS/CYS/g' ${collective%.*}_r_prep1.pdb                 #CYSTEINESULFONIC_ACID
sed -i 's/DAL/ALA/g' ${collective%.*}_r_prep1.pdb                 #D-ALANINE
sed -i 's/DAR/ARG/g' ${collective%.*}_r_prep1.pdb                 #D-ARGININE
sed -i 's/DSG/ASN/g' ${collective%.*}_r_prep1.pdb                 #D-ASPARAGINE
sed -i 's/DSP/ASP/g' ${collective%.*}_r_prep1.pdb                 #D-ASPARTATE
sed -i 's/DCY/CYS/g' ${collective%.*}_r_prep1.pdb                 #D-CYSTEINE
sed -i 's/CRO/CRO/g' ${collective%.*}_r_prep1.pdb                 #DECARBOXY(PARAHYDROXYBENZYLIDENE-IMIDAZOLIDINONE)THREONINE
sed -i 's/DGL/GLU/g' ${collective%.*}_r_prep1.pdb                 #D-GLUTAMATE
sed -i 's/DGN/GLN/g' ${collective%.*}_r_prep1.pdb                 #D-GLUTAMINE
sed -i 's/DHI/HIS/g' ${collective%.*}_r_prep1.pdb                 #D-HISTIDINE
sed -i 's/DIL/ILE/g' ${collective%.*}_r_prep1.pdb                 #D-ISOLEUCINE
sed -i 's/DIV/VAL/g' ${collective%.*}_r_prep1.pdb                 #D-ISOVALINE
sed -i 's/DLE/LEU/g' ${collective%.*}_r_prep1.pdb                 #D-LEUCINE
sed -i 's/DLY/LYS/g' ${collective%.*}_r_prep1.pdb                 #D-LYSINE
sed -i 's/DPN/PHE/g' ${collective%.*}_r_prep1.pdb                 #D-PHENYLALANINE
sed -i 's/DPR/PRO/g' ${collective%.*}_r_prep1.pdb                 #D-PROLINE
sed -i 's/DSN/SER/g' ${collective%.*}_r_prep1.pdb                 #D-SERINE
sed -i 's/DTH/THR/g' ${collective%.*}_r_prep1.pdb                 #D-THREONINE
sed -i 's/DTR/DTR/g' ${collective%.*}_r_prep1.pdb                 #D-TRYPTOPHANE
sed -i 's/DTY/TYR/g' ${collective%.*}_r_prep1.pdb                 #D-TYROSINE
sed -i 's/DVA/VAL/g' ${collective%.*}_r_prep1.pdb                 #D-VALINE
sed -i 's/CGU/GLU/g' ${collective%.*}_r_prep1.pdb                 #GAMMA-CARBOXY-GLUTAMIC_ACID
sed -i 's/KCX/LYS/g' ${collective%.*}_r_prep1.pdb                 #LYSINE_NZ-CARBOXYLIC_ACID
sed -i 's/LLP/LYS/g' ${collective%.*}_r_prep1.pdb                 #LYSINE-PYRIDOXAL-5'-PHOSPHATE
sed -i 's/CXM/MET/g' ${collective%.*}_r_prep1.pdb                 #N-CARBOXYMETHIONINE
sed -i 's/FME/MET/g' ${collective%.*}_r_prep1.pdb                 #N-FORMYLMETHIONINE
sed -i 's/MLE/LEU/g' ${collective%.*}_r_prep1.pdb                 #N-METHYLLEUCINE
sed -i 's/MVA/VAL/g' ${collective%.*}_r_prep1.pdb                 #N-METHYLVALINE
sed -i 's/NLE/LEU/g' ${collective%.*}_r_prep1.pdb                 #NORLEUCINE
sed -i 's/PTR/TYR/g' ${collective%.*}_r_prep1.pdb                 #O-PHOSPHOTYROSINE
sed -i 's/ORN/ALA/g' ${collective%.*}_r_prep1.pdb                 #ORNITHINE
sed -i 's/SEP/SER/g' ${collective%.*}_r_prep1.pdb                 #PHOSPHOSERINE
sed -i 's/TPO/THR/g' ${collective%.*}_r_prep1.pdb                 #PHOSPHOTHREONINE
sed -i 's/PCA/GLU/g' ${collective%.*}_r_prep1.pdb                 #PYROGLUTAMIC_ACID
sed -i 's/SAR/GLY/g' ${collective%.*}_r_prep1.pdb                 #SARCOSINE
sed -i 's/CEA/CYS/g' ${collective%.*}_r_prep1.pdb                 #S-HYDROXY-CYSTEINE
sed -i 's/CSO/CYS/g' ${collective%.*}_r_prep1.pdb                 #S-HYDROXYCYSTEINE
sed -i 's/CSS/CYS/g' ${collective%.*}_r_prep1.pdb                 #S-MERCAPTOCYSTEINE
sed -i 's/CSX/CYS/g' ${collective%.*}_r_prep1.pdb                 #S-OXY_CYSTEINE
sed -i 's/CME/CYS/g' ${collective%.*}_r_prep1.pdb                 #S,S-(2-HYDROXYETHYL)THIOCYSTEINE
sed -i 's/TYS/TYR/g' ${collective%.*}_r_prep1.pdb                 #SULFONATED_TYROSINE
sed -i 's/TPQ/PHE/g' ${collective%.*}_r_prep1.pdb                 #TOPO-QUINONE
sed -i 's/STY/TYR/g' ${collective%.*}_r_prep1.pdb                 #TYROSINE-O-SULPHONIC_ACID
sed -i 's/CCS/CYS/g' ${collective%.*}_r_prep1.pdb                 #https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py
sed -i 's/CALA/ALA/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CARG/ARG/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CASN/ASN/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CASP/ASP/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CCYS/CYS/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CCYX/CYX/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CGLN/GLN/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CGLU/GLU/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CGLY/GLY/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CHID/HID/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CHIE/HIE/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CHIP/HIP/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CHYP/HYP/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CILE/ILE/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CLEU/LEU/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CLYS/LYS/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CMET/MET/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CPHE/PHE/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CPRO/PRO/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CSER/SER/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CTHR/THR/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CTRP/TRP/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CTYR/TYR/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CVAL/VAL/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NALA/ALA/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NARG/ARG/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NASN/ASN/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NASP/ASP/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NCYS/CYS/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NCYX/CYX/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NGLN/GLN/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NGLU/GLU/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NGLY/GLY/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NHID/HID/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NHIE/HIE/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NHIP/HIP/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NILE/ILE/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NLEU/LEU/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NLYS/LYS/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NMET/MET/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NPHE/PHE/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NPRO/PRO/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NSER/SER/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NTHR/THR/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NTRP/TRP/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NTYR/TYR/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NVAL/VAL/g' ${collective%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/DAS/ASP/g' ${collective%.*}_r_prep1.pdb                 #GUESS https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py
sed -i 's/CAF/CYS/g' ${collective%.*}_r_prep1.pdb                 #GUESS https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py ## S-DIMETHYLARSINOYL-CYSTEINE
sed -i 's/CAS/CYS/g' ${collective%.*}_r_prep1.pdb                 #GUESS https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py ## S-(DIMETHYLARSENIC)CYSTEINE
#sed -i 's/AIB/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/ALA/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                  ALA
#sed -i 's/ALM/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/AYA/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/BNN/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/CHG/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/CSD/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/DAL/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/DHA/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/DNP/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/FLA/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/HAC/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/PRR/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/MAA/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/TIH/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/TPQ/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/0CS/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    0CS ALA  3-[(S)-HYDROPEROXYSULFINYL]-L-ALANINE
#sed -i 's/2BU/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    2BU ADE
#sed -i 's/2OP/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    2OP (2S  2-HYDROXYPROPANAL
#sed -i 's/4F3/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    4F3 ALA  CYCLIZED
#sed -i 's/AA4/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    AA4 ALA  2-AMINO-5-HYDROXYPENTANOIC ACID
#sed -i 's/ABA/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    ABA ALA  ALPHA-AMINOBUTYRIC ACID
#sed -i 's/AHO/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    AHO ALA  N-ACETYL-N-HYDROXY-L-ORNITHINE
#sed -i 's/AHP/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    AHP ALA  2-AMINO-HEPTANOIC ACID
#sed -i 's/AIB/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    AIB ALA  ALPHA-AMINOISOBUTYRIC ACID
#sed -i 's/ALA/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    ALA ALA
#sed -i 's/ALC/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    ALC ALA  2-AMINO-3-CYCLOHEXYL-PROPIONIC ACID
#sed -i 's/ALM/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    ALM ALA  1-METHYL-ALANINAL
#sed -i 's/ALN/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    ALN ALA  NAPHTHALEN-2-YL-3-ALANINE
#sed -i 's/ALS/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    ALS ALA  2-AMINO-3-OXO-4-SULFO-BUTYRIC ACID
#sed -i 's/ALT/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    ALT ALA  THIOALANINE
#sed -i 's/AP7/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    AP7 ADE
#sed -i 's/APH/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    APH ALA  P-AMIDINOPHENYL-3-ALANINE
#sed -i 's/AYA/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    AYA ALA  N-ACETYLALANINE
#sed -i 's/AYG/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    AYG ALA
#sed -i 's/B2A/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    B2A ALA  ALANINE BORONIC ACID
#sed -i 's/B3A/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    B3A ALA  (3S)-3-AMINOBUTANOIC ACID
#sed -i 's/BAL/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    BAL ALA  BETA-ALANINE
#sed -i 's/BNN/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    BNN ALA  ACETYL-P-AMIDINOPHENYLALANINE
#sed -i 's/C12/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    C12 ALA
#sed -i 's/C99/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    C99 ALA
#sed -i 's/CAB/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    CAB ALA  4-CARBOXY-4-AMINOBUTANAL
#sed -i 's/CH6/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    CH6 ALA
#sed -i 's/CH7/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    CH7 ALA
#sed -i 's/CLB/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    CLB ALA
#sed -i 's/CLD/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    CLD ALA
#sed -i 's/CLV/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    CLV ALA
#sed -i 's/CQR/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    CQR ALA
#sed -i 's/CR2/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    CR2 ALA  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/CR5/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    CR5 ALA
#sed -i 's/CR7/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    CR7 ALA
#sed -i 's/CR8/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    CR8 ALA
#sed -i 's/CRK/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    CRK ALA
#sed -i 's/CRW/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    CRW ALA
#sed -i 's/CRX/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    CRX ALA
#sed -i 's/CSI/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    CSI ALA
#sed -i 's/CSY/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    CSY ALA  MODIFIED TYROSINE COMPLEX
#sed -i 's/CWR/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    CWR ALA
#sed -i 's/DAB/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    DAB ALA  24-DIAMINOBUTYRIC ACID
#sed -i 's/DAL/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    DAL ALA  D-ALANINE
#sed -i 's/DAM/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    DAM ALA  N-METHYL-ALPHA-BETA-DEHYDROALANINE
#sed -i 's/DBU/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    DBU ALA  (2E)-2-AMINOBUT-2-ENOIC ACID
#sed -i 's/DBZ/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    DBZ ALA  3-(BENZOYLAMINO)-L-ALANINE
#sed -i 's/DHA/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    DHA ALA  2-AMINO-ACRYLIC ACID
#sed -i 's/DPP/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    DPP ALA  DIAMMINOPROPANOIC ACID
#sed -i 's/FGL/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    FGL ALA  2-AMINOPROPANEDIOIC ACID
#sed -i 's/HHK/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    HHK ALA  (2S)-28-DIAMINOOCTANOIC ACID
#sed -i 's/HMF/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    HMF ALA  2-AMINO-4-PHENYL-BUTYRIC ACID
#sed -i 's/IAM/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    IAM ALA  4-[(ISOPROPYLAMINO)METHYL]PHENYLALANINE
#sed -i 's/IGL/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    IGL ALA  ALPHA-AMINO-2-INDANACETIC ACID
#sed -i 's/KYN/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    KYN ALA  KYNURENINE
#sed -i 's/LAL/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    LAL ALA  NN-DIMETHYL-L-ALANINE
#sed -i 's/MAA/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    MAA ALA  N-METHYLALANINE
#sed -i 's/MDO/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    MDO ALA
#sed -i 's/MFC/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    MFC ALA  CYCLIZED
#sed -i 's/NAL/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    NAL ALA  BETA-(2-NAPHTHYL)-ALANINE
#sed -i 's/NAM/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    NAM ALA  NAM NAPTHYLAMINOALANINE
#sed -i 's/NCB/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    NCB ALA  CHEMICAL MODIFICATION
#sed -i 's/NRQ/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    NRQ ALA
#sed -i 's/NYC/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    NYC ALA
#sed -i 's/ORN/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    ORN ALA  ORNITHINE
#sed -i 's/PIA/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    PIA ALA  FUSION OF ALA 65 TYR 66 GLY 67
#sed -i 's/PRR/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    PRR ALA  3-(METHYL-PYRIDINIUM)ALANINE
#sed -i 's/PYA/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    PYA ALA  3-(110-PHENANTHROL-2-YL)-L-ALANINE
#sed -i 's/PYC/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    PYC ALA  PYRROLE-2-CARBOXYLATE
#sed -i 's/PYT/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    PYT ALA  MODIFIED ALANINE
#sed -i 's/RC7/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    RC7 ALA
#sed -i 's/SEC/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    SEC ALA  2-AMINO-3-SELENINO-PROPIONIC ACID
#sed -i 's/SIC/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    SIC ALA
#sed -i 's/SUI/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    SUI ALA
#sed -i 's/TIH/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    TIH ALA  BETA(2-THIENYL)ALANINE
#sed -i 's/TPQ/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    TPQ ALA  245-TRIHYDROXYPHENYLALANINE
#sed -i 's/UMA/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    UMA ALA
#sed -i 's/X9Q/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    X9Q ALA
#sed -i 's/XXY/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    XXY ALA
#sed -i 's/XYG/ALA/g' ${collective%.*}_r_prep1.pdb                 ###                    XYG ALA
#sed -i 's/BCS/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/BUC/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/C5C/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/C6C/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CCS/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CEA/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CME/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSO/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSP/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSS/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSX/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSW/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CY1/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CY3/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CYG/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i 's/CYM/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CYS/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  CYS
#sed -i 's/CYQ/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/DCY/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/EFC/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/OCS/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/PEC/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/PR3/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SCH/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SCS/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SCY/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SHC/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SMC/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SOC/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/5CS/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    5CS CYS
#sed -i 's/AGT/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    AGT CYS  AGMATINE-CYSTEINE ADDUCT
#sed -i 's/BBC/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    BBC CYS
#sed -i 's/BCS/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    BCS CYS  BENZYLCYSTEINE
#sed -i 's/BCX/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    BCX CYS  BETA-3-CYSTEINE
#sed -i 's/BPE/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    BPE CYS
#sed -i 's/BUC/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    BUC CYS  SS-BUTYLTHIOCYSTEINE
#sed -i 's/C3Y/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    C3Y CYS  MODIFIED CYSTEINE
#sed -i 's/C5C/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    C5C CYS  S-CYCLOPENTYL THIOCYSTEINE
#sed -i 's/C6C/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    C6C CYS  S-CYCLOHEXYL THIOCYSTEINE
#sed -i 's/CAF/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CAF CYS  S-DIMETHYLARSINOYL-CYSTEINE
#sed -i 's/CAS/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CAS CYS  S-(DIMETHYLARSENIC)CYSTEINE
#sed -i 's/CCS/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CCS CYS  CARBOXYMETHYLATED CYSTEINE
#sed -i 's/CME/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CME CYS  MODIFIED CYSTEINE
#sed -i 's/CML/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CML CYS
#sed -i 's/CMT/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CMT CYS  O-METHYLCYSTEINE
#sed -i 's/CS1/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CS1 CYS  S-(2-ANILINYL-SULFANYL)-CYSTEINE
#sed -i 's/CS3/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CS3 CYS
#sed -i 's/CS4/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CS4 CYS
#sed -i 's/CSA/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CSA CYS  S-ACETONYLCYSTEIN
#sed -i 's/CSB/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CSB CYS  CYS BOUND TO LEAD ION
#sed -i 's/CSD/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CSD CYS  3-SULFINOALANINE
#sed -i 's/CSE/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CSE CYS  SELENOCYSTEINE
#sed -i 's/CSO/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CSO CYS  INE S-HYDROXYCYSTEINE
#sed -i 's/CSR/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CSR CYS  S-ARSONOCYSTEINE
#sed -i 's/CSS/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CSS CYS  13-THIAZOLE-4-CARBOXYLIC ACID
#sed -i 's/CSU/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CSU CYS  CYSTEINE-S-SULFONIC ACID
#sed -i 's/CSW/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CSW CYS  CYSTEINE-S-DIOXIDE
#sed -i 's/CSX/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CSX CYS  OXOCYSTEINE
#sed -i 's/CSZ/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CSZ CYS  S-SELANYL CYSTEINE
#sed -i 's/CY0/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CY0 CYS  MODIFIED CYSTEINE
#sed -i 's/CY1/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CY1 CYS  ACETAMIDOMETHYLCYSTEINE
#sed -i 's/CY3/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CY3 CYS  2-AMINO-3-MERCAPTO-PROPIONAMIDE
#sed -i 's/CY4/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CY4 CYS  S-BUTYRYL-CYSTEIN
#sed -i 's/CY7/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CY7 CYS  MODIFIED CYSTEINE
#sed -i 's/CYD/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CYD CYS
#sed -i 's/CYF/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CYF CYS  FLUORESCEIN LABELLED CYS380 (P14)
#sed -i 's/CYG/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CYG CYS
#sed -i 's/CYQ/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CYQ CYS
#sed -i 's/CYR/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CYR CYS
#sed -i 's/CYS/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CYS CYS
#sed -i 's/CZ2/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CZ2 CYS  S-(DIHYDROXYARSINO)CYSTEINE
#sed -i 's/CZZ/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CZZ CYS  THIARSAHYDROXY-CYSTEINE
#sed -i 's/DCY/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    DCY CYS  D-CYSTEINE
#sed -i 's/DYS/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    DYS CYS
#sed -i 's/EFC/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    EFC CYS  SS-(2-FLUOROETHYL)THIOCYSTEINE
#sed -i 's/FOE/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    FOE CYS
#sed -i 's/GT9/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    GT9 CYS  SG ALKYLATED
#sed -i 's/GYC/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    GYC CYS
#sed -i 's/HTI/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    HTI CYS
#sed -i 's/KOR/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    KOR CYS  MODIFIED CYSTEINE
#sed -i 's/M0H/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    M0H CYS  S-(HYDROXYMETHYL)-L-CYSTEINE
#sed -i 's/MCS/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    MCS CYS  MALONYLCYSTEINE
#sed -i 's/NPH/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    NPH CYS
#sed -i 's/NYS/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    NYS CYS
#sed -i 's/OCS/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    OCS CYS  CYSTEINE SULFONIC ACID
#sed -i 's/OCY/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    OCY CYS  HYDROXYETHYLCYSTEINE
#sed -i 's/P1L/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    P1L CYS  S-PALMITOYL CYSTEINE
#sed -i 's/PBB/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    PBB CYS  S-(4-BROMOBENZYL)CYSTEINE
#sed -i 's/PEC/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    PEC CYS  SS-PENTYLTHIOCYSTEINE
#sed -i 's/PR3/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    PR3 CYS  INE DTT-CYSTEINE
#sed -i 's/PYX/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    PYX CYS  S-[S-THIOPYRIDOXAMINYL]CYSTEINE
#sed -i 's/R1A/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    R1A CYS
#sed -i 's/R1B/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    R1B CYS
#sed -i 's/R1F/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    R1F CYS
#sed -i 's/R7A/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    R7A CYS
#sed -i 's/RCY/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    RCY CYS
#sed -i 's/SAH/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    SAH CYS  S-ADENOSYL-L-HOMOCYSTEINE
#sed -i 's/SC2/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    SC2 CYS  N-ACETYL-L-CYSTEINE
#sed -i 's/SCH/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    SCH CYS  S-METHYL THIOCYSTEINE GROUP
#sed -i 's/SCS/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    SCS CYS  MODIFIED CYSTEINE
#sed -i 's/SCY/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    SCY CYS  CETYLATED CYSTEINE
#sed -i 's/SHC/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    SHC CYS  S-HEXYLCYSTEINE
#sed -i 's/SMC/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    SMC CYS  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/SNC/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    SNC CYS  S-NITROSO CYSTEINE
#sed -i 's/SOC/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    SOC CYS  DIOXYSELENOCYSTEINE
#sed -i 's/TEE/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    TEE CYS  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/TNB/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    TNB CYS  S-(236-TRINITROPHENYL)CYSTEINE
#sed -i 's/TYX/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    TYX CYS  S-(2-ANILINO-2-OXOETHYL)-L-CYSTEINE
#sed -i 's/YCM/CYS/g' ${collective%.*}_r_prep1.pdb                 ###                    YCM CYS  S-(2-AMINO-2-OXOETHYL)-L-CYSTEINE
#sed -i 's/2AS/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASA/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASB/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASK/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASL/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASP/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                  ASP
#sed -i 's/ASQ/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/BHD/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/DAS/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/DSP/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/3MD/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    3MD ASP  2S3S-3-METHYLASPARTIC ACID
#sed -i 's/A0A/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    A0A ASP  ASPARTYL-FORMYL MIXED ANHYDRIDE
#sed -i 's/ACB/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    ACB ASP  3-METHYL-ASPARTIC ACID
#sed -i 's/AKL/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    AKL ASP  3-AMINO-5-CHLORO-4-OXOPENTANOIC ACID
#sed -i 's/ASA/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    ASA ASP  ASPARTIC ALDEHYDE
#sed -i 's/ASB/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    ASB ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
#sed -i 's/ASI/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    ASI ASP  L-ISO-ASPARTATE
#sed -i 's/ASK/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    ASK ASP  DEHYDROXYMETHYLASPARTIC ACID
#sed -i 's/ASL/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    ASL ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
#sed -i 's/ASP/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    ASP ASP
#sed -i 's/B3D/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    B3D ASP  3-AMINOPENTANEDIOIC ACID
#sed -i 's/BFD/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    BFD ASP  ASPARTATE BERYLLIUM FLUORIDE
#sed -i 's/BHD/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    BHD ASP  BETA-HYDROXYASPARTIC ACID
#sed -i 's/DAS/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    DAS ASP  D-ASPARTIC ACID
#sed -i 's/DMK/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    DMK ASP  DIMETHYL ASPARTIC ACID
#sed -i 's/IAS/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    IAS ASP  ASPARTYL GROUP
#sed -i 's/OHS/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    OHS ASP  O-(CARBOXYSULFANYL)-4-OXO-L-HOMOSERINE
#sed -i 's/OXX/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    OXX ASP  OXALYL-ASPARTYL ANHYDRIDE
#sed -i 's/PHD/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    PHD ASP  2-AMINO-4-OXO-4-PHOSPHONOOXY-BUTYRIC ACID
#sed -i 's/SNN/ASP/g' ${collective%.*}_r_prep1.pdb                 ###                    SNN ASP  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/5HP/GLU/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/CGU/GLU/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/DGL/GLU/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/GGL/GLU/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/GLU/GLU/g' ${collective%.*}_r_prep1.pdb                 ###                  GLU
#sed -i 's/GMA/GLU/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/PCA/GLU/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/AB7/GLU/g' ${collective%.*}_r_prep1.pdb                 ###                    AB7 GLU  ALPHA-AMINOBUTYRIC ACID
#sed -i 's/AR4/GLU/g' ${collective%.*}_r_prep1.pdb                 ###                    AR4 GLU
#sed -i 's/B3E/GLU/g' ${collective%.*}_r_prep1.pdb                 ###                    B3E GLU  (3S)-3-AMINOHEXANEDIOIC ACID
#sed -i 's/CGU/GLU/g' ${collective%.*}_r_prep1.pdb                 ###                    CGU GLU  CARBOXYLATION OF THE CG ATOM
#sed -i 's/DGL/GLU/g' ${collective%.*}_r_prep1.pdb                 ###                    DGL GLU  D-GLU
#sed -i 's/GLU/GLU/g' ${collective%.*}_r_prep1.pdb                 ###                    GLU GLU
#sed -i 's/GMA/GLU/g' ${collective%.*}_r_prep1.pdb                 ###                    GMA GLU  1-AMIDO-GLUTAMIC ACID
#sed -i 's/ILG/GLU/g' ${collective%.*}_r_prep1.pdb                 ###                    ILG GLU  GLU LINKED TO NEXT RESIDUE VIA CG
#sed -i 's/LME/GLU/g' ${collective%.*}_r_prep1.pdb                 ###                    LME GLU  (3R)-3-METHYL-L-GLUTAMIC ACID
#sed -i 's/MEG/GLU/g' ${collective%.*}_r_prep1.pdb                 ###                    MEG GLU  (2S3R)-3-METHYL-GLUTAMIC ACID
#sed -i 's/DAH/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/DPN/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/HPQ/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/PHE/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                  PHE
#sed -i 's/PHI/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/PHL/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/1PA/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    1PA PHE  PHENYLMETHYLACETIC ACID ALANINE
#sed -i 's/23F/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    23F PHE  (2Z)-2-AMINO-3-PHENYLACRYLIC ACID
#sed -i 's/4PH/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    4PH PHE  4-METHYL-L-PHENYLALANINE
#sed -i 's/B2F/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    B2F PHE  PHENYLALANINE BORONIC ACID
#sed -i 's/BIF/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    BIF PHE
#sed -i 's/CHS/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    CHS PHE  4-AMINO-5-CYCLOHEXYL-3-HYDROXY-PENTANOIC AC
#sed -i 's/DAH/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    DAH PHE  34-DIHYDROXYDAHNYLALANINE
#sed -i 's/DPH/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    DPH PHE  DEAMINO-METHYL-PHENYLALANINE
#sed -i 's/DPN/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    DPN PHE  D-CONFIGURATION
#sed -i 's/FCL/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    FCL PHE  3-CHLORO-L-PHENYLALANINE
#sed -i 's/FOG/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    FOG PHE  PHENYLALANINOYL-[1-HYDROXY]-2-PROPYLENE
#sed -i 's/FRF/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    FRF PHE  PHE FOLLOWED BY REDUCED PHE
#sed -i 's/HPE/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    HPE PHE  HOMOPHENYLALANINE
#sed -i 's/HPH/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    HPH PHE  PHENYLALANINOL GROUP
#sed -i 's/HPQ/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    HPQ PHE  HOMOPHENYLALANINYLMETHANE
#sed -i 's/MEA/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    MEA PHE  N-METHYLPHENYLALANINE
#sed -i 's/MTY/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    MTY PHE  3-HYDROXYPHENYLALANINE
#sed -i 's/NFA/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    NFA PHE  MODIFIED PHENYLALANINE
#sed -i 's/PBF/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    PBF PHE  PARA-(BENZOYL)-PHENYLALANINE
#sed -i 's/PCS/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    PCS PHE  PHENYLALANYLMETHYLCHLORIDE
#sed -i 's/PF5/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    PF5 PHE  23456-PENTAFLUORO-L-PHENYLALANINE
#sed -i 's/PFF/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    PFF PHE  4-FLUORO-L-PHENYLALANINE
#sed -i 's/PHA/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    PHA PHE  PHENYLALANINAL
#sed -i 's/PHE/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    PHE PHE
#sed -i 's/PHI/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    PHI PHE  IODO-PHENYLALANINE
#sed -i 's/PHL/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    PHL PHE  L-PHENYLALANINOL
#sed -i 's/PHM/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    PHM PHE  PHENYLALANYLMETHANE
#sed -i 's/PM3/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    PM3 PHE
#sed -i 's/PPN/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    PPN PHE  THE LIGAND IS A PARA-NITRO-PHENYLALANINE
#sed -i 's/PRQ/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    PRQ PHE  PHENYLALANINE
#sed -i 's/PSA/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    PSA PHE
#sed -i 's/SMF/PHE/g' ${collective%.*}_r_prep1.pdb                 ###                    SMF PHE  4-SULFOMETHYL-L-PHENYLALANINE
#sed -i 's/GL3/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/GLY/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                  GLY
#sed -i 's/GLZ/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/GSC/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/MPQ/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/MSA/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/NMC/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/SAR/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/ACY/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                    ACY GLY  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/CHG/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                    CHG GLY  CYCLOHEXYL GLYCINE
#sed -i 's/CHP/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                    CHP GLY  3-CHLORO-4-HYDROXYPHENYLGLYCINE
#sed -i 's/GHP/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                    GHP GLY  4-HYDROXYPHENYLGLYCINE
#sed -i 's/GL3/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                    GL3 GLY  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/GLY/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                    GLY GLY
#sed -i 's/GLZ/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                    GLZ GLY  AMINO-ACETALDEHYDE
#sed -i 's/GYS/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                    GYS GLY
#sed -i 's/IPG/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                    IPG GLY  N-ISOPROPYL GLYCINE
#sed -i 's/MEU/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                    MEU GLY  O-METHYL-GLYCINE
#sed -i 's/MPQ/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                    MPQ GLY  N-METHYL-ALPHA-PHENYL-GLYCINE
#sed -i 's/MSA/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                    MSA GLY  (2-S-METHYL) SARCOSINE
#sed -i 's/NMC/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                    NMC GLY  N-CYCLOPROPYLMETHYL GLYCINE
#sed -i 's/PG9/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                    PG9 GLY  D-PHENYLGLYCINE
#sed -i 's/SAR/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                    SAR GLY  SARCOSINE
#sed -i 's/SHP/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                    SHP GLY  (4-HYDROXYMALTOSEPHENYL)GLYCINE
#sed -i 's/TBG/GLY/g' ${collective%.*}_r_prep1.pdb                 ###                    TBG GLY  T-BUTYL GLYCINE
#sed -i 's/3AH/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/DHI/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/HIC/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/HIS/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                  HIS
#sed -i 's/MHS/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/NEM/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/NEP/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/HID/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                  single delta N protonation
#sed -i 's/HIE/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                  single epsilon N protonation
#sed -i 's/3AH/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                    3AH HIS
#sed -i 's/DDE/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                    DDE HIS
#sed -i 's/DHI/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                    DHI HIS  D-HISTIDINE
#sed -i 's/HIA/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                    HIA HIS  L-HISTIDINE AMIDE
#sed -i 's/HIC/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                    HIC HIS  4-METHYL-HISTIDINE
#sed -i 's/HIP/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                    HIP HIS  ND1-PHOSPHONOHISTIDINE...or commonly used doubly protonated state
#sed -i 's/HIQ/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                    HIQ HIS  MODIFIED HISTIDINE
#sed -i 's/HIS/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                    HIS HIS
#sed -i 's/HSO/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                    HSO HIS  HISTIDINOL
#sed -i 's/MHS/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                    MHS HIS  1-N-METHYLHISTIDINE
#sed -i 's/NEP/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                    NEP HIS  N1-PHOSPHONOHISTIDINE
#sed -i 's/NZH/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                    NZH HIS
#sed -i 's/OHI/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                    OHI HIS  3-(2-OXO-2H-IMIDAZOL-4-YL)-L-ALANINE
#sed -i 's/PSH/HIS/g' ${collective%.*}_r_prep1.pdb                 ###                    PSH HIS  1-THIOPHOSPHONO-L-HISTIDINE
#sed -i 's/DIL/ILE/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ILE
#sed -i 's/IIL/ILE/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ILE
#sed -i 's/ILE/ILE/g' ${collective%.*}_r_prep1.pdb                 ###                  ILE
#sed -i 's/B2I/ILE/g' ${collective%.*}_r_prep1.pdb                 ###                    B2I ILE  ISOLEUCINE BORONIC ACID
#sed -i 's/DIL/ILE/g' ${collective%.*}_r_prep1.pdb                 ###                    DIL ILE  D-ISOLEUCINE
#sed -i 's/IIL/ILE/g' ${collective%.*}_r_prep1.pdb                 ###                    IIL ILE  ISO-ISOLEUCINE
#sed -i 's/ILE/ILE/g' ${collective%.*}_r_prep1.pdb                 ###                    ILE ILE
#sed -i 's/ILX/ILE/g' ${collective%.*}_r_prep1.pdb                 ###                    ILX ILE  45-DIHYDROXYISOLEUCINE
#sed -i 's/IML/ILE/g' ${collective%.*}_r_prep1.pdb                 ###                    IML ILE  N-METHYLATED
#sed -i 's/ALY/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/DLY/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/KCX/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/LLP/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/LLY/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/LYM/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/LYS/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                  LYS
#sed -i 's/LYZ/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/MLY/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/SHR/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/TRG/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/6CL/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    6CL LYS  6-CARBOXYLYSINE
#sed -i 's/ALY/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    ALY LYS  N(6)-ACETYLLYSINE
#sed -i 's/API/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    API LYS  26-DIAMINOPIMELIC ACID
#sed -i 's/APK/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    APK LYS
#sed -i 's/AZK/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    AZK LYS  (2S)-2-AMINO-6-TRIAZANYLHEXAN-1-OL
#sed -i 's/B3K/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    B3K LYS  (3S)-37-DIAMINOHEPTANOIC ACID
#sed -i 's/BLY/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    BLY LYS  LYSINE BORONIC ACID
#sed -i 's/C1X/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    C1X LYS  MODIFIED LYSINE
#sed -i 's/CLG/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CLG LYS
#sed -i 's/CLH/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CLH LYS
#sed -i 's/CYJ/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    CYJ LYS  MODIFIED LYSINE
#sed -i 's/DLS/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    DLS LYS  DI-ACETYL-LYSINE
#sed -i 's/DLY/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    DLY LYS  D-LYSINE
#sed -i 's/DNL/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    DNL LYS  6-AMINO-HEXANAL
#sed -i 's/FHL/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    FHL LYS  MODIFIED LYSINE
#sed -i 's/GPL/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    GPL LYS  LYSINE GUANOSINE-5-MONOPHOSPHATE
#sed -i 's/IT1/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    IT1 LYS
#sed -i 's/KCX/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    KCX LYS  CARBAMOYLATED LYSINE
#sed -i 's/KGC/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    KGC LYS
#sed -i 's/KST/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    KST LYS  N~6~-(5-CARBOXY-3-THIENYL)-L-LYSINE
#sed -i 's/LA2/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    LA2 LYS
#sed -i 's/LCK/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    LCK LYS
#sed -i 's/LCX/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    LCX LYS  CARBAMYLATED LYSINE
#sed -i 's/LDH/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    LDH LYS  N~6~-ETHYL-L-LYSINE
#sed -i 's/LET/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    LET LYS  ODIFIED LYSINE
#sed -i 's/LLP/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    LLP LYS
#sed -i 's/LLY/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    LLY LYS  NZ-(DICARBOXYMETHYL)LYSINE
#sed -i 's/LSO/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    LSO LYS  MODIFIED LYSINE
#sed -i 's/LYM/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    LYM LYS  DEOXY-METHYL-LYSINE
#sed -i 's/LYN/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    LYN LYS  26-DIAMINO-HEXANOIC ACID AMIDE
#sed -i 's/LYP/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    LYP LYS  N~6~-METHYL-N~6~-PROPYL-L-LYSINE
#sed -i 's/LYR/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    LYR LYS  MODIFIED LYSINE
#sed -i 's/LYS/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    LYS LYS
#sed -i 's/LYX/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    LYX LYS  N-(2-COENZYME A)-PROPANOYL-LYSINE
#sed -i 's/LYZ/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    LYZ LYS  5-HYDROXYLYSINE
#sed -i 's/M2L/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    M2L LYS
#sed -i 's/M3L/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    M3L LYS  N-TRIMETHYLLYSINE
#sed -i 's/MCL/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    MCL LYS  NZ-(1-CARBOXYETHYL)-LYSINE
#sed -i 's/MLY/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    MLY LYS  METHYLATED LYSINE
#sed -i 's/MLZ/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    MLZ LYS  N-METHYL-LYSINE
#sed -i 's/OBS/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    OBS LYS  MODIFIED LYSINE
#sed -i 's/SLZ/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    SLZ LYS  L-THIALYSINE
#sed -i 's/XX1/LYS/g' ${collective%.*}_r_prep1.pdb                 ###                    XX1 LYS  N~6~-7H-PURIN-6-YL-L-LYSINE
#sed -i 's/BUG/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/CLE/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/DLE/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/LEU/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                  LEU
#sed -i 's/MLE/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/NLE/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/NLN/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/NLP/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/1LU/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    1LU LEU  4-METHYL-PENTANOIC ACID-2-OXYL GROUP
#sed -i 's/2ML/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    2ML LEU  2-METHYLLEUCINE
#sed -i 's/BLE/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    BLE LEU  LEUCINE BORONIC ACID
#sed -i 's/BUG/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    BUG LEU  TERT-LEUCYL AMINE
#sed -i 's/CLE/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    CLE LEU  LEUCINE AMIDE
#sed -i 's/DCL/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    DCL LEU  2-AMINO-4-METHYL-PENTANYL GROUP
#sed -i 's/DLE/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    DLE LEU  D-LEUCINE
#sed -i 's/DNE/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    DNE LEU  D-NORLEUCINE
#sed -i 's/DNG/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    DNG LEU  N-FORMYL-D-NORLEUCINE
#sed -i 's/DNM/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    DNM LEU  D-N-METHYL NORLEUCINE
#sed -i 's/FLE/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    FLE LEU  FUROYL-LEUCINE
#sed -i 's/HLU/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    HLU LEU  BETA-HYDROXYLEUCINE
#sed -i 's/LED/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    LED LEU  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/LEF/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    LEF LEU  2-5-FLUOROLEUCINE
#sed -i 's/LEU/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    LEU LEU
#sed -i 's/LNT/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    LNT LEU
#sed -i 's/MHL/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    MHL LEU  N-METHYLATED HYDROXY
#sed -i 's/MLE/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    MLE LEU  N-METHYLATED
#sed -i 's/MLL/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    MLL LEU  METHYL L-LEUCINATE
#sed -i 's/MNL/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    MNL LEU  4N-DIMETHYLNORLEUCINE
#sed -i 's/NLE/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    NLE LEU  NORLEUCINE
#sed -i 's/NLN/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    NLN LEU  NORLEUCINE AMIDE
#sed -i 's/NLO/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    NLO LEU  O-METHYL-L-NORLEUCINE
#sed -i 's/PLE/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    PLE LEU  LEUCINE PHOSPHINIC ACID
#sed -i 's/PPH/LEU/g' ${collective%.*}_r_prep1.pdb                 ###                    PPH LEU  PHENYLALANINE PHOSPHINIC ACID
#sed -i 's/CXM/MET/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
#sed -i 's/FME/MET/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
#sed -i 's/MET/MET/g' ${collective%.*}_r_prep1.pdb                 ###                  MET
#sed -i 's/MSE/MET/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
#sed -i 's/OMT/MET/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
#sed -i 's/AME/MET/g' ${collective%.*}_r_prep1.pdb                 ###                    AME MET  ACETYLATED METHIONINE
#sed -i 's/CXM/MET/g' ${collective%.*}_r_prep1.pdb                 ###                    CXM MET  N-CARBOXYMETHIONINE
#sed -i 's/ESC/MET/g' ${collective%.*}_r_prep1.pdb                 ###                    ESC MET  2-AMINO-4-ETHYL SULFANYL BUTYRIC ACID
#sed -i 's/FME/MET/g' ${collective%.*}_r_prep1.pdb                 ###                    FME MET  FORMYL-METHIONINE
#sed -i 's/FOR/MET/g' ${collective%.*}_r_prep1.pdb                 ###                    FOR MET
#sed -i 's/MET/MET/g' ${collective%.*}_r_prep1.pdb                 ###                    MET MET
#sed -i 's/MHO/MET/g' ${collective%.*}_r_prep1.pdb                 ###                    MHO MET  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/MME/MET/g' ${collective%.*}_r_prep1.pdb                 ###                    MME MET  N-METHYL METHIONINE
#sed -i 's/MSE/MET/g' ${collective%.*}_r_prep1.pdb                 ###                    MSE MET  ELENOMETHIONINE
#sed -i 's/MSO/MET/g' ${collective%.*}_r_prep1.pdb                 ###                    MSO MET  METHIONINE SULFOXIDE
#sed -i 's/OMT/MET/g' ${collective%.*}_r_prep1.pdb                 ###                    OMT MET  METHIONINE SULFONE
#sed -i 's/SME/MET/g' ${collective%.*}_r_prep1.pdb                 ###                    SME MET  METHIONINE SULFOXIDE
#sed -i 's/ASN/ASN/g' ${collective%.*}_r_prep1.pdb                 ###                  ASN
#sed -i 's/MEN/ASN/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASN
#sed -i 's/AFA/ASN/g' ${collective%.*}_r_prep1.pdb                 ###                    AFA ASN  N-[7-METHYL-OCT-24-DIENOYL]ASPARAGINE
#sed -i 's/AHB/ASN/g' ${collective%.*}_r_prep1.pdb                 ###                    AHB ASN  BETA-HYDROXYASPARAGINE
#sed -i 's/ASN/ASN/g' ${collective%.*}_r_prep1.pdb                 ###                    ASN ASN
#sed -i 's/B3X/ASN/g' ${collective%.*}_r_prep1.pdb                 ###                    B3X ASN  (3S)-35-DIAMINO-5-OXOPENTANOIC ACID
#sed -i 's/DMH/ASN/g' ${collective%.*}_r_prep1.pdb                 ###                    DMH ASN  N4N4-DIMETHYL-ASPARAGINE
#sed -i 's/DSG/ASN/g' ${collective%.*}_r_prep1.pdb                 ###                    DSG ASN  D-ASPARAGINE
#sed -i 's/MEN/ASN/g' ${collective%.*}_r_prep1.pdb                 ###                    MEN ASN  GAMMA METHYL ASPARAGINE
#sed -i 's/DPR/PRO/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PRO
#sed -i 's/PRO/PRO/g' ${collective%.*}_r_prep1.pdb                 ###                  PRO
#sed -i 's/1AB/PRO/g' ${collective%.*}_r_prep1.pdb                 ###                    1AB PRO  14-DIDEOXY-14-IMINO-D-ARABINITOL
#sed -i 's/2MT/PRO/g' ${collective%.*}_r_prep1.pdb                 ###                    2MT PRO
#sed -i 's/4FB/PRO/g' ${collective%.*}_r_prep1.pdb                 ###                    4FB PRO  (4S)-4-FLUORO-L-PROLINE
#sed -i 's/DPL/PRO/g' ${collective%.*}_r_prep1.pdb                 ###                    DPL PRO  4-OXOPROLINE
#sed -i 's/DPR/PRO/g' ${collective%.*}_r_prep1.pdb                 ###                    DPR PRO  D-PROLINE
#sed -i 's/H5M/PRO/g' ${collective%.*}_r_prep1.pdb                 ###                    H5M PRO  TRANS-3-HYDROXY-5-METHYLPROLINE
#sed -i 's/HY3/PRO/g' ${collective%.*}_r_prep1.pdb                 ###                    HY3 PRO  3-HYDROXYPROLINE
#sed -i 's/HYP/PRO/g' ${collective%.*}_r_prep1.pdb                 ###                    HYP PRO  4-HYDROXYPROLINE
#sed -i 's/LPD/PRO/g' ${collective%.*}_r_prep1.pdb                 ###                    LPD PRO  L-PROLINAMIDE
#sed -i 's/P2Y/PRO/g' ${collective%.*}_r_prep1.pdb                 ###                    P2Y PRO  (2S)-PYRROLIDIN-2-YLMETHYLAMINE
#sed -i 's/PCA/PRO/g' ${collective%.*}_r_prep1.pdb                 ###                    PCA PRO  5-OXOPROLINE
#sed -i 's/POM/PRO/g' ${collective%.*}_r_prep1.pdb                 ###                    POM PRO  CIS-5-METHYL-4-OXOPROLINE
#sed -i 's/PRO/PRO/g' ${collective%.*}_r_prep1.pdb                 ###                    PRO PRO
#sed -i 's/PRS/PRO/g' ${collective%.*}_r_prep1.pdb                 ###                    PRS PRO  THIOPROLINE
#sed -i 's/DGN/GLN/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLN
#sed -i 's/GLN/GLN/g' ${collective%.*}_r_prep1.pdb                 ###                  GLN
#sed -i 's/DGN/GLN/g' ${collective%.*}_r_prep1.pdb                 ###                    DGN GLN  D-GLUTAMINE
#sed -i 's/GHG/GLN/g' ${collective%.*}_r_prep1.pdb                 ###                    GHG GLN  GAMMA-HYDROXY-GLUTAMINE
#sed -i 's/GLH/GLN/g' ${collective%.*}_r_prep1.pdb                 ###                    GLH GLN
#sed -i 's/GLN/GLN/g' ${collective%.*}_r_prep1.pdb                 ###                    GLN GLN
#sed -i 's/MGN/GLN/g' ${collective%.*}_r_prep1.pdb                 ###                    MGN GLN  2-METHYL-GLUTAMINE
#sed -i 's/ACL/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/AGM/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/ARG/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                  ARG
#sed -i 's/ARM/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/DAR/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/HAR/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/HMR/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/2MR/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    2MR ARG  N3 N4-DIMETHYLARGININE
#sed -i 's/AAR/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    AAR ARG  ARGININEAMIDE
#sed -i 's/ACL/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    ACL ARG  DEOXY-CHLOROMETHYL-ARGININE
#sed -i 's/AGM/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    AGM ARG  4-METHYL-ARGININE
#sed -i 's/ALG/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    ALG ARG  GUANIDINOBUTYRYL GROUP
#sed -i 's/AR2/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    AR2 ARG  ARGINYL-BENZOTHIAZOLE-6-CARBOXYLIC ACID
#sed -i 's/ARG/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    ARG ARG
#sed -i 's/ARM/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    ARM ARG  DEOXY-METHYL-ARGININE
#sed -i 's/ARO/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    ARO ARG  C-GAMMA-HYDROXY ARGININE
#sed -i 's/BOR/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    BOR ARG
#sed -i 's/CIR/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    CIR ARG  CITRULLINE
#sed -i 's/DA2/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    DA2 ARG  MODIFIED ARGININE
#sed -i 's/DAR/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    DAR ARG  D-ARGININE
#sed -i 's/HMR/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    HMR ARG  BETA-HOMOARGININE
#sed -i 's/HRG/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    HRG ARG  L-HOMOARGININE
#sed -i 's/MAI/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    MAI ARG  DEOXO-METHYLARGININE
#sed -i 's/MGG/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    MGG ARG  MODIFIED D-ARGININE
#sed -i 's/NMM/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    NMM ARG  MODIFIED ARGININE
#sed -i 's/OPR/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    OPR ARG  C-(3-OXOPROPYL)ARGININE
#sed -i 's/ORQ/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    ORQ ARG  N~5~-ACETYL-L-ORNITHINE
#sed -i 's/TYZ/ARG/g' ${collective%.*}_r_prep1.pdb                 ###                    TYZ ARG  PARA ACETAMIDO BENZOIC ACID
#sed -i 's/DSN/SER/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/MIS/SER/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/OAS/SER/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SAC/SER/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SEL/SER/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SEP/SER/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SER/SER/g' ${collective%.*}_r_prep1.pdb                 ###                  SER
#sed -i 's/SET/SER/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SVA/SER/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/B3S/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    B3S SER  (3R)-3-AMINO-4-HYDROXYBUTANOIC ACID
#sed -i 's/BG1/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    BG1 SER
#sed -i 's/DHL/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    DHL SER  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/DSE/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    DSE SER  D-SERINE N-METHYLATED
#sed -i 's/DSN/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    DSN SER  D-SERINE
#sed -i 's/FGP/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    FGP SER
#sed -i 's/GVL/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    GVL SER  SERINE MODIFED WITH PHOSPHOPANTETHEINE
#sed -i 's/HSE/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    HSE SER  L-HOMOSERINE
#sed -i 's/HSL/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    HSL SER  HOMOSERINE LACTONE
#sed -i 's/MC1/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    MC1 SER  METHICILLIN ACYL-SERINE
#sed -i 's/MIS/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    MIS SER  MODIFIED SERINE
#sed -i 's/N10/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    N10 SER  O-[(HEXYLAMINO)CARBONYL]-L-SERINE
#sed -i 's/NC1/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    NC1 SER  NITROCEFIN ACYL-SERINE
#sed -i 's/OAS/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    OAS SER  O-ACETYLSERINE
#sed -i 's/OSE/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    OSE SER  O-SULFO-L-SERINE
#sed -i 's/PG1/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    PG1 SER  BENZYLPENICILLOYL-ACYLATED SERINE
#sed -i 's/PYR/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    PYR SER  CHEMICALLY MODIFIED
#sed -i 's/S1H/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    S1H SER  1-HEXADECANOSULFONYL-O-L-SERINE
#sed -i 's/SAC/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SAC SER  N-ACETYL-SERINE
#sed -i 's/SBD/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SBD SER
#sed -i 's/SBG/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SBG SER  MODIFIED SERINE
#sed -i 's/SBL/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SBL SER
#sed -i 's/SDP/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SDP SER
#sed -i 's/SEB/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SEB SER  O-BENZYLSULFONYL-SERINE
#sed -i 's/SEL/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SEL SER  2-AMINO-13-PROPANEDIOL
#sed -i 's/SEP/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SEP SER  E PHOSPHOSERINE
#sed -i 's/SER/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SER SER
#sed -i 's/SET/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SET SER  AMINOSERINE
#sed -i 's/SGB/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SGB SER  MODIFIED SERINE
#sed -i 's/SGR/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SGR SER  MODIFIED SERINE
#sed -i 's/SOY/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SOY SER  OXACILLOYL-ACYLATED SERINE
#sed -i 's/SUN/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SUN SER  TABUN CONJUGATED SERINE
#sed -i 's/SVA/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SVA SER  SERINE VANADATE
#sed -i 's/SVV/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SVV SER  MODIFIED SERINE
#sed -i 's/SVX/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SVX SER  MODIFIED SERINE
#sed -i 's/SVY/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SVY SER  MODIFIED SERINE
#sed -i 's/SVZ/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SVZ SER  MODIFIED SERINE
#sed -i 's/SXE/SER/g' ${collective%.*}_r_prep1.pdb                 ###                    SXE SER  MODIFIED SERINE
#sed -i 's/ALO/THR/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
#sed -i 's/BMT/THR/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
#sed -i 's/DTH/THR/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
#sed -i 's/THR/THR/g' ${collective%.*}_r_prep1.pdb                 ###                  THR
#sed -i 's/TPO/THR/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
#sed -i 's/AEI/THR/g' ${collective%.*}_r_prep1.pdb                 ###                    AEI THR  ACYLATED THR
#sed -i 's/ALO/THR/g' ${collective%.*}_r_prep1.pdb                 ###                    ALO THR  ALLO-THREONINE
#sed -i 's/BMT/THR/g' ${collective%.*}_r_prep1.pdb                 ###                    BMT THR
#sed -i 's/CRO/THR/g' ${collective%.*}_r_prep1.pdb                 ###                    CRO THR  CYCLIZED
#sed -i 's/CTH/THR/g' ${collective%.*}_r_prep1.pdb                 ###                    CTH THR  4-CHLOROTHREONINE
#sed -i 's/DTH/THR/g' ${collective%.*}_r_prep1.pdb                 ###                    DTH THR  D-THREONINE
#sed -i 's/OLT/THR/g' ${collective%.*}_r_prep1.pdb                 ###                    OLT THR  O-METHYL-L-THREONINE
#sed -i 's/TBM/THR/g' ${collective%.*}_r_prep1.pdb                 ###                    TBM THR
#sed -i 's/TH5/THR/g' ${collective%.*}_r_prep1.pdb                 ###                    TH5 THR  O-ACETYL-L-THREONINE
#sed -i 's/THC/THR/g' ${collective%.*}_r_prep1.pdb                 ###                    THC THR  N-METHYLCARBONYLTHREONINE
#sed -i 's/THR/THR/g' ${collective%.*}_r_prep1.pdb                 ###                    THR THR
#sed -i 's/TMD/THR/g' ${collective%.*}_r_prep1.pdb                 ###                    TMD THR  N-METHYLATED EPSILON C ALKYLATED
#sed -i 's/TPO/THR/g' ${collective%.*}_r_prep1.pdb                 ###                    TPO THR  HOSPHOTHREONINE
#sed -i 's/DIV/VAL/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS VAL
#sed -i 's/DVA/VAL/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS VAL
#sed -i 's/MVA/VAL/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS VAL
#sed -i 's/VAL/VAL/g' ${collective%.*}_r_prep1.pdb                 ###                  VAL
#sed -i 's/B2V/VAL/g' ${collective%.*}_r_prep1.pdb                 ###                    B2V VAL  VALINE BORONIC ACID
#sed -i 's/DIV/VAL/g' ${collective%.*}_r_prep1.pdb                 ###                    DIV VAL  D-ISOVALINE
#sed -i 's/DVA/VAL/g' ${collective%.*}_r_prep1.pdb                 ###                    DVA VAL  D-VALINE
#sed -i 's/MNV/VAL/g' ${collective%.*}_r_prep1.pdb                 ###                    MNV VAL  N-METHYL-C-AMINO VALINE
#sed -i 's/MVA/VAL/g' ${collective%.*}_r_prep1.pdb                 ###                    MVA VAL  N-METHYLATED
#sed -i 's/NVA/VAL/g' ${collective%.*}_r_prep1.pdb                 ###                    NVA VAL  NORVALINE
#sed -i 's/VAD/VAL/g' ${collective%.*}_r_prep1.pdb                 ###                    VAD VAL  DEAMINOHYDROXYVALINE
#sed -i 's/VAF/VAL/g' ${collective%.*}_r_prep1.pdb                 ###                    VAF VAL  METHYLVALINE
#sed -i 's/VAL/VAL/g' ${collective%.*}_r_prep1.pdb                 ###                    VAL VAL
#sed -i 's/VDL/VAL/g' ${collective%.*}_r_prep1.pdb                 ###                    VDL VAL  (2R3R)-23-DIAMINOBUTANOIC ACID
#sed -i 's/VLL/VAL/g' ${collective%.*}_r_prep1.pdb                 ###                    VLL VAL  (2S)-23-DIAMINOBUTANOIC ACID
#sed -i 's/VME/VAL/g' ${collective%.*}_r_prep1.pdb                 ###                    VME VAL  O- METHYLVALINE
#sed -i 's/DTR/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/HTR/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/LTR/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/TPL/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/TRO/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/TRP/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                  TRP
#sed -i 's/BTR/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    BTR TRP  6-BROMO-TRYPTOPHAN
#sed -i 's/1TQ/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    1TQ TRP  6-(FORMYLAMINO)-7-HYDROXY-L-TRYPTOPHAN
#sed -i 's/23S/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    23S TRP  MODIFIED TRYPTOPHAN
#sed -i 's/32S/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    32S TRP  MODIFIED TRYPTOPHAN
#sed -i 's/32T/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    32T TRP  MODIFIED TRYPTOPHAN
#sed -i 's/4DP/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    4DP TRP
#sed -i 's/4FW/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    4FW TRP  4-FLUOROTRYPTOPHANE
#sed -i 's/4HT/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    4HT TRP  4-HYDROXYTRYPTOPHAN
#sed -i 's/4IN/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    4IN TRP  4-AMINO-L-TRYPTOPHAN
#sed -i 's/6CW/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    6CW TRP  6-CHLORO-L-TRYPTOPHAN
#sed -i 's/DTR/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    DTR TRP  D-TRYPTOPHAN
#sed -i 's/FTR/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    FTR TRP  FLUOROTRYPTOPHANE
#sed -i 's/HTR/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    HTR TRP  BETA-HYDROXYTRYPTOPHANE
#sed -i 's/PAT/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    PAT TRP  ALPHA-PHOSPHONO-TRYPTOPHAN
#sed -i 's/TOX/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    TOX TRP
#sed -i 's/TPL/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    TPL TRP  TRYTOPHANOL
#sed -i 's/TQQ/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    TQQ TRP
#sed -i 's/TRF/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    TRF TRP  N1-FORMYL-TRYPTOPHAN
#sed -i 's/TRN/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    TRN TRP  AZA-TRYPTOPHAN
#sed -i 's/TRO/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    TRO TRP  2-HYDROXY-TRYPTOPHAN
#sed -i 's/TRP/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    TRP TRP
#sed -i 's/TRQ/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    TRQ TRP
#sed -i 's/TRW/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    TRW TRP
#sed -i 's/TRX/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    TRX TRP  6-HYDROXYTRYPTOPHAN
#sed -i 's/TTQ/TRP/g' ${collective%.*}_r_prep1.pdb                 ###                    TTQ TRP  6-AMINO-7-HYDROXY-L-TRYPTOPHAN
#sed -i 's/DTY/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/IYR/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/PAQ/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/PTR/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/STY/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/TYB/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/TYQ/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/TYR/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/TYS/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                  TYR
#sed -i 's/TYY/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/1TY/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    1TY TYR
#sed -i 's/2TY/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    2TY TYR
#sed -i 's/3TY/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    3TY TYR  MODIFIED TYROSINE
#sed -i 's/B3Y/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    B3Y TYR
#sed -i 's/CRQ/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    CRQ TYR
#sed -i 's/DBY/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    DBY TYR  35 DIBROMOTYROSINE
#sed -i 's/DPQ/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    DPQ TYR  TYROSINE DERIVATIVE
#sed -i 's/DTY/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    DTY TYR  D-TYROSINE
#sed -i 's/ESB/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    ESB TYR
#sed -i 's/FLT/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    FLT TYR  FLUOROMALONYL TYROSINE
#sed -i 's/FTY/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    FTY TYR  DEOXY-DIFLUOROMETHELENE-PHOSPHOTYROSINE
#sed -i 's/IYR/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    IYR TYR  3-IODO-TYROSINE
#sed -i 's/MBQ/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    MBQ TYR
#sed -i 's/NIY/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    NIY TYR  META-NITRO-TYROSINE
#sed -i 's/NBQ/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    NBQ TYR
#sed -i 's/OTY/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    OTY TYR
#sed -i 's/PAQ/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    PAQ TYR  SEE REMARK 999
#sed -i 's/PTH/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    PTH TYR  METHYLENE-HYDROXY-PHOSPHOTYROSINE
#sed -i 's/PTM/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    PTM TYR  ALPHA-METHYL-O-PHOSPHOTYROSINE
#sed -i 's/PTR/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    PTR TYR  O-PHOSPHOTYROSINE
#sed -i 's/TCQ/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    TCQ TYR  MODIFIED TYROSINE
#sed -i 's/TTS/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    TTS TYR
#sed -i 's/TY2/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    TY2 TYR  3-AMINO-L-TYROSINE
#sed -i 's/TY3/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    TY3 TYR  3-HYDROXY-L-TYROSINE
#sed -i 's/TYB/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    TYB TYR  TYROSINAL
#sed -i 's/TYC/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    TYC TYR  L-TYROSINAMIDE
#sed -i 's/TYI/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    TYI TYR  35-DIIODOTYROSINE
#sed -i 's/TYN/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    TYN TYR  ADDUCT AT HYDROXY GROUP
#sed -i 's/TYO/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    TYO TYR
#sed -i 's/TYQ/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    TYQ TYR  AMINOQUINOL FORM OF TOPA QUINONONE
#sed -i 's/TYR/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    TYR TYR
#sed -i 's/TYS/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    TYS TYR  INE SULPHONATED TYROSINE
#sed -i 's/TYT/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    TYT TYR
#sed -i 's/TYY/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    TYY TYR  IMINOQUINONE FORM OF TOPA QUINONONE
#sed -i 's/YOF/TYR/g' ${collective%.*}_r_prep1.pdb                 ###                    YOF TYR  3-FLUOROTYROSINE


#################################################################################################################################################################



pdb4amber -i ${collective%.*}_r_prep1.pdb -o ${collective%.*}_r_prep11.pdb
sed -i '/HETATM.*ZN\|ZN.*HETATM/d' ${collective%.*}_ligand.pdb
sed -i '/HETATM.*NI\|NI.*HETATM/d' ${collective%.*}_ligand.pdb
sed -i '/HETATM.*MG\|MG.*HETATM/d' ${collective%.*}_ligand.pdb
sed -i '/HETATM.*CA\|CA.*HETATM/d' ${collective%.*}_ligand.pdb
sed -i '/HETATM.*MN\|MN.*HETATM/d' ${collective%.*}_ligand.pdb

sed -i 's/ATOM  /HETATM/g' ${collective%.*}_ligand.pdb
sed -i 's/CL/Cl/g' ${collective%.*}_ligand.pdb
sed -i 's/BR/Br/g' ${collective%.*}_ligand.pdb

#### EDITED
antechamber -fi pdb -i ${collective%.*}_ligand.pdb -fo mol2 -o ${collective%.*}_ligand_oprot_prep1_gas.mol2 -c gas -pf y

antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc 0 -pf y
else
	:
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -1 -pf y
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +1 -pf y
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -2 -pf y
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +2 -pf y
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -3 -pf y
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +3 -pf y
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +4 -pf y
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -4 -pf y
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -5 -pf y
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +5 -pf y
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
	#pdb4amber -i ${collective%.*}_ligand.pdb -o ${collective%.*}_ligand_oprot_prep1.pdb --reduce
	#antechamber -fi pdb -i ${collective%.*}_ligand_oprot_prep1.pdb -fo mol2 -o ${collective%.*}_ligand_oprot_prep1_gas.mol2 -c gas -pf y
	antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -dr n
else
	:
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc 0 -pf y -dr n
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -1 -pf y -dr n
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +1 -pf y -dr n
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -2 -pf y -dr n
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +2 -pf y -dr n
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -3 -pf y -dr n
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +3 -pf y -dr n
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +4 -pf y -dr n
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -4 -pf y -dr n
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc -5 -pf y -dr n
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi mol2 -i ${collective%.*}_ligand_oprot_prep1_gas.mol2 -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -nc +5 -pf y -dr n
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas.mol2 ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas.mol2 ]; then
        antechamber -fi pdb -i ${collective%.*}_ligand.pdb -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${collective%.*}_ligand.pdb -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 0 -dr n
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${collective%.*}_ligand.pdb -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 1 -dr n
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${collective%.*}_ligand.pdb -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -1 -dr n
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${collective%.*}_ligand.pdb -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 2 -dr n
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${collective%.*}_ligand.pdb -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -2 -dr n
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${collective%.*}_ligand.pdb -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc 3 -dr n
else
        :
fi

if [ ! -f ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ] || [ ! -s ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi ]; then
        antechamber -fi pdb -i ${collective%.*}_ligand.pdb -fo prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -c bcc -pf y -nc -3 -dr n
else
        :
fi


parmchk2 -f prepi -i ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi -o ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${collective%.*}_r_prep11.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${collective%.*}_r_prep2.prmtop ${collective%.*}_r_prep2.inpcrd
quit
@


cp ${collective%.*}_r_prep11.pdb ${collective%.*}_r_prep00012.pdb


declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 5);
do
if [ ! -f ${collective%.*}_r_prep2.inpcrd ] || [ ! -s ${collective%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Paramter expansion susbtring extraction ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis000${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis000${count}.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis000${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
a=${A:1:3}
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${collective%.*}_r_prep009.pdb >> ${collective%.*}_r_prep010.pdb
sed -e "/${a}.*${B}.*H/d" ${collective%.*}_r_prep000${pcount}.pdb >> ${collective%.*}_r_prep000${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${collective%.*}_r_prep000${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${collective%.*}_r_prep2.prmtop ${collective%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;

cp ${collective%.*}_r_prep00017.pdb ${collective%.*}_r_prep1112.pdb

declare -i pcount=12
declare -i p2count=13

for i in $(seq 35);
do

if [ ! -f ${collective%.*}_r_prep2.inpcrd ] || [ ! -s ${collective%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis11${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis11${count}.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis11${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

sed -e "/${A}.*${B}.*H/d" ${collective%.*}_r_prep11${pcount}.pdb >> ${collective%.*}_r_prep11${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${collective%.*}_r_prep11${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${collective%.*}_r_prep2.prmtop ${collective%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :

fi;
done;
cp ${collective%.*}_r_prep1145.pdb ${collective%.*}_r_prep11112.pdb

declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 5);
do
if [ ! -f ${collective%.*}_r_prep2.inpcrd ] || [ ! -s ${collective%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Paramter expansion susbtring extraction ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis111${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis111${count}.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis111${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
a=${A:1:3}
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${collective%.*}_r_prep009.pdb >> ${collective%.*}_r_prep010.pdb
sed -e "/${a}.*${B}.*H/d" ${collective%.*}_r_prep111${pcount}.pdb >> ${collective%.*}_r_prep111${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${collective%.*}_r_prep111${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${collective%.*}_r_prep2.prmtop ${collective%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;
cp ${collective%.*}_r_prep11117.pdb ${collective%.*}_r_prep0012.pdb

declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 35);
do
if [ ! -f ${collective%.*}_r_prep2.inpcrd ] || [ ! -s ${collective%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - BACKBONE removal in progress ...";
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis00${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis00${count}.out);
del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis00${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"
C="$(cut -d' ' -f1 <<<$del2)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${collective%.*}_r_prep00${pcount}.pdb >> ${collective%.*}_r_prep00${p2count}.pdb
sed -e "/${C}.*${A}.*${B}/d" ${collective%.*}_r_prep00${pcount}.pdb >> ${collective%.*}_r_prep00${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${collective%.*}_r_prep00${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${collective%.*}_r_prep2.prmtop ${collective%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;


if [ ! -f ${collective%.*}_r_prep2.inpcrd ] || [ ! -s ${collective%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - O1P removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis001.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis001.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis001.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${collective%.*}_r_prep0047.pdb >> ${collective%.*}_r_prep001.pdb
sed -e "/O1P.*${A}.*${B}/d" ${collective%.*}_r_prep0047.pdb >> ${collective%.*}_r_prep001.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${collective%.*}_r_prep001.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${collective%.*}_r_prep2.prmtop ${collective%.*}_r_prep2.inpcrd
quit
@
fi;


if [ ! -f ${collective%.*}_r_prep2.inpcrd ] || [ ! -s ${collective%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - O2P removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis003.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis003.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis003.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${collective%.*}_r_prep12.pdb >> ${collective%.*}_r_prep002.pdb
sed -e "/O2P.*${A}.*${B}/d" ${collective%.*}_r_prep001.pdb >> ${collective%.*}_r_prep002.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${collective%.*}_r_prep002.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${collective%.*}_r_prep2.prmtop ${collective%.*}_r_prep2.inpcrd
quit
@
fi;



if [ ! -f ${collective%.*}_r_prep2.inpcrd ] || [ ! -s ${collective%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - O3P removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis005.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis005.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis005.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${collective%.*}_r_prep12.pdb >> ${collective%.*}_r_prep003.pdb
sed -e "/O3P.*${A}.*${B}/d" ${collective%.*}_r_prep002.pdb >> ${collective%.*}_r_prep003.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${collective%.*}_r_prep003.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${collective%.*}_r_prep2.prmtop ${collective%.*}_r_prep2.inpcrd
quit
@
fi;



if [ ! -f ${collective%.*}_r_prep2.inpcrd ] || [ ! -s ${collective%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis007.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis007.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis007.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${collective%.*}_r_prep003.pdb >> ${collective%.*}_r_prep004.pdb
sed -e "/${A}.*${B}.*P/d" ${collective%.*}_r_prep003.pdb >> ${collective%.*}_r_prep004.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${collective%.*}_r_prep004.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${collective%.*}_r_prep2.prmtop ${collective%.*}_r_prep2.inpcrd
quit
@
fi;


if [ ! -f ${collective%.*}_r_prep2.inpcrd ] || [ ! -s ${collective%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis009.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis009.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis009.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${collective%.*}_r_prep004.pdb >> ${collective%.*}_r_prep005.pdb
sed -e "/${A}.*${B}.*P/d" ${collective%.*}_r_prep004.pdb >> ${collective%.*}_r_prep005.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${collective%.*}_r_prep005.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${collective%.*}_r_prep2.prmtop ${collective%.*}_r_prep2.inpcrd
quit
@
fi;


if [ ! -f ${collective%.*}_r_prep2.inpcrd ] || [ ! -s ${collective%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis011.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis011.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis011.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${collective%.*}_r_prep005.pdb >> ${collective%.*}_r_prep006.pdb
sed -e "/${A}.*${B}.*P/d" ${collective%.*}_r_prep005.pdb >> ${collective%.*}_r_prep006.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${collective%.*}_r_prep006.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${collective%.*}_r_prep2.prmtop ${collective%.*}_r_prep2.inpcrd
quit
@
fi;


if [ ! -f ${collective%.*}_r_prep2.inpcrd ] || [ ! -s ${collective%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis013.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis013.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis013.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${collective%.*}_r_prep006.pdb >> ${collective%.*}_r_prep007.pdb
sed -e "/${A}.*${B}.*H/d" ${collective%.*}_r_prep006.pdb >> ${collective%.*}_r_prep007.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${collective%.*}_r_prep007.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${collective%.*}_r_prep2.prmtop ${collective%.*}_r_prep2.inpcrd
quit
@
fi;


if [ ! -f ${collective%.*}_r_prep2.inpcrd ] || [ ! -s ${collective%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis015.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis015.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis015.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${collective%.*}_r_prep007.pdb >> ${collective%.*}_r_prep008.pdb
sed -e "/${A}.*${B}.*H/d" ${collective%.*}_r_prep007.pdb >> ${collective%.*}_r_prep008.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${collective%.*}_r_prep008.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${collective%.*}_r_prep2.prmtop ${collective%.*}_r_prep2.inpcrd
quit
@
fi;


if [ ! -f ${collective%.*}_r_prep2.inpcrd ] || [ ! -s ${collective%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis017.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis017.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis017.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${collective%.*}_r_prep008.pdb >> ${collective%.*}_r_prep009.pdb
sed -e "/${A}.*${B}.*H/d" ${collective%.*}_r_prep008.pdb >> ${collective%.*}_r_prep009.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${collective%.*}_r_prep009.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${collective%.*}_r_prep2.prmtop ${collective%.*}_r_prep2.inpcrd
quit
@
fi;


if [ ! -f ${collective%.*}_r_prep2.inpcrd ] || [ ! -s ${collective%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Paramter expansion susbtring extraction ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis019.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis019.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis019.out);
A="$(cut -d' ' -f1 <<<$del)"
a=${A:1:3}
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${collective%.*}_r_prep009.pdb >> ${collective%.*}_r_prep010.pdb
sed -e "/${a}.*${B}.*H/d" ${collective%.*}_r_prep009.pdb >> ${collective%.*}_r_prep010.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${collective%.*}_r_prep010.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${collective%.*}_r_prep2.prmtop ${collective%.*}_r_prep2.inpcrd
quit
@
fi;


if [ ! -f ${collective%.*}_r_prep2.inpcrd ] || [ ! -s ${collective%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - HIS removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis021.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis021.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis021.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${collective%.*}_r_prep010.pdb >> ${collective%.*}_r_prep011.pdb
sed -e "/HIS.*${B}.*H/d" ${collective%.*}_r_prep010.pdb >> ${collective%.*}_r_prep011.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${collective%.*}_r_prep011.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${collective%.*}_r_prep2.prmtop ${collective%.*}_r_prep2.inpcrd
quit
@
fi;


declare -i count=23
declare -i pcount=11
declare -i p2count=12

for i in $(seq 35);
do
if [ ! -f ${collective%.*}_r_prep2.inpcrd ] || [ ! -s ${collective%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - BACKBONE removal in progress ...";
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis0${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis0${count}.out);
del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis0${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"
C="$(cut -d' ' -f1 <<<$del2)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${collective%.*}_r_prep0${pcount}.pdb >> ${collective%.*}_r_prep0${p2count}.pdb
sed -e "/${C}.*${A}.*${B}/d" ${collective%.*}_r_prep0${pcount}.pdb >> ${collective%.*}_r_prep0${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${collective%.*}_r_prep0${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${collective%.*}_r_prep2.prmtop ${collective%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
	:
fi;
done;

if [ ! -f ${collective%.*}_r_prep2.inpcrd ] || [ ! -s ${collective%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P, O1P, O2P, O3P removal in progress ..."

sed -i '/ATOM.*\<P\>/d' ${collective%.*}_r_prep046.pdb
sed -i '/ATOM.*\<O1P\>/d' ${collective%.*}_r_prep046.pdb
sed -i '/ATOM.*\<O2P\>/d' ${collective%.*}_r_prep046.pdb
sed -i '/ATOM.*\<O3P\>/d' ${collective%.*}_r_prep046.pdb



tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${collective%.*}_r_prep046.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${collective%.*}_r_prep2.prmtop ${collective%.*}_r_prep2.inpcrd
quit
@
fi;

declare -i count=1
declare -i pcount=46
declare -i p2count=47

for i in $(seq 8);
do
if [ ! -f ${collective%.*}_r_prep2.inpcrd ] || [ ! -s ${collective%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Paramter expansion susbtring extraction ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis22${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis22${count}.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis22${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
a=${A:1:3}
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${collective%.*}_r_prep009.pdb >> ${collective%.*}_r_prep010.pdb
sed -e "/${a}.*${B}.*H/d" ${collective%.*}_r_prep0${pcount}.pdb >> ${collective%.*}_r_prep0${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${collective%.*}_r_prep0${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${collective%.*}_r_prep2.prmtop ${collective%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;



ambpdb -p ${collective%.*}_r_prep2.prmtop -c ${collective%.*}_r_prep2.inpcrd -mol2 -sybyl > ${collective%.*}_r_prep3.mol2

for model in model*pdb; do


rm -f *trial* *prep*pdb *ligand* sqm.* ANTECHAMBER* ATOMTYPE* *log *noHET* *sed* *sslink* *txt *diagnosis

cp $model ${collective%.*}_prep1.pdb

pdb4amber -i ${model%.*}_prep1.pdb -o ${model%.*}_r_prep1.pdb

sed -i '/ACE/d' ${model%.*}_r_prep1.pdb
sed -i '/NME/d' ${model%.*}_r_prep1.pdb

sed -i '/HETATM/s/CL/Cl/g' ${model%.*}_r_prep1.pdb
sed -i '/HETATM/s/BR/Br/g' ${model%.*}_r_prep1.pdb

############################################################################################################################################################
sed -i 's/CSD/CYS/g' ${model%.*}_r_prep1.pdb                 #3-SULFINOALANINE
sed -i 's/HYP/PRO/g' ${model%.*}_r_prep1.pdb                 #4-HYDROXYPROLINE
sed -i 's/BMT/THR/g' ${model%.*}_r_prep1.pdb                 #4-METHYL-4-[(E)-2-BUTENYL]-4,N-METHYL-THREONINE
sed -i 's/5HP/GLU/g' ${model%.*}_r_prep1.pdb                 #5-HYDROXYPROLINE
sed -i 's/ABA/ALA/g' ${model%.*}_r_prep1.pdb                 #ALPHA-AMINOBUTYRIC_ACID
sed -i 's/AIB/ALA/g' ${model%.*}_r_prep1.pdb                 #ALPHA-AMINOISOBUTYRIC_ACID
sed -i 's/CSW/CYS/g' ${model%.*}_r_prep1.pdb                 #CYSTEINE-S-DIOXIDE
sed -i 's/OCS/CYS/g' ${model%.*}_r_prep1.pdb                 #CYSTEINESULFONIC_ACID
sed -i 's/DAL/ALA/g' ${model%.*}_r_prep1.pdb                 #D-ALANINE
sed -i 's/DAR/ARG/g' ${model%.*}_r_prep1.pdb                 #D-ARGININE
sed -i 's/DSG/ASN/g' ${model%.*}_r_prep1.pdb                 #D-ASPARAGINE
sed -i 's/DSP/ASP/g' ${model%.*}_r_prep1.pdb                 #D-ASPARTATE
sed -i 's/DCY/CYS/g' ${model%.*}_r_prep1.pdb                 #D-CYSTEINE
sed -i 's/CRO/CRO/g' ${model%.*}_r_prep1.pdb                 #DECARBOXY(PARAHYDROXYBENZYLIDENE-IMIDAZOLIDINONE)THREONINE
sed -i 's/DGL/GLU/g' ${model%.*}_r_prep1.pdb                 #D-GLUTAMATE
sed -i 's/DGN/GLN/g' ${model%.*}_r_prep1.pdb                 #D-GLUTAMINE
sed -i 's/DHI/HIS/g' ${model%.*}_r_prep1.pdb                 #D-HISTIDINE
sed -i 's/DIL/ILE/g' ${model%.*}_r_prep1.pdb                 #D-ISOLEUCINE
sed -i 's/DIV/VAL/g' ${model%.*}_r_prep1.pdb                 #D-ISOVALINE
sed -i 's/DLE/LEU/g' ${model%.*}_r_prep1.pdb                 #D-LEUCINE
sed -i 's/DLY/LYS/g' ${model%.*}_r_prep1.pdb                 #D-LYSINE
sed -i 's/DPN/PHE/g' ${model%.*}_r_prep1.pdb                 #D-PHENYLALANINE
sed -i 's/DPR/PRO/g' ${model%.*}_r_prep1.pdb                 #D-PROLINE
sed -i 's/DSN/SER/g' ${model%.*}_r_prep1.pdb                 #D-SERINE
sed -i 's/DTH/THR/g' ${model%.*}_r_prep1.pdb                 #D-THREONINE
sed -i 's/DTR/DTR/g' ${model%.*}_r_prep1.pdb                 #D-TRYPTOPHANE
sed -i 's/DTY/TYR/g' ${model%.*}_r_prep1.pdb                 #D-TYROSINE
sed -i 's/DVA/VAL/g' ${model%.*}_r_prep1.pdb                 #D-VALINE
sed -i 's/CGU/GLU/g' ${model%.*}_r_prep1.pdb                 #GAMMA-CARBOXY-GLUTAMIC_ACID
sed -i 's/KCX/LYS/g' ${model%.*}_r_prep1.pdb                 #LYSINE_NZ-CARBOXYLIC_ACID
sed -i 's/LLP/LYS/g' ${model%.*}_r_prep1.pdb                 #LYSINE-PYRIDOXAL-5'-PHOSPHATE
sed -i 's/CXM/MET/g' ${model%.*}_r_prep1.pdb                 #N-CARBOXYMETHIONINE
sed -i 's/FME/MET/g' ${model%.*}_r_prep1.pdb                 #N-FORMYLMETHIONINE
sed -i 's/MLE/LEU/g' ${model%.*}_r_prep1.pdb                 #N-METHYLLEUCINE
sed -i 's/MVA/VAL/g' ${model%.*}_r_prep1.pdb                 #N-METHYLVALINE
sed -i 's/NLE/LEU/g' ${model%.*}_r_prep1.pdb                 #NORLEUCINE
sed -i 's/PTR/TYR/g' ${model%.*}_r_prep1.pdb                 #O-PHOSPHOTYROSINE
sed -i 's/ORN/ALA/g' ${model%.*}_r_prep1.pdb                 #ORNITHINE
sed -i 's/SEP/SER/g' ${model%.*}_r_prep1.pdb                 #PHOSPHOSERINE
sed -i 's/TPO/THR/g' ${model%.*}_r_prep1.pdb                 #PHOSPHOTHREONINE
sed -i 's/PCA/GLU/g' ${model%.*}_r_prep1.pdb                 #PYROGLUTAMIC_ACID
sed -i 's/SAR/GLY/g' ${model%.*}_r_prep1.pdb                 #SARCOSINE
sed -i 's/CEA/CYS/g' ${model%.*}_r_prep1.pdb                 #S-HYDROXY-CYSTEINE
sed -i 's/CSO/CYS/g' ${model%.*}_r_prep1.pdb                 #S-HYDROXYCYSTEINE
sed -i 's/CSS/CYS/g' ${model%.*}_r_prep1.pdb                 #S-MERCAPTOCYSTEINE
sed -i 's/CSX/CYS/g' ${model%.*}_r_prep1.pdb                 #S-OXY_CYSTEINE
sed -i 's/CME/CYS/g' ${model%.*}_r_prep1.pdb                 #S,S-(2-HYDROXYETHYL)THIOCYSTEINE
sed -i 's/TYS/TYR/g' ${model%.*}_r_prep1.pdb                 #SULFONATED_TYROSINE
sed -i 's/TPQ/PHE/g' ${model%.*}_r_prep1.pdb                 #TOPO-QUINONE
sed -i 's/STY/TYR/g' ${model%.*}_r_prep1.pdb                 #TYROSINE-O-SULPHONIC_ACID
sed -i 's/CCS/CYS/g' ${model%.*}_r_prep1.pdb                 #https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py
sed -i 's/CALA/ALA/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CARG/ARG/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CASN/ASN/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CASP/ASP/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CCYS/CYS/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CCYX/CYX/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CGLN/GLN/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CGLU/GLU/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CGLY/GLY/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CHID/HID/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CHIE/HIE/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CHIP/HIP/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CHYP/HYP/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CILE/ILE/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CLEU/LEU/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CLYS/LYS/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CMET/MET/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CPHE/PHE/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CPRO/PRO/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CSER/SER/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CTHR/THR/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CTRP/TRP/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CTYR/TYR/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/CVAL/VAL/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NALA/ALA/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NARG/ARG/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NASN/ASN/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NASP/ASP/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NCYS/CYS/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NCYX/CYX/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NGLN/GLN/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NGLU/GLU/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NGLY/GLY/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NHID/HID/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NHIE/HIE/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NHIP/HIP/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NILE/ILE/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NLEU/LEU/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NLYS/LYS/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NMET/MET/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NPHE/PHE/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NPRO/PRO/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NSER/SER/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NTHR/THR/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NTRP/TRP/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NTYR/TYR/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/NVAL/VAL/g' ${model%.*}_r_prep1.pdb                #GUESS http://archive.ambermd.org/201909/0130.html
sed -i 's/DAS/ASP/g' ${model%.*}_r_prep1.pdb                 #GUESS https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py
sed -i 's/CAF/CYS/g' ${model%.*}_r_prep1.pdb                 #GUESS https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py ## S-DIMETHYLARSINOYL-CYSTEINE
sed -i 's/CAS/CYS/g' ${model%.*}_r_prep1.pdb                 #GUESS https://github.com/EvanBaugh/miscellaneous_scripts/blob/master/process_pdb.py ## S-(DIMETHYLARSENIC)CYSTEINE
#sed -i 's/AIB/ALA/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/ALA/ALA/g' ${model%.*}_r_prep1.pdb                 ###                  ALA
#sed -i 's/ALM/ALA/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/AYA/ALA/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/BNN/ALA/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/CHG/ALA/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/CSD/ALA/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/DAL/ALA/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/DHA/ALA/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/DNP/ALA/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/FLA/ALA/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/HAC/ALA/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/PRR/ALA/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/MAA/ALA/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/TIH/ALA/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/TPQ/ALA/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ALA
#sed -i 's/0CS/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    0CS ALA  3-[(S)-HYDROPEROXYSULFINYL]-L-ALANINE
#sed -i 's/2BU/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    2BU ADE
#sed -i 's/2OP/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    2OP (2S  2-HYDROXYPROPANAL
#sed -i 's/4F3/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    4F3 ALA  CYCLIZED
#sed -i 's/AA4/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    AA4 ALA  2-AMINO-5-HYDROXYPENTANOIC ACID
#sed -i 's/ABA/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    ABA ALA  ALPHA-AMINOBUTYRIC ACID
#sed -i 's/AHO/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    AHO ALA  N-ACETYL-N-HYDROXY-L-ORNITHINE
#sed -i 's/AHP/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    AHP ALA  2-AMINO-HEPTANOIC ACID
#sed -i 's/AIB/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    AIB ALA  ALPHA-AMINOISOBUTYRIC ACID
#sed -i 's/ALA/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    ALA ALA
#sed -i 's/ALC/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    ALC ALA  2-AMINO-3-CYCLOHEXYL-PROPIONIC ACID
#sed -i 's/ALM/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    ALM ALA  1-METHYL-ALANINAL
#sed -i 's/ALN/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    ALN ALA  NAPHTHALEN-2-YL-3-ALANINE
#sed -i 's/ALS/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    ALS ALA  2-AMINO-3-OXO-4-SULFO-BUTYRIC ACID
#sed -i 's/ALT/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    ALT ALA  THIOALANINE
#sed -i 's/AP7/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    AP7 ADE
#sed -i 's/APH/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    APH ALA  P-AMIDINOPHENYL-3-ALANINE
#sed -i 's/AYA/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    AYA ALA  N-ACETYLALANINE
#sed -i 's/AYG/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    AYG ALA
#sed -i 's/B2A/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    B2A ALA  ALANINE BORONIC ACID
#sed -i 's/B3A/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    B3A ALA  (3S)-3-AMINOBUTANOIC ACID
#sed -i 's/BAL/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    BAL ALA  BETA-ALANINE
#sed -i 's/BNN/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    BNN ALA  ACETYL-P-AMIDINOPHENYLALANINE
#sed -i 's/C12/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    C12 ALA
#sed -i 's/C99/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    C99 ALA
#sed -i 's/CAB/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    CAB ALA  4-CARBOXY-4-AMINOBUTANAL
#sed -i 's/CH6/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    CH6 ALA
#sed -i 's/CH7/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    CH7 ALA
#sed -i 's/CLB/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    CLB ALA
#sed -i 's/CLD/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    CLD ALA
#sed -i 's/CLV/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    CLV ALA
#sed -i 's/CQR/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    CQR ALA
#sed -i 's/CR2/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    CR2 ALA  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/CR5/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    CR5 ALA
#sed -i 's/CR7/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    CR7 ALA
#sed -i 's/CR8/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    CR8 ALA
#sed -i 's/CRK/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    CRK ALA
#sed -i 's/CRW/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    CRW ALA
#sed -i 's/CRX/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    CRX ALA
#sed -i 's/CSI/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    CSI ALA
#sed -i 's/CSY/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    CSY ALA  MODIFIED TYROSINE COMPLEX
#sed -i 's/CWR/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    CWR ALA
#sed -i 's/DAB/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    DAB ALA  24-DIAMINOBUTYRIC ACID
#sed -i 's/DAL/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    DAL ALA  D-ALANINE
#sed -i 's/DAM/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    DAM ALA  N-METHYL-ALPHA-BETA-DEHYDROALANINE
#sed -i 's/DBU/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    DBU ALA  (2E)-2-AMINOBUT-2-ENOIC ACID
#sed -i 's/DBZ/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    DBZ ALA  3-(BENZOYLAMINO)-L-ALANINE
#sed -i 's/DHA/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    DHA ALA  2-AMINO-ACRYLIC ACID
#sed -i 's/DPP/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    DPP ALA  DIAMMINOPROPANOIC ACID
#sed -i 's/FGL/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    FGL ALA  2-AMINOPROPANEDIOIC ACID
#sed -i 's/HHK/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    HHK ALA  (2S)-28-DIAMINOOCTANOIC ACID
#sed -i 's/HMF/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    HMF ALA  2-AMINO-4-PHENYL-BUTYRIC ACID
#sed -i 's/IAM/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    IAM ALA  4-[(ISOPROPYLAMINO)METHYL]PHENYLALANINE
#sed -i 's/IGL/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    IGL ALA  ALPHA-AMINO-2-INDANACETIC ACID
#sed -i 's/KYN/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    KYN ALA  KYNURENINE
#sed -i 's/LAL/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    LAL ALA  NN-DIMETHYL-L-ALANINE
#sed -i 's/MAA/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    MAA ALA  N-METHYLALANINE
#sed -i 's/MDO/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    MDO ALA
#sed -i 's/MFC/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    MFC ALA  CYCLIZED
#sed -i 's/NAL/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    NAL ALA  BETA-(2-NAPHTHYL)-ALANINE
#sed -i 's/NAM/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    NAM ALA  NAM NAPTHYLAMINOALANINE
#sed -i 's/NCB/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    NCB ALA  CHEMICAL MODIFICATION
#sed -i 's/NRQ/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    NRQ ALA
#sed -i 's/NYC/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    NYC ALA
#sed -i 's/ORN/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    ORN ALA  ORNITHINE
#sed -i 's/PIA/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    PIA ALA  FUSION OF ALA 65 TYR 66 GLY 67
#sed -i 's/PRR/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    PRR ALA  3-(METHYL-PYRIDINIUM)ALANINE
#sed -i 's/PYA/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    PYA ALA  3-(110-PHENANTHROL-2-YL)-L-ALANINE
#sed -i 's/PYC/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    PYC ALA  PYRROLE-2-CARBOXYLATE
#sed -i 's/PYT/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    PYT ALA  MODIFIED ALANINE
#sed -i 's/RC7/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    RC7 ALA
#sed -i 's/SEC/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    SEC ALA  2-AMINO-3-SELENINO-PROPIONIC ACID
#sed -i 's/SIC/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    SIC ALA
#sed -i 's/SUI/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    SUI ALA
#sed -i 's/TIH/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    TIH ALA  BETA(2-THIENYL)ALANINE
#sed -i 's/TPQ/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    TPQ ALA  245-TRIHYDROXYPHENYLALANINE
#sed -i 's/UMA/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    UMA ALA
#sed -i 's/X9Q/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    X9Q ALA
#sed -i 's/XXY/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    XXY ALA
#sed -i 's/XYG/ALA/g' ${model%.*}_r_prep1.pdb                 ###                    XYG ALA
#sed -i 's/BCS/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/BUC/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/C5C/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/C6C/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CCS/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CEA/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CME/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSO/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSP/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSS/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSX/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CSW/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CY1/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CY3/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CYG/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
sed -i 's/CYM/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/CYS/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  CYS
#sed -i 's/CYQ/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/DCY/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/EFC/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/OCS/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/PEC/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/PR3/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SCH/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SCS/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SCY/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SHC/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SMC/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/SOC/CYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS CYS
#sed -i 's/5CS/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    5CS CYS
#sed -i 's/AGT/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    AGT CYS  AGMATINE-CYSTEINE ADDUCT
#sed -i 's/BBC/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    BBC CYS
#sed -i 's/BCS/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    BCS CYS  BENZYLCYSTEINE
#sed -i 's/BCX/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    BCX CYS  BETA-3-CYSTEINE
#sed -i 's/BPE/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    BPE CYS
#sed -i 's/BUC/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    BUC CYS  SS-BUTYLTHIOCYSTEINE
#sed -i 's/C3Y/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    C3Y CYS  MODIFIED CYSTEINE
#sed -i 's/C5C/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    C5C CYS  S-CYCLOPENTYL THIOCYSTEINE
#sed -i 's/C6C/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    C6C CYS  S-CYCLOHEXYL THIOCYSTEINE
#sed -i 's/CAF/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CAF CYS  S-DIMETHYLARSINOYL-CYSTEINE
#sed -i 's/CAS/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CAS CYS  S-(DIMETHYLARSENIC)CYSTEINE
#sed -i 's/CCS/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CCS CYS  CARBOXYMETHYLATED CYSTEINE
#sed -i 's/CME/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CME CYS  MODIFIED CYSTEINE
#sed -i 's/CML/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CML CYS
#sed -i 's/CMT/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CMT CYS  O-METHYLCYSTEINE
#sed -i 's/CS1/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CS1 CYS  S-(2-ANILINYL-SULFANYL)-CYSTEINE
#sed -i 's/CS3/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CS3 CYS
#sed -i 's/CS4/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CS4 CYS
#sed -i 's/CSA/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CSA CYS  S-ACETONYLCYSTEIN
#sed -i 's/CSB/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CSB CYS  CYS BOUND TO LEAD ION
#sed -i 's/CSD/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CSD CYS  3-SULFINOALANINE
#sed -i 's/CSE/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CSE CYS  SELENOCYSTEINE
#sed -i 's/CSO/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CSO CYS  INE S-HYDROXYCYSTEINE
#sed -i 's/CSR/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CSR CYS  S-ARSONOCYSTEINE
#sed -i 's/CSS/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CSS CYS  13-THIAZOLE-4-CARBOXYLIC ACID
#sed -i 's/CSU/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CSU CYS  CYSTEINE-S-SULFONIC ACID
#sed -i 's/CSW/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CSW CYS  CYSTEINE-S-DIOXIDE
#sed -i 's/CSX/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CSX CYS  OXOCYSTEINE
#sed -i 's/CSZ/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CSZ CYS  S-SELANYL CYSTEINE
#sed -i 's/CY0/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CY0 CYS  MODIFIED CYSTEINE
#sed -i 's/CY1/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CY1 CYS  ACETAMIDOMETHYLCYSTEINE
#sed -i 's/CY3/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CY3 CYS  2-AMINO-3-MERCAPTO-PROPIONAMIDE
#sed -i 's/CY4/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CY4 CYS  S-BUTYRYL-CYSTEIN
#sed -i 's/CY7/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CY7 CYS  MODIFIED CYSTEINE
#sed -i 's/CYD/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CYD CYS
#sed -i 's/CYF/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CYF CYS  FLUORESCEIN LABELLED CYS380 (P14)
#sed -i 's/CYG/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CYG CYS
#sed -i 's/CYQ/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CYQ CYS
#sed -i 's/CYR/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CYR CYS
#sed -i 's/CYS/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CYS CYS
#sed -i 's/CZ2/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CZ2 CYS  S-(DIHYDROXYARSINO)CYSTEINE
#sed -i 's/CZZ/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    CZZ CYS  THIARSAHYDROXY-CYSTEINE
#sed -i 's/DCY/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    DCY CYS  D-CYSTEINE
#sed -i 's/DYS/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    DYS CYS
#sed -i 's/EFC/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    EFC CYS  SS-(2-FLUOROETHYL)THIOCYSTEINE
#sed -i 's/FOE/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    FOE CYS
#sed -i 's/GT9/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    GT9 CYS  SG ALKYLATED
#sed -i 's/GYC/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    GYC CYS
#sed -i 's/HTI/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    HTI CYS
#sed -i 's/KOR/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    KOR CYS  MODIFIED CYSTEINE
#sed -i 's/M0H/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    M0H CYS  S-(HYDROXYMETHYL)-L-CYSTEINE
#sed -i 's/MCS/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    MCS CYS  MALONYLCYSTEINE
#sed -i 's/NPH/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    NPH CYS
#sed -i 's/NYS/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    NYS CYS
#sed -i 's/OCS/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    OCS CYS  CYSTEINE SULFONIC ACID
#sed -i 's/OCY/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    OCY CYS  HYDROXYETHYLCYSTEINE
#sed -i 's/P1L/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    P1L CYS  S-PALMITOYL CYSTEINE
#sed -i 's/PBB/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    PBB CYS  S-(4-BROMOBENZYL)CYSTEINE
#sed -i 's/PEC/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    PEC CYS  SS-PENTYLTHIOCYSTEINE
#sed -i 's/PR3/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    PR3 CYS  INE DTT-CYSTEINE
#sed -i 's/PYX/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    PYX CYS  S-[S-THIOPYRIDOXAMINYL]CYSTEINE
#sed -i 's/R1A/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    R1A CYS
#sed -i 's/R1B/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    R1B CYS
#sed -i 's/R1F/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    R1F CYS
#sed -i 's/R7A/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    R7A CYS
#sed -i 's/RCY/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    RCY CYS
#sed -i 's/SAH/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    SAH CYS  S-ADENOSYL-L-HOMOCYSTEINE
#sed -i 's/SC2/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    SC2 CYS  N-ACETYL-L-CYSTEINE
#sed -i 's/SCH/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    SCH CYS  S-METHYL THIOCYSTEINE GROUP
#sed -i 's/SCS/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    SCS CYS  MODIFIED CYSTEINE
#sed -i 's/SCY/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    SCY CYS  CETYLATED CYSTEINE
#sed -i 's/SHC/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    SHC CYS  S-HEXYLCYSTEINE
#sed -i 's/SMC/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    SMC CYS  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/SNC/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    SNC CYS  S-NITROSO CYSTEINE
#sed -i 's/SOC/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    SOC CYS  DIOXYSELENOCYSTEINE
#sed -i 's/TEE/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    TEE CYS  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/TNB/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    TNB CYS  S-(236-TRINITROPHENYL)CYSTEINE
#sed -i 's/TYX/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    TYX CYS  S-(2-ANILINO-2-OXOETHYL)-L-CYSTEINE
#sed -i 's/YCM/CYS/g' ${model%.*}_r_prep1.pdb                 ###                    YCM CYS  S-(2-AMINO-2-OXOETHYL)-L-CYSTEINE
#sed -i 's/2AS/ASP/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASA/ASP/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASB/ASP/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASK/ASP/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASL/ASP/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/ASP/ASP/g' ${model%.*}_r_prep1.pdb                 ###                  ASP
#sed -i 's/ASQ/ASP/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/BHD/ASP/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/DAS/ASP/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/DSP/ASP/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASP
#sed -i 's/3MD/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    3MD ASP  2S3S-3-METHYLASPARTIC ACID
#sed -i 's/A0A/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    A0A ASP  ASPARTYL-FORMYL MIXED ANHYDRIDE
#sed -i 's/ACB/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    ACB ASP  3-METHYL-ASPARTIC ACID
#sed -i 's/AKL/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    AKL ASP  3-AMINO-5-CHLORO-4-OXOPENTANOIC ACID
#sed -i 's/ASA/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    ASA ASP  ASPARTIC ALDEHYDE
#sed -i 's/ASB/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    ASB ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
#sed -i 's/ASI/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    ASI ASP  L-ISO-ASPARTATE
#sed -i 's/ASK/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    ASK ASP  DEHYDROXYMETHYLASPARTIC ACID
#sed -i 's/ASL/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    ASL ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
#sed -i 's/ASP/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    ASP ASP
#sed -i 's/B3D/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    B3D ASP  3-AMINOPENTANEDIOIC ACID
#sed -i 's/BFD/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    BFD ASP  ASPARTATE BERYLLIUM FLUORIDE
#sed -i 's/BHD/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    BHD ASP  BETA-HYDROXYASPARTIC ACID
#sed -i 's/DAS/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    DAS ASP  D-ASPARTIC ACID
#sed -i 's/DMK/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    DMK ASP  DIMETHYL ASPARTIC ACID
#sed -i 's/IAS/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    IAS ASP  ASPARTYL GROUP
#sed -i 's/OHS/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    OHS ASP  O-(CARBOXYSULFANYL)-4-OXO-L-HOMOSERINE
#sed -i 's/OXX/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    OXX ASP  OXALYL-ASPARTYL ANHYDRIDE
#sed -i 's/PHD/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    PHD ASP  2-AMINO-4-OXO-4-PHOSPHONOOXY-BUTYRIC ACID
#sed -i 's/SNN/ASP/g' ${model%.*}_r_prep1.pdb                 ###                    SNN ASP  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/5HP/GLU/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/CGU/GLU/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/DGL/GLU/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/GGL/GLU/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/GLU/GLU/g' ${model%.*}_r_prep1.pdb                 ###                  GLU
#sed -i 's/GMA/GLU/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/PCA/GLU/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLU
#sed -i 's/AB7/GLU/g' ${model%.*}_r_prep1.pdb                 ###                    AB7 GLU  ALPHA-AMINOBUTYRIC ACID
#sed -i 's/AR4/GLU/g' ${model%.*}_r_prep1.pdb                 ###                    AR4 GLU
#sed -i 's/B3E/GLU/g' ${model%.*}_r_prep1.pdb                 ###                    B3E GLU  (3S)-3-AMINOHEXANEDIOIC ACID
#sed -i 's/CGU/GLU/g' ${model%.*}_r_prep1.pdb                 ###                    CGU GLU  CARBOXYLATION OF THE CG ATOM
#sed -i 's/DGL/GLU/g' ${model%.*}_r_prep1.pdb                 ###                    DGL GLU  D-GLU
#sed -i 's/GLU/GLU/g' ${model%.*}_r_prep1.pdb                 ###                    GLU GLU
#sed -i 's/GMA/GLU/g' ${model%.*}_r_prep1.pdb                 ###                    GMA GLU  1-AMIDO-GLUTAMIC ACID
#sed -i 's/ILG/GLU/g' ${model%.*}_r_prep1.pdb                 ###                    ILG GLU  GLU LINKED TO NEXT RESIDUE VIA CG
#sed -i 's/LME/GLU/g' ${model%.*}_r_prep1.pdb                 ###                    LME GLU  (3R)-3-METHYL-L-GLUTAMIC ACID
#sed -i 's/MEG/GLU/g' ${model%.*}_r_prep1.pdb                 ###                    MEG GLU  (2S3R)-3-METHYL-GLUTAMIC ACID
#sed -i 's/DAH/PHE/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/DPN/PHE/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/HPQ/PHE/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/PHE/PHE/g' ${model%.*}_r_prep1.pdb                 ###                  PHE
#sed -i 's/PHI/PHE/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/PHL/PHE/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PHE
#sed -i 's/1PA/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    1PA PHE  PHENYLMETHYLACETIC ACID ALANINE
#sed -i 's/23F/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    23F PHE  (2Z)-2-AMINO-3-PHENYLACRYLIC ACID
#sed -i 's/4PH/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    4PH PHE  4-METHYL-L-PHENYLALANINE
#sed -i 's/B2F/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    B2F PHE  PHENYLALANINE BORONIC ACID
#sed -i 's/BIF/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    BIF PHE
#sed -i 's/CHS/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    CHS PHE  4-AMINO-5-CYCLOHEXYL-3-HYDROXY-PENTANOIC AC
#sed -i 's/DAH/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    DAH PHE  34-DIHYDROXYDAHNYLALANINE
#sed -i 's/DPH/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    DPH PHE  DEAMINO-METHYL-PHENYLALANINE
#sed -i 's/DPN/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    DPN PHE  D-CONFIGURATION
#sed -i 's/FCL/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    FCL PHE  3-CHLORO-L-PHENYLALANINE
#sed -i 's/FOG/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    FOG PHE  PHENYLALANINOYL-[1-HYDROXY]-2-PROPYLENE
#sed -i 's/FRF/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    FRF PHE  PHE FOLLOWED BY REDUCED PHE
#sed -i 's/HPE/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    HPE PHE  HOMOPHENYLALANINE
#sed -i 's/HPH/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    HPH PHE  PHENYLALANINOL GROUP
#sed -i 's/HPQ/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    HPQ PHE  HOMOPHENYLALANINYLMETHANE
#sed -i 's/MEA/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    MEA PHE  N-METHYLPHENYLALANINE
#sed -i 's/MTY/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    MTY PHE  3-HYDROXYPHENYLALANINE
#sed -i 's/NFA/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    NFA PHE  MODIFIED PHENYLALANINE
#sed -i 's/PBF/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    PBF PHE  PARA-(BENZOYL)-PHENYLALANINE
#sed -i 's/PCS/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    PCS PHE  PHENYLALANYLMETHYLCHLORIDE
#sed -i 's/PF5/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    PF5 PHE  23456-PENTAFLUORO-L-PHENYLALANINE
#sed -i 's/PFF/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    PFF PHE  4-FLUORO-L-PHENYLALANINE
#sed -i 's/PHA/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    PHA PHE  PHENYLALANINAL
#sed -i 's/PHE/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    PHE PHE
#sed -i 's/PHI/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    PHI PHE  IODO-PHENYLALANINE
#sed -i 's/PHL/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    PHL PHE  L-PHENYLALANINOL
#sed -i 's/PHM/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    PHM PHE  PHENYLALANYLMETHANE
#sed -i 's/PM3/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    PM3 PHE
#sed -i 's/PPN/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    PPN PHE  THE LIGAND IS A PARA-NITRO-PHENYLALANINE
#sed -i 's/PRQ/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    PRQ PHE  PHENYLALANINE
#sed -i 's/PSA/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    PSA PHE
#sed -i 's/SMF/PHE/g' ${model%.*}_r_prep1.pdb                 ###                    SMF PHE  4-SULFOMETHYL-L-PHENYLALANINE
#sed -i 's/GL3/GLY/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/GLY/GLY/g' ${model%.*}_r_prep1.pdb                 ###                  GLY
#sed -i 's/GLZ/GLY/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/GSC/GLY/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/MPQ/GLY/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/MSA/GLY/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/NMC/GLY/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/SAR/GLY/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLY
#sed -i 's/ACY/GLY/g' ${model%.*}_r_prep1.pdb                 ###                    ACY GLY  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/CHG/GLY/g' ${model%.*}_r_prep1.pdb                 ###                    CHG GLY  CYCLOHEXYL GLYCINE
#sed -i 's/CHP/GLY/g' ${model%.*}_r_prep1.pdb                 ###                    CHP GLY  3-CHLORO-4-HYDROXYPHENYLGLYCINE
#sed -i 's/GHP/GLY/g' ${model%.*}_r_prep1.pdb                 ###                    GHP GLY  4-HYDROXYPHENYLGLYCINE
#sed -i 's/GL3/GLY/g' ${model%.*}_r_prep1.pdb                 ###                    GL3 GLY  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/GLY/GLY/g' ${model%.*}_r_prep1.pdb                 ###                    GLY GLY
#sed -i 's/GLZ/GLY/g' ${model%.*}_r_prep1.pdb                 ###                    GLZ GLY  AMINO-ACETALDEHYDE
#sed -i 's/GYS/GLY/g' ${model%.*}_r_prep1.pdb                 ###                    GYS GLY
#sed -i 's/IPG/GLY/g' ${model%.*}_r_prep1.pdb                 ###                    IPG GLY  N-ISOPROPYL GLYCINE
#sed -i 's/MEU/GLY/g' ${model%.*}_r_prep1.pdb                 ###                    MEU GLY  O-METHYL-GLYCINE
#sed -i 's/MPQ/GLY/g' ${model%.*}_r_prep1.pdb                 ###                    MPQ GLY  N-METHYL-ALPHA-PHENYL-GLYCINE
#sed -i 's/MSA/GLY/g' ${model%.*}_r_prep1.pdb                 ###                    MSA GLY  (2-S-METHYL) SARCOSINE
#sed -i 's/NMC/GLY/g' ${model%.*}_r_prep1.pdb                 ###                    NMC GLY  N-CYCLOPROPYLMETHYL GLYCINE
#sed -i 's/PG9/GLY/g' ${model%.*}_r_prep1.pdb                 ###                    PG9 GLY  D-PHENYLGLYCINE
#sed -i 's/SAR/GLY/g' ${model%.*}_r_prep1.pdb                 ###                    SAR GLY  SARCOSINE
#sed -i 's/SHP/GLY/g' ${model%.*}_r_prep1.pdb                 ###                    SHP GLY  (4-HYDROXYMALTOSEPHENYL)GLYCINE
#sed -i 's/TBG/GLY/g' ${model%.*}_r_prep1.pdb                 ###                    TBG GLY  T-BUTYL GLYCINE
#sed -i 's/3AH/HIS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/DHI/HIS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/HIC/HIS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/HIS/HIS/g' ${model%.*}_r_prep1.pdb                 ###                  HIS
#sed -i 's/MHS/HIS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/NEM/HIS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/NEP/HIS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS HIS
#sed -i 's/HID/HIS/g' ${model%.*}_r_prep1.pdb                 ###                  single delta N protonation
#sed -i 's/HIE/HIS/g' ${model%.*}_r_prep1.pdb                 ###                  single epsilon N protonation
#sed -i 's/3AH/HIS/g' ${model%.*}_r_prep1.pdb                 ###                    3AH HIS
#sed -i 's/DDE/HIS/g' ${model%.*}_r_prep1.pdb                 ###                    DDE HIS
#sed -i 's/DHI/HIS/g' ${model%.*}_r_prep1.pdb                 ###                    DHI HIS  D-HISTIDINE
#sed -i 's/HIA/HIS/g' ${model%.*}_r_prep1.pdb                 ###                    HIA HIS  L-HISTIDINE AMIDE
#sed -i 's/HIC/HIS/g' ${model%.*}_r_prep1.pdb                 ###                    HIC HIS  4-METHYL-HISTIDINE
#sed -i 's/HIP/HIS/g' ${model%.*}_r_prep1.pdb                 ###                    HIP HIS  ND1-PHOSPHONOHISTIDINE...or commonly used doubly protonated state
#sed -i 's/HIQ/HIS/g' ${model%.*}_r_prep1.pdb                 ###                    HIQ HIS  MODIFIED HISTIDINE
#sed -i 's/HIS/HIS/g' ${model%.*}_r_prep1.pdb                 ###                    HIS HIS
#sed -i 's/HSO/HIS/g' ${model%.*}_r_prep1.pdb                 ###                    HSO HIS  HISTIDINOL
#sed -i 's/MHS/HIS/g' ${model%.*}_r_prep1.pdb                 ###                    MHS HIS  1-N-METHYLHISTIDINE
#sed -i 's/NEP/HIS/g' ${model%.*}_r_prep1.pdb                 ###                    NEP HIS  N1-PHOSPHONOHISTIDINE
#sed -i 's/NZH/HIS/g' ${model%.*}_r_prep1.pdb                 ###                    NZH HIS
#sed -i 's/OHI/HIS/g' ${model%.*}_r_prep1.pdb                 ###                    OHI HIS  3-(2-OXO-2H-IMIDAZOL-4-YL)-L-ALANINE
#sed -i 's/PSH/HIS/g' ${model%.*}_r_prep1.pdb                 ###                    PSH HIS  1-THIOPHOSPHONO-L-HISTIDINE
#sed -i 's/DIL/ILE/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ILE
#sed -i 's/IIL/ILE/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ILE
#sed -i 's/ILE/ILE/g' ${model%.*}_r_prep1.pdb                 ###                  ILE
#sed -i 's/B2I/ILE/g' ${model%.*}_r_prep1.pdb                 ###                    B2I ILE  ISOLEUCINE BORONIC ACID
#sed -i 's/DIL/ILE/g' ${model%.*}_r_prep1.pdb                 ###                    DIL ILE  D-ISOLEUCINE
#sed -i 's/IIL/ILE/g' ${model%.*}_r_prep1.pdb                 ###                    IIL ILE  ISO-ISOLEUCINE
#sed -i 's/ILE/ILE/g' ${model%.*}_r_prep1.pdb                 ###                    ILE ILE
#sed -i 's/ILX/ILE/g' ${model%.*}_r_prep1.pdb                 ###                    ILX ILE  45-DIHYDROXYISOLEUCINE
#sed -i 's/IML/ILE/g' ${model%.*}_r_prep1.pdb                 ###                    IML ILE  N-METHYLATED
#sed -i 's/ALY/LYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/DLY/LYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/KCX/LYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/LLP/LYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/LLY/LYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/LYM/LYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/LYS/LYS/g' ${model%.*}_r_prep1.pdb                 ###                  LYS
#sed -i 's/LYZ/LYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/MLY/LYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/SHR/LYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/TRG/LYS/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LYS
#sed -i 's/6CL/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    6CL LYS  6-CARBOXYLYSINE
#sed -i 's/ALY/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    ALY LYS  N(6)-ACETYLLYSINE
#sed -i 's/API/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    API LYS  26-DIAMINOPIMELIC ACID
#sed -i 's/APK/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    APK LYS
#sed -i 's/AZK/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    AZK LYS  (2S)-2-AMINO-6-TRIAZANYLHEXAN-1-OL
#sed -i 's/B3K/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    B3K LYS  (3S)-37-DIAMINOHEPTANOIC ACID
#sed -i 's/BLY/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    BLY LYS  LYSINE BORONIC ACID
#sed -i 's/C1X/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    C1X LYS  MODIFIED LYSINE
#sed -i 's/CLG/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    CLG LYS
#sed -i 's/CLH/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    CLH LYS
#sed -i 's/CYJ/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    CYJ LYS  MODIFIED LYSINE
#sed -i 's/DLS/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    DLS LYS  DI-ACETYL-LYSINE
#sed -i 's/DLY/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    DLY LYS  D-LYSINE
#sed -i 's/DNL/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    DNL LYS  6-AMINO-HEXANAL
#sed -i 's/FHL/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    FHL LYS  MODIFIED LYSINE
#sed -i 's/GPL/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    GPL LYS  LYSINE GUANOSINE-5-MONOPHOSPHATE
#sed -i 's/IT1/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    IT1 LYS
#sed -i 's/KCX/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    KCX LYS  CARBAMOYLATED LYSINE
#sed -i 's/KGC/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    KGC LYS
#sed -i 's/KST/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    KST LYS  N~6~-(5-CARBOXY-3-THIENYL)-L-LYSINE
#sed -i 's/LA2/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    LA2 LYS
#sed -i 's/LCK/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    LCK LYS
#sed -i 's/LCX/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    LCX LYS  CARBAMYLATED LYSINE
#sed -i 's/LDH/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    LDH LYS  N~6~-ETHYL-L-LYSINE
#sed -i 's/LET/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    LET LYS  ODIFIED LYSINE
#sed -i 's/LLP/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    LLP LYS
#sed -i 's/LLY/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    LLY LYS  NZ-(DICARBOXYMETHYL)LYSINE
#sed -i 's/LSO/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    LSO LYS  MODIFIED LYSINE
#sed -i 's/LYM/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    LYM LYS  DEOXY-METHYL-LYSINE
#sed -i 's/LYN/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    LYN LYS  26-DIAMINO-HEXANOIC ACID AMIDE
#sed -i 's/LYP/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    LYP LYS  N~6~-METHYL-N~6~-PROPYL-L-LYSINE
#sed -i 's/LYR/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    LYR LYS  MODIFIED LYSINE
#sed -i 's/LYS/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    LYS LYS
#sed -i 's/LYX/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    LYX LYS  N-(2-COENZYME A)-PROPANOYL-LYSINE
#sed -i 's/LYZ/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    LYZ LYS  5-HYDROXYLYSINE
#sed -i 's/M2L/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    M2L LYS
#sed -i 's/M3L/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    M3L LYS  N-TRIMETHYLLYSINE
#sed -i 's/MCL/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    MCL LYS  NZ-(1-CARBOXYETHYL)-LYSINE
#sed -i 's/MLY/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    MLY LYS  METHYLATED LYSINE
#sed -i 's/MLZ/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    MLZ LYS  N-METHYL-LYSINE
#sed -i 's/OBS/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    OBS LYS  MODIFIED LYSINE
#sed -i 's/SLZ/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    SLZ LYS  L-THIALYSINE
#sed -i 's/XX1/LYS/g' ${model%.*}_r_prep1.pdb                 ###                    XX1 LYS  N~6~-7H-PURIN-6-YL-L-LYSINE
#sed -i 's/BUG/LEU/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/CLE/LEU/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/DLE/LEU/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/LEU/LEU/g' ${model%.*}_r_prep1.pdb                 ###                  LEU
#sed -i 's/MLE/LEU/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/NLE/LEU/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/NLN/LEU/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/NLP/LEU/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS LEU
#sed -i 's/1LU/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    1LU LEU  4-METHYL-PENTANOIC ACID-2-OXYL GROUP
#sed -i 's/2ML/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    2ML LEU  2-METHYLLEUCINE
#sed -i 's/BLE/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    BLE LEU  LEUCINE BORONIC ACID
#sed -i 's/BUG/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    BUG LEU  TERT-LEUCYL AMINE
#sed -i 's/CLE/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    CLE LEU  LEUCINE AMIDE
#sed -i 's/DCL/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    DCL LEU  2-AMINO-4-METHYL-PENTANYL GROUP
#sed -i 's/DLE/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    DLE LEU  D-LEUCINE
#sed -i 's/DNE/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    DNE LEU  D-NORLEUCINE
#sed -i 's/DNG/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    DNG LEU  N-FORMYL-D-NORLEUCINE
#sed -i 's/DNM/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    DNM LEU  D-N-METHYL NORLEUCINE
#sed -i 's/FLE/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    FLE LEU  FUROYL-LEUCINE
#sed -i 's/HLU/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    HLU LEU  BETA-HYDROXYLEUCINE
#sed -i 's/LED/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    LED LEU  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/LEF/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    LEF LEU  2-5-FLUOROLEUCINE
#sed -i 's/LEU/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    LEU LEU
#sed -i 's/LNT/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    LNT LEU
#sed -i 's/MHL/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    MHL LEU  N-METHYLATED HYDROXY
#sed -i 's/MLE/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    MLE LEU  N-METHYLATED
#sed -i 's/MLL/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    MLL LEU  METHYL L-LEUCINATE
#sed -i 's/MNL/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    MNL LEU  4N-DIMETHYLNORLEUCINE
#sed -i 's/NLE/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    NLE LEU  NORLEUCINE
#sed -i 's/NLN/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    NLN LEU  NORLEUCINE AMIDE
#sed -i 's/NLO/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    NLO LEU  O-METHYL-L-NORLEUCINE
#sed -i 's/PLE/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    PLE LEU  LEUCINE PHOSPHINIC ACID
#sed -i 's/PPH/LEU/g' ${model%.*}_r_prep1.pdb                 ###                    PPH LEU  PHENYLALANINE PHOSPHINIC ACID
#sed -i 's/CXM/MET/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
#sed -i 's/FME/MET/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
#sed -i 's/MET/MET/g' ${model%.*}_r_prep1.pdb                 ###                  MET
#sed -i 's/MSE/MET/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
#sed -i 's/OMT/MET/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS MET
#sed -i 's/AME/MET/g' ${model%.*}_r_prep1.pdb                 ###                    AME MET  ACETYLATED METHIONINE
#sed -i 's/CXM/MET/g' ${model%.*}_r_prep1.pdb                 ###                    CXM MET  N-CARBOXYMETHIONINE
#sed -i 's/ESC/MET/g' ${model%.*}_r_prep1.pdb                 ###                    ESC MET  2-AMINO-4-ETHYL SULFANYL BUTYRIC ACID
#sed -i 's/FME/MET/g' ${model%.*}_r_prep1.pdb                 ###                    FME MET  FORMYL-METHIONINE
#sed -i 's/FOR/MET/g' ${model%.*}_r_prep1.pdb                 ###                    FOR MET
#sed -i 's/MET/MET/g' ${model%.*}_r_prep1.pdb                 ###                    MET MET
#sed -i 's/MHO/MET/g' ${model%.*}_r_prep1.pdb                 ###                    MHO MET  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/MME/MET/g' ${model%.*}_r_prep1.pdb                 ###                    MME MET  N-METHYL METHIONINE
#sed -i 's/MSE/MET/g' ${model%.*}_r_prep1.pdb                 ###                    MSE MET  ELENOMETHIONINE
#sed -i 's/MSO/MET/g' ${model%.*}_r_prep1.pdb                 ###                    MSO MET  METHIONINE SULFOXIDE
#sed -i 's/OMT/MET/g' ${model%.*}_r_prep1.pdb                 ###                    OMT MET  METHIONINE SULFONE
#sed -i 's/SME/MET/g' ${model%.*}_r_prep1.pdb                 ###                    SME MET  METHIONINE SULFOXIDE
#sed -i 's/ASN/ASN/g' ${model%.*}_r_prep1.pdb                 ###                  ASN
#sed -i 's/MEN/ASN/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ASN
#sed -i 's/AFA/ASN/g' ${model%.*}_r_prep1.pdb                 ###                    AFA ASN  N-[7-METHYL-OCT-24-DIENOYL]ASPARAGINE
#sed -i 's/AHB/ASN/g' ${model%.*}_r_prep1.pdb                 ###                    AHB ASN  BETA-HYDROXYASPARAGINE
#sed -i 's/ASN/ASN/g' ${model%.*}_r_prep1.pdb                 ###                    ASN ASN
#sed -i 's/B3X/ASN/g' ${model%.*}_r_prep1.pdb                 ###                    B3X ASN  (3S)-35-DIAMINO-5-OXOPENTANOIC ACID
#sed -i 's/DMH/ASN/g' ${model%.*}_r_prep1.pdb                 ###                    DMH ASN  N4N4-DIMETHYL-ASPARAGINE
#sed -i 's/DSG/ASN/g' ${model%.*}_r_prep1.pdb                 ###                    DSG ASN  D-ASPARAGINE
#sed -i 's/MEN/ASN/g' ${model%.*}_r_prep1.pdb                 ###                    MEN ASN  GAMMA METHYL ASPARAGINE
#sed -i 's/DPR/PRO/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS PRO
#sed -i 's/PRO/PRO/g' ${model%.*}_r_prep1.pdb                 ###                  PRO
#sed -i 's/1AB/PRO/g' ${model%.*}_r_prep1.pdb                 ###                    1AB PRO  14-DIDEOXY-14-IMINO-D-ARABINITOL
#sed -i 's/2MT/PRO/g' ${model%.*}_r_prep1.pdb                 ###                    2MT PRO
#sed -i 's/4FB/PRO/g' ${model%.*}_r_prep1.pdb                 ###                    4FB PRO  (4S)-4-FLUORO-L-PROLINE
#sed -i 's/DPL/PRO/g' ${model%.*}_r_prep1.pdb                 ###                    DPL PRO  4-OXOPROLINE
#sed -i 's/DPR/PRO/g' ${model%.*}_r_prep1.pdb                 ###                    DPR PRO  D-PROLINE
#sed -i 's/H5M/PRO/g' ${model%.*}_r_prep1.pdb                 ###                    H5M PRO  TRANS-3-HYDROXY-5-METHYLPROLINE
#sed -i 's/HY3/PRO/g' ${model%.*}_r_prep1.pdb                 ###                    HY3 PRO  3-HYDROXYPROLINE
#sed -i 's/HYP/PRO/g' ${model%.*}_r_prep1.pdb                 ###                    HYP PRO  4-HYDROXYPROLINE
#sed -i 's/LPD/PRO/g' ${model%.*}_r_prep1.pdb                 ###                    LPD PRO  L-PROLINAMIDE
#sed -i 's/P2Y/PRO/g' ${model%.*}_r_prep1.pdb                 ###                    P2Y PRO  (2S)-PYRROLIDIN-2-YLMETHYLAMINE
#sed -i 's/PCA/PRO/g' ${model%.*}_r_prep1.pdb                 ###                    PCA PRO  5-OXOPROLINE
#sed -i 's/POM/PRO/g' ${model%.*}_r_prep1.pdb                 ###                    POM PRO  CIS-5-METHYL-4-OXOPROLINE
#sed -i 's/PRO/PRO/g' ${model%.*}_r_prep1.pdb                 ###                    PRO PRO
#sed -i 's/PRS/PRO/g' ${model%.*}_r_prep1.pdb                 ###                    PRS PRO  THIOPROLINE
#sed -i 's/DGN/GLN/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS GLN
#sed -i 's/GLN/GLN/g' ${model%.*}_r_prep1.pdb                 ###                  GLN
#sed -i 's/DGN/GLN/g' ${model%.*}_r_prep1.pdb                 ###                    DGN GLN  D-GLUTAMINE
#sed -i 's/GHG/GLN/g' ${model%.*}_r_prep1.pdb                 ###                    GHG GLN  GAMMA-HYDROXY-GLUTAMINE
#sed -i 's/GLH/GLN/g' ${model%.*}_r_prep1.pdb                 ###                    GLH GLN
#sed -i 's/GLN/GLN/g' ${model%.*}_r_prep1.pdb                 ###                    GLN GLN
#sed -i 's/MGN/GLN/g' ${model%.*}_r_prep1.pdb                 ###                    MGN GLN  2-METHYL-GLUTAMINE
#sed -i 's/ACL/ARG/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/AGM/ARG/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/ARG/ARG/g' ${model%.*}_r_prep1.pdb                 ###                  ARG
#sed -i 's/ARM/ARG/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/DAR/ARG/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/HAR/ARG/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/HMR/ARG/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS ARG
#sed -i 's/2MR/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    2MR ARG  N3 N4-DIMETHYLARGININE
#sed -i 's/AAR/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    AAR ARG  ARGININEAMIDE
#sed -i 's/ACL/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    ACL ARG  DEOXY-CHLOROMETHYL-ARGININE
#sed -i 's/AGM/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    AGM ARG  4-METHYL-ARGININE
#sed -i 's/ALG/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    ALG ARG  GUANIDINOBUTYRYL GROUP
#sed -i 's/AR2/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    AR2 ARG  ARGINYL-BENZOTHIAZOLE-6-CARBOXYLIC ACID
#sed -i 's/ARG/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    ARG ARG
#sed -i 's/ARM/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    ARM ARG  DEOXY-METHYL-ARGININE
#sed -i 's/ARO/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    ARO ARG  C-GAMMA-HYDROXY ARGININE
#sed -i 's/BOR/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    BOR ARG
#sed -i 's/CIR/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    CIR ARG  CITRULLINE
#sed -i 's/DA2/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    DA2 ARG  MODIFIED ARGININE
#sed -i 's/DAR/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    DAR ARG  D-ARGININE
#sed -i 's/HMR/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    HMR ARG  BETA-HOMOARGININE
#sed -i 's/HRG/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    HRG ARG  L-HOMOARGININE
#sed -i 's/MAI/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    MAI ARG  DEOXO-METHYLARGININE
#sed -i 's/MGG/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    MGG ARG  MODIFIED D-ARGININE
#sed -i 's/NMM/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    NMM ARG  MODIFIED ARGININE
#sed -i 's/OPR/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    OPR ARG  C-(3-OXOPROPYL)ARGININE
#sed -i 's/ORQ/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    ORQ ARG  N~5~-ACETYL-L-ORNITHINE
#sed -i 's/TYZ/ARG/g' ${model%.*}_r_prep1.pdb                 ###                    TYZ ARG  PARA ACETAMIDO BENZOIC ACID
#sed -i 's/DSN/SER/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/MIS/SER/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/OAS/SER/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SAC/SER/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SEL/SER/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SEP/SER/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SER/SER/g' ${model%.*}_r_prep1.pdb                 ###                  SER
#sed -i 's/SET/SER/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/SVA/SER/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS SER
#sed -i 's/B3S/SER/g' ${model%.*}_r_prep1.pdb                 ###                    B3S SER  (3R)-3-AMINO-4-HYDROXYBUTANOIC ACID
#sed -i 's/BG1/SER/g' ${model%.*}_r_prep1.pdb                 ###                    BG1 SER
#sed -i 's/DHL/SER/g' ${model%.*}_r_prep1.pdb                 ###                    DHL SER  POST-TRANSLATIONAL MODIFICATION
#sed -i 's/DSE/SER/g' ${model%.*}_r_prep1.pdb                 ###                    DSE SER  D-SERINE N-METHYLATED
#sed -i 's/DSN/SER/g' ${model%.*}_r_prep1.pdb                 ###                    DSN SER  D-SERINE
#sed -i 's/FGP/SER/g' ${model%.*}_r_prep1.pdb                 ###                    FGP SER
#sed -i 's/GVL/SER/g' ${model%.*}_r_prep1.pdb                 ###                    GVL SER  SERINE MODIFED WITH PHOSPHOPANTETHEINE
#sed -i 's/HSE/SER/g' ${model%.*}_r_prep1.pdb                 ###                    HSE SER  L-HOMOSERINE
#sed -i 's/HSL/SER/g' ${model%.*}_r_prep1.pdb                 ###                    HSL SER  HOMOSERINE LACTONE
#sed -i 's/MC1/SER/g' ${model%.*}_r_prep1.pdb                 ###                    MC1 SER  METHICILLIN ACYL-SERINE
#sed -i 's/MIS/SER/g' ${model%.*}_r_prep1.pdb                 ###                    MIS SER  MODIFIED SERINE
#sed -i 's/N10/SER/g' ${model%.*}_r_prep1.pdb                 ###                    N10 SER  O-[(HEXYLAMINO)CARBONYL]-L-SERINE
#sed -i 's/NC1/SER/g' ${model%.*}_r_prep1.pdb                 ###                    NC1 SER  NITROCEFIN ACYL-SERINE
#sed -i 's/OAS/SER/g' ${model%.*}_r_prep1.pdb                 ###                    OAS SER  O-ACETYLSERINE
#sed -i 's/OSE/SER/g' ${model%.*}_r_prep1.pdb                 ###                    OSE SER  O-SULFO-L-SERINE
#sed -i 's/PG1/SER/g' ${model%.*}_r_prep1.pdb                 ###                    PG1 SER  BENZYLPENICILLOYL-ACYLATED SERINE
#sed -i 's/PYR/SER/g' ${model%.*}_r_prep1.pdb                 ###                    PYR SER  CHEMICALLY MODIFIED
#sed -i 's/S1H/SER/g' ${model%.*}_r_prep1.pdb                 ###                    S1H SER  1-HEXADECANOSULFONYL-O-L-SERINE
#sed -i 's/SAC/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SAC SER  N-ACETYL-SERINE
#sed -i 's/SBD/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SBD SER
#sed -i 's/SBG/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SBG SER  MODIFIED SERINE
#sed -i 's/SBL/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SBL SER
#sed -i 's/SDP/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SDP SER
#sed -i 's/SEB/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SEB SER  O-BENZYLSULFONYL-SERINE
#sed -i 's/SEL/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SEL SER  2-AMINO-13-PROPANEDIOL
#sed -i 's/SEP/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SEP SER  E PHOSPHOSERINE
#sed -i 's/SER/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SER SER
#sed -i 's/SET/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SET SER  AMINOSERINE
#sed -i 's/SGB/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SGB SER  MODIFIED SERINE
#sed -i 's/SGR/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SGR SER  MODIFIED SERINE
#sed -i 's/SOY/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SOY SER  OXACILLOYL-ACYLATED SERINE
#sed -i 's/SUN/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SUN SER  TABUN CONJUGATED SERINE
#sed -i 's/SVA/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SVA SER  SERINE VANADATE
#sed -i 's/SVV/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SVV SER  MODIFIED SERINE
#sed -i 's/SVX/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SVX SER  MODIFIED SERINE
#sed -i 's/SVY/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SVY SER  MODIFIED SERINE
#sed -i 's/SVZ/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SVZ SER  MODIFIED SERINE
#sed -i 's/SXE/SER/g' ${model%.*}_r_prep1.pdb                 ###                    SXE SER  MODIFIED SERINE
#sed -i 's/ALO/THR/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
#sed -i 's/BMT/THR/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
#sed -i 's/DTH/THR/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
#sed -i 's/THR/THR/g' ${model%.*}_r_prep1.pdb                 ###                  THR
#sed -i 's/TPO/THR/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS THR
#sed -i 's/AEI/THR/g' ${model%.*}_r_prep1.pdb                 ###                    AEI THR  ACYLATED THR
#sed -i 's/ALO/THR/g' ${model%.*}_r_prep1.pdb                 ###                    ALO THR  ALLO-THREONINE
#sed -i 's/BMT/THR/g' ${model%.*}_r_prep1.pdb                 ###                    BMT THR
#sed -i 's/CRO/THR/g' ${model%.*}_r_prep1.pdb                 ###                    CRO THR  CYCLIZED
#sed -i 's/CTH/THR/g' ${model%.*}_r_prep1.pdb                 ###                    CTH THR  4-CHLOROTHREONINE
#sed -i 's/DTH/THR/g' ${model%.*}_r_prep1.pdb                 ###                    DTH THR  D-THREONINE
#sed -i 's/OLT/THR/g' ${model%.*}_r_prep1.pdb                 ###                    OLT THR  O-METHYL-L-THREONINE
#sed -i 's/TBM/THR/g' ${model%.*}_r_prep1.pdb                 ###                    TBM THR
#sed -i 's/TH5/THR/g' ${model%.*}_r_prep1.pdb                 ###                    TH5 THR  O-ACETYL-L-THREONINE
#sed -i 's/THC/THR/g' ${model%.*}_r_prep1.pdb                 ###                    THC THR  N-METHYLCARBONYLTHREONINE
#sed -i 's/THR/THR/g' ${model%.*}_r_prep1.pdb                 ###                    THR THR
#sed -i 's/TMD/THR/g' ${model%.*}_r_prep1.pdb                 ###                    TMD THR  N-METHYLATED EPSILON C ALKYLATED
#sed -i 's/TPO/THR/g' ${model%.*}_r_prep1.pdb                 ###                    TPO THR  HOSPHOTHREONINE
#sed -i 's/DIV/VAL/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS VAL
#sed -i 's/DVA/VAL/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS VAL
#sed -i 's/MVA/VAL/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS VAL
#sed -i 's/VAL/VAL/g' ${model%.*}_r_prep1.pdb                 ###                  VAL
#sed -i 's/B2V/VAL/g' ${model%.*}_r_prep1.pdb                 ###                    B2V VAL  VALINE BORONIC ACID
#sed -i 's/DIV/VAL/g' ${model%.*}_r_prep1.pdb                 ###                    DIV VAL  D-ISOVALINE
#sed -i 's/DVA/VAL/g' ${model%.*}_r_prep1.pdb                 ###                    DVA VAL  D-VALINE
#sed -i 's/MNV/VAL/g' ${model%.*}_r_prep1.pdb                 ###                    MNV VAL  N-METHYL-C-AMINO VALINE
#sed -i 's/MVA/VAL/g' ${model%.*}_r_prep1.pdb                 ###                    MVA VAL  N-METHYLATED
#sed -i 's/NVA/VAL/g' ${model%.*}_r_prep1.pdb                 ###                    NVA VAL  NORVALINE
#sed -i 's/VAD/VAL/g' ${model%.*}_r_prep1.pdb                 ###                    VAD VAL  DEAMINOHYDROXYVALINE
#sed -i 's/VAF/VAL/g' ${model%.*}_r_prep1.pdb                 ###                    VAF VAL  METHYLVALINE
#sed -i 's/VAL/VAL/g' ${model%.*}_r_prep1.pdb                 ###                    VAL VAL
#sed -i 's/VDL/VAL/g' ${model%.*}_r_prep1.pdb                 ###                    VDL VAL  (2R3R)-23-DIAMINOBUTANOIC ACID
#sed -i 's/VLL/VAL/g' ${model%.*}_r_prep1.pdb                 ###                    VLL VAL  (2S)-23-DIAMINOBUTANOIC ACID
#sed -i 's/VME/VAL/g' ${model%.*}_r_prep1.pdb                 ###                    VME VAL  O- METHYLVALINE
#sed -i 's/DTR/TRP/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/HTR/TRP/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/LTR/TRP/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/TPL/TRP/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/TRO/TRP/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TRP
#sed -i 's/TRP/TRP/g' ${model%.*}_r_prep1.pdb                 ###                  TRP
#sed -i 's/BTR/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    BTR TRP  6-BROMO-TRYPTOPHAN
#sed -i 's/1TQ/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    1TQ TRP  6-(FORMYLAMINO)-7-HYDROXY-L-TRYPTOPHAN
#sed -i 's/23S/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    23S TRP  MODIFIED TRYPTOPHAN
#sed -i 's/32S/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    32S TRP  MODIFIED TRYPTOPHAN
#sed -i 's/32T/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    32T TRP  MODIFIED TRYPTOPHAN
#sed -i 's/4DP/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    4DP TRP
#sed -i 's/4FW/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    4FW TRP  4-FLUOROTRYPTOPHANE
#sed -i 's/4HT/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    4HT TRP  4-HYDROXYTRYPTOPHAN
#sed -i 's/4IN/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    4IN TRP  4-AMINO-L-TRYPTOPHAN
#sed -i 's/6CW/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    6CW TRP  6-CHLORO-L-TRYPTOPHAN
#sed -i 's/DTR/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    DTR TRP  D-TRYPTOPHAN
#sed -i 's/FTR/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    FTR TRP  FLUOROTRYPTOPHANE
#sed -i 's/HTR/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    HTR TRP  BETA-HYDROXYTRYPTOPHANE
#sed -i 's/PAT/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    PAT TRP  ALPHA-PHOSPHONO-TRYPTOPHAN
#sed -i 's/TOX/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    TOX TRP
#sed -i 's/TPL/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    TPL TRP  TRYTOPHANOL
#sed -i 's/TQQ/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    TQQ TRP
#sed -i 's/TRF/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    TRF TRP  N1-FORMYL-TRYPTOPHAN
#sed -i 's/TRN/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    TRN TRP  AZA-TRYPTOPHAN
#sed -i 's/TRO/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    TRO TRP  2-HYDROXY-TRYPTOPHAN
#sed -i 's/TRP/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    TRP TRP
#sed -i 's/TRQ/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    TRQ TRP
#sed -i 's/TRW/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    TRW TRP
#sed -i 's/TRX/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    TRX TRP  6-HYDROXYTRYPTOPHAN
#sed -i 's/TTQ/TRP/g' ${model%.*}_r_prep1.pdb                 ###                    TTQ TRP  6-AMINO-7-HYDROXY-L-TRYPTOPHAN
#sed -i 's/DTY/TYR/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/IYR/TYR/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/PAQ/TYR/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/PTR/TYR/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/STY/TYR/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/TYB/TYR/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/TYQ/TYR/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/TYR/TYR/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/TYS/TYR/g' ${model%.*}_r_prep1.pdb                 ###                  TYR
#sed -i 's/TYY/TYR/g' ${model%.*}_r_prep1.pdb                 ###                  HETEROATOM THAT MAY BE TREATED AS TYR
#sed -i 's/1TY/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    1TY TYR
#sed -i 's/2TY/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    2TY TYR
#sed -i 's/3TY/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    3TY TYR  MODIFIED TYROSINE
#sed -i 's/B3Y/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    B3Y TYR
#sed -i 's/CRQ/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    CRQ TYR
#sed -i 's/DBY/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    DBY TYR  35 DIBROMOTYROSINE
#sed -i 's/DPQ/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    DPQ TYR  TYROSINE DERIVATIVE
#sed -i 's/DTY/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    DTY TYR  D-TYROSINE
#sed -i 's/ESB/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    ESB TYR
#sed -i 's/FLT/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    FLT TYR  FLUOROMALONYL TYROSINE
#sed -i 's/FTY/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    FTY TYR  DEOXY-DIFLUOROMETHELENE-PHOSPHOTYROSINE
#sed -i 's/IYR/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    IYR TYR  3-IODO-TYROSINE
#sed -i 's/MBQ/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    MBQ TYR
#sed -i 's/NIY/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    NIY TYR  META-NITRO-TYROSINE
#sed -i 's/NBQ/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    NBQ TYR
#sed -i 's/OTY/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    OTY TYR
#sed -i 's/PAQ/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    PAQ TYR  SEE REMARK 999
#sed -i 's/PTH/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    PTH TYR  METHYLENE-HYDROXY-PHOSPHOTYROSINE
#sed -i 's/PTM/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    PTM TYR  ALPHA-METHYL-O-PHOSPHOTYROSINE
#sed -i 's/PTR/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    PTR TYR  O-PHOSPHOTYROSINE
#sed -i 's/TCQ/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    TCQ TYR  MODIFIED TYROSINE
#sed -i 's/TTS/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    TTS TYR
#sed -i 's/TY2/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    TY2 TYR  3-AMINO-L-TYROSINE
#sed -i 's/TY3/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    TY3 TYR  3-HYDROXY-L-TYROSINE
#sed -i 's/TYB/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    TYB TYR  TYROSINAL
#sed -i 's/TYC/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    TYC TYR  L-TYROSINAMIDE
#sed -i 's/TYI/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    TYI TYR  35-DIIODOTYROSINE
#sed -i 's/TYN/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    TYN TYR  ADDUCT AT HYDROXY GROUP
#sed -i 's/TYO/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    TYO TYR
#sed -i 's/TYQ/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    TYQ TYR  AMINOQUINOL FORM OF TOPA QUINONONE
#sed -i 's/TYR/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    TYR TYR
#sed -i 's/TYS/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    TYS TYR  INE SULPHONATED TYROSINE
#sed -i 's/TYT/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    TYT TYR
#sed -i 's/TYY/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    TYY TYR  IMINOQUINONE FORM OF TOPA QUINONONE
#sed -i 's/YOF/TYR/g' ${model%.*}_r_prep1.pdb                 ###                    YOF TYR  3-FLUOROTYROSINE


#################################################################################################################################################################



pdb4amber -i ${model%.*}_r_prep1.pdb -o ${model%.*}_r_prep11.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${model%.*}_r_prep11.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${model%.*}_r_prep2.prmtop ${model%.*}_r_prep2.inpcrd
quit
@


cp ${model%.*}_r_prep11.pdb ${model%.*}_r_prep00012.pdb


declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 5);
do
if [ ! -f ${model%.*}_r_prep2.inpcrd ] || [ ! -s ${model%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Paramter expansion susbtring extraction ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis000${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis000${count}.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis000${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
a=${A:1:3}
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${model%.*}_r_prep009.pdb >> ${model%.*}_r_prep010.pdb
sed -e "/${a}.*${B}.*H/d" ${model%.*}_r_prep000${pcount}.pdb >> ${model%.*}_r_prep000${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${model%.*}_r_prep000${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${model%.*}_r_prep2.prmtop ${model%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;

cp ${model%.*}_r_prep00017.pdb ${model%.*}_r_prep1112.pdb

declare -i pcount=12
declare -i p2count=13

for i in $(seq 35);
do

if [ ! -f ${model%.*}_r_prep2.inpcrd ] || [ ! -s ${model%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis11${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis11${count}.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis11${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

sed -e "/${A}.*${B}.*H/d" ${model%.*}_r_prep11${pcount}.pdb >> ${model%.*}_r_prep11${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${model%.*}_r_prep11${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${model%.*}_r_prep2.prmtop ${model%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :

fi;
done;
cp ${model%.*}_r_prep1145.pdb ${model%.*}_r_prep11112.pdb

declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 5);
do
if [ ! -f ${model%.*}_r_prep2.inpcrd ] || [ ! -s ${model%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Paramter expansion susbtring extraction ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis111${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis111${count}.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis111${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
a=${A:1:3}
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${model%.*}_r_prep009.pdb >> ${model%.*}_r_prep010.pdb
sed -e "/${a}.*${B}.*H/d" ${model%.*}_r_prep111${pcount}.pdb >> ${model%.*}_r_prep111${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${model%.*}_r_prep111${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${model%.*}_r_prep2.prmtop ${model%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;
cp ${model%.*}_r_prep11117.pdb ${model%.*}_r_prep0012.pdb

declare -i count=1
declare -i pcount=12
declare -i p2count=13

for i in $(seq 35);
do
if [ ! -f ${model%.*}_r_prep2.inpcrd ] || [ ! -s ${model%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - BACKBONE removal in progress ...";
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis00${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis00${count}.out);
del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis00${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"
C="$(cut -d' ' -f1 <<<$del2)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${model%.*}_r_prep00${pcount}.pdb >> ${model%.*}_r_prep00${p2count}.pdb
sed -e "/${C}.*${A}.*${B}/d" ${model%.*}_r_prep00${pcount}.pdb >> ${model%.*}_r_prep00${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${model%.*}_r_prep00${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${model%.*}_r_prep2.prmtop ${model%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;


if [ ! -f ${model%.*}_r_prep2.inpcrd ] || [ ! -s ${model%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - O1P removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis001.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis001.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis001.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${model%.*}_r_prep0047.pdb >> ${model%.*}_r_prep001.pdb
sed -e "/O1P.*${A}.*${B}/d" ${model%.*}_r_prep0047.pdb >> ${model%.*}_r_prep001.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${model%.*}_r_prep001.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${model%.*}_r_prep2.prmtop ${model%.*}_r_prep2.inpcrd
quit
@
fi;


if [ ! -f ${model%.*}_r_prep2.inpcrd ] || [ ! -s ${model%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - O2P removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis003.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis003.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis003.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${model%.*}_r_prep12.pdb >> ${model%.*}_r_prep002.pdb
sed -e "/O2P.*${A}.*${B}/d" ${model%.*}_r_prep001.pdb >> ${model%.*}_r_prep002.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${model%.*}_r_prep002.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${model%.*}_r_prep2.prmtop ${model%.*}_r_prep2.inpcrd
quit
@
fi;



if [ ! -f ${model%.*}_r_prep2.inpcrd ] || [ ! -s ${model%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - O3P removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis005.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis005.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis005.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${model%.*}_r_prep12.pdb >> ${model%.*}_r_prep003.pdb
sed -e "/O3P.*${A}.*${B}/d" ${model%.*}_r_prep002.pdb >> ${model%.*}_r_prep003.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${model%.*}_r_prep003.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${model%.*}_r_prep2.prmtop ${model%.*}_r_prep2.inpcrd
quit
@
fi;



if [ ! -f ${model%.*}_r_prep2.inpcrd ] || [ ! -s ${model%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis007.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis007.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis007.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${model%.*}_r_prep003.pdb >> ${model%.*}_r_prep004.pdb
sed -e "/${A}.*${B}.*P/d" ${model%.*}_r_prep003.pdb >> ${model%.*}_r_prep004.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${model%.*}_r_prep004.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${model%.*}_r_prep2.prmtop ${model%.*}_r_prep2.inpcrd
quit
@
fi;


if [ ! -f ${model%.*}_r_prep2.inpcrd ] || [ ! -s ${model%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis009.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis009.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis009.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${model%.*}_r_prep004.pdb >> ${model%.*}_r_prep005.pdb
sed -e "/${A}.*${B}.*P/d" ${model%.*}_r_prep004.pdb >> ${model%.*}_r_prep005.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${model%.*}_r_prep005.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${model%.*}_r_prep2.prmtop ${model%.*}_r_prep2.inpcrd
quit
@
fi;


if [ ! -f ${model%.*}_r_prep2.inpcrd ] || [ ! -s ${model%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis011.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis011.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis011.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/P/d'  ${model%.*}_r_prep005.pdb >> ${model%.*}_r_prep006.pdb
sed -e "/${A}.*${B}.*P/d" ${model%.*}_r_prep005.pdb >> ${model%.*}_r_prep006.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${model%.*}_r_prep006.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${model%.*}_r_prep2.prmtop ${model%.*}_r_prep2.inpcrd
quit
@
fi;


if [ ! -f ${model%.*}_r_prep2.inpcrd ] || [ ! -s ${model%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis013.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis013.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis013.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${model%.*}_r_prep006.pdb >> ${model%.*}_r_prep007.pdb
sed -e "/${A}.*${B}.*H/d" ${model%.*}_r_prep006.pdb >> ${model%.*}_r_prep007.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${model%.*}_r_prep007.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${model%.*}_r_prep2.prmtop ${model%.*}_r_prep2.inpcrd
quit
@
fi;


if [ ! -f ${model%.*}_r_prep2.inpcrd ] || [ ! -s ${model%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis015.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis015.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis015.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${model%.*}_r_prep007.pdb >> ${model%.*}_r_prep008.pdb
sed -e "/${A}.*${B}.*H/d" ${model%.*}_r_prep007.pdb >> ${model%.*}_r_prep008.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${model%.*}_r_prep008.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${model%.*}_r_prep2.prmtop ${model%.*}_r_prep2.inpcrd
quit
@
fi;


if [ ! -f ${model%.*}_r_prep2.inpcrd ] || [ ! -s ${model%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Hyd removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis017.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis017.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis017.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${model%.*}_r_prep008.pdb >> ${model%.*}_r_prep009.pdb
sed -e "/${A}.*${B}.*H/d" ${model%.*}_r_prep008.pdb >> ${model%.*}_r_prep009.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${model%.*}_r_prep009.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${model%.*}_r_prep2.prmtop ${model%.*}_r_prep2.inpcrd
quit
@
fi;


if [ ! -f ${model%.*}_r_prep2.inpcrd ] || [ ! -s ${model%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Paramter expansion susbtring extraction ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis019.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis019.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis019.out);
A="$(cut -d' ' -f1 <<<$del)"
a=${A:1:3}
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${model%.*}_r_prep009.pdb >> ${model%.*}_r_prep010.pdb
sed -e "/${a}.*${B}.*H/d" ${model%.*}_r_prep009.pdb >> ${model%.*}_r_prep010.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${model%.*}_r_prep010.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${model%.*}_r_prep2.prmtop ${model%.*}_r_prep2.inpcrd
quit
@
fi;


if [ ! -f ${model%.*}_r_prep2.inpcrd ] || [ ! -s ${model%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - HIS removal in progress ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis021.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis021.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis021.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${model%.*}_r_prep010.pdb >> ${model%.*}_r_prep011.pdb
sed -e "/HIS.*${B}.*H/d" ${model%.*}_r_prep010.pdb >> ${model%.*}_r_prep011.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${model%.*}_r_prep011.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${model%.*}_r_prep2.prmtop ${model%.*}_r_prep2.inpcrd
quit
@
fi;


declare -i count=23
declare -i pcount=11
declare -i p2count=12

for i in $(seq 35);
do
if [ ! -f ${model%.*}_r_prep2.inpcrd ] || [ ! -s ${model%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - BACKBONE removal in progress ...";
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis0${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis0${count}.out);
del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis0${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
B="$(cut -d' ' -f2 <<<$del)"
C="$(cut -d' ' -f1 <<<$del2)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${model%.*}_r_prep0${pcount}.pdb >> ${model%.*}_r_prep0${p2count}.pdb
sed -e "/${C}.*${A}.*${B}/d" ${model%.*}_r_prep0${pcount}.pdb >> ${model%.*}_r_prep0${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${model%.*}_r_prep0${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${model%.*}_r_prep2.prmtop ${model%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
	:
fi;
done;

if [ ! -f ${model%.*}_r_prep2.inpcrd ] || [ ! -s ${model%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - P, O1P, O2P, O3P removal in progress ..."

sed -i '/ATOM.*\<P\>/d' ${model%.*}_r_prep046.pdb
sed -i '/ATOM.*\<O1P\>/d' ${model%.*}_r_prep046.pdb
sed -i '/ATOM.*\<O2P\>/d' ${model%.*}_r_prep046.pdb
sed -i '/ATOM.*\<O3P\>/d' ${model%.*}_r_prep046.pdb



tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${model%.*}_r_prep046.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${model%.*}_r_prep2.prmtop ${model%.*}_r_prep2.inpcrd
quit
@
fi;

declare -i count=1
declare -i pcount=46
declare -i p2count=47

for i in $(seq 8);
do
if [ ! -f ${model%.*}_r_prep2.inpcrd ] || [ ! -s ${model%.*}_r_prep2.inpcrd ]; then
echo "Incompatibility detected with the amino acid library in ff14SB - Paramter expansion susbtring extraction ..."
grep "^FATAL" ./slurm*.out | tail -1 > diagnosis22${count}.out;
del=$(awk -F"[<>]" 'NR==1 {print $2}' diagnosis22${count}.out);
#del2=$(awk -F"[<>]" 'NR==1 {print $4}' diagnosis22${count}.out);
A="$(cut -d' ' -f1 <<<$del)"
a=${A:1:3}
B="$(cut -d' ' -f2 <<<$del)"

#sed -e '/${A}/!b' -e '/${B}/!b' -e '/H/d'  ${model%.*}_r_prep009.pdb >> ${model%.*}_r_prep010.pdb
sed -e "/${a}.*${B}.*H/d" ${model%.*}_r_prep0${pcount}.pdb >> ${model%.*}_r_prep0${p2count}.pdb

tleap -f - <<@
source leaprc.protein.ff14SB         # aminoacid dictionary
source leaprc.gaff                   # small molecule atom types
source leaprc.water.tip3p            # ions and water
loadamberprep ${collective%.*}_ligand_oprot_prep1_gas_bcc.prepi    # load cofactor parameters
loadamberparams ${collective%.*}_ligand_oprot_prep1_gas_bcc.frcmod
m1 = loadpdb ${model%.*}_r_prep0${p2count}.pdb        # load receptor with ion and cofactor
check m1
charge m1
solvatebox m1 TIP3PBOX 8 #Solvate the complex with a cubic water box
addions m1 Cl- 0 #Add Cl- ions to neutralize the system
saveamberparm m1 ${model%.*}_r_prep2.prmtop ${model%.*}_r_prep2.inpcrd
quit
@

let count++
let pcount++
let p2count++

else
        :
fi;
done;



ambpdb -p ${model%.*}_r_prep2.prmtop -c ${model%.*}_r_prep2.inpcrd -mol2 -sybyl > ${model%.*}_r_prep3.mol2

done;
