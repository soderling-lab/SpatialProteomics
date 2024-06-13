import pandas as pd
import os
import numpy as np
import sys
import argparse

cwd = os.getcwd()
parent = os.path.dirname(cwd)
gparent = os.path.dirname(parent)

def is_first_char_digit(s):
    return s[0].isdigit() if s else False

def ensure_dirs_exists(path):
    if "." in path:
        path = os.path.dirname(path)
    if not os.path.exists(path):
        os.makedirs(path)
    return
# 112744 Fxn10_Control Normalized
#### NEED TO CHANGE FOR EACH DIFF DATASET
def scrapeGenotype(sample_name):
    s1 = sample_name.split()[1]
    s2 = s1.split(f'_')[-1]
    return s2

def scrapeBioFraction(sample_name):
    s1 = sample_name.split()[1]
    s2 = s1.split(f'_')[0]
    return s2

infile = f'~/Documents/SoderlingLab/ProjectHELA/10415/KinSub10415_unbiased_normalized_041124_AnnotatedHumanKinome.csv'
ORGANISM = F'Homo sapiens'
data = pd.read_csv(infile, low_memory=False)
## CV is control not needed. 

to_drop = [x for x in data.columns if any(sub in x.lower() for sub in ['cv', 'spqc'])]

data = data.drop(columns=to_drop)
print(f'data columns you are keeping {data.columns}')
## DATA COLS START W NUM. SO FOLLOWING FILTER APPLIED.
pdata_cols = [x for x in data.columns if is_first_char_digit(x) is True] 

tdata = pd.melt(data, id_vars=['Accession', 'Description', 'Human_Kinase', 'Genes', 'Detected_Imputed'],
        value_vars=pdata_cols, var_name='Mixture', value_name='Intensity')


## Intensity is raw vals 
## Required: Total_Intensity, Rel_Intensity, Abunadnce
## total intensity is sum of accession across ALL (mut,tg) samples given. 
### Abundance log2(Rel_Intensity)
tdata['Total_Intensity'] = (tdata.groupby('Accession')['Intensity']).transform('sum')
tdata['Relative_Intensity'] = tdata['Intensity']/tdata['Total_Intensity']
tdata['Abundance'] = np.log2(tdata['Intensity'])

cols2change = {'Accession': 'Protein', 'Description': 'Function',
              'Genes': 'Gene'}
tdata = tdata.rename(columns = cols2change)

cols2check = {'Origin': f'{ORGANISM}', 'Genotype': scrapeGenotype, 
              'BioFraction': scrapeBioFraction}
for x, v in cols2check.items():
    if x not in tdata.columns:
        if callable(v):
            tdata[x] = tdata['Mixture'].apply(v)
        else:
            tdata[x] = v

lk = f'{cwd}/transformeddata/{os.path.split(os.path.dirname(infile))[-1]}/Transformed_{os.path.basename(infile)}'
ensure_dirs_exists(lk)
print('saved to: ', lk)
tdata.to_csv(lk)


## Spatial proteomics needs verbatim cols: Protein, Function, Gene, Origin, Mixture, Genotype, Biofraction