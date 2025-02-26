import pandas as pd
import os
import numpy as np
import sys
import argparse
import re 

cwd = os.getcwd()
parent = os.path.dirname(cwd)
gparent = os.path.dirname(parent)


def main():
    infile = f'/home/poojaparameswaran/Documents/SoderlingLab/SpatialProteomics/rawdata/LOPIT_SNCA_young.xlsx'
    ORGANISM = F'Mus Musculus'
    MIXTURES = True
    if 'csv' in infile:
        data = pd.read_csv(infile, low_memory=False)
    elif 'xlsx' in infile:
        file = pd.ExcelFile(infile)
        SHEETNAME= next(x for x in file.sheet_names if 'normalized' in x.lower())
        data = pd.read_excel(infile, sheet_name=SHEETNAME)
    ## CV is control not needed.
    to_drop = [x for x in data.columns if any(sub in x.lower() for sub in ['cv', 'spqc'])]
    data = data.drop(columns=to_drop)
    data = data[data['Description'].apply(is_mouse, args=(ORGANISM,))]
    ## DATA COLS START W NUM. SO FOLLOWING FILTER APPLIED.
    pdata_cols = [x for x in data.columns if is_first_char_digit(x) is True]
    print(f'Data will be elongated by a factor of {len(pdata_cols)} since each condition is exploded. {data["Accession"].nunique()} proteins')
    tdata = pd.melt(data, id_vars=[x for x in data.columns if not is_first_char_digit(x)],
            value_vars=pdata_cols, var_name='Identity', value_name='Intensity')
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
                tdata[x] = tdata['Identity'].apply(v)
            else:
                tdata[x] = v
    if MIXTURES:
        mixdict = {}
        mixtures_mapped = extract_mixture(pdata_cols)
        tdata[f'Mixture'] = tdata['Identity'].str.extract(r'(\S+)\s+')[0].astype(int).map(mixtures_mapped)
    tdata = tdata.drop(columns=[x for x in tdata.columns if x not in ['Protein', 'Function', 'Identity', 'Intensity', 
                                                                      'Total_Intensity', 'Relative_Intensity', 'BioFraction',
                                                                      'Genotype', 'Mixture',
                                                                      'Abundance']])
    ## the raw data is labeled under Intensity.
    tdata['Abundance'].fillna(0.001, inplace=True)
    lk = f'/home/poojaparameswaran/Documents/SoderlingLab/SpatialProteomics/transformeddata/Transformed_{os.path.basename(infile).split(".")[0]}.csv'
    ensure_dirs_exists(lk)
    print(f'TransformedData Shape {tdata.shape}, rows should be {len(pdata_cols) * data["Accession"].nunique()} (cross product of condition cols x protein count)')
    tdata.to_csv(lk, index=False)
    print('saved to: ', lk) 
    ## Spatial proteomics needs verbatim cols: Protein, Function, Gene, Origin, Mixture, Genotype, Biofraction
    return

def extract_mixture(data_cols):
    mixture_ranges = [int(re.search(r'(\S+)\s+', x).group(1)) for x in data_cols]
    d = {x: f'M{i//14}' for i, x in enumerate(mixture_ranges)}
    return d

def is_mouse(description, organism_to_pick):
    org = re.search('OS=(.*?)\sOX=', description).group(1) if re.search('OS=(.*?)\sOX=', description) \
         else ''
    if all(x in org.lower() for x in organism_to_pick.lower().split()):
        return True
    return False

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
def scrapeGenotype(sample_name): # '73724 wt_F7 Normalized'
    s1 = re.search('^\S+\s+(\S+)_.+$', sample_name).group(1)
    return s1

def scrapeBioFraction(sample_name):
    s1 = re.search(r'^\S+\s+\S+_(\S+).+$', sample_name).group(1)
    return s1


if __name__ == '__main__':
    main()