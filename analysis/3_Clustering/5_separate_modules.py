import pandas as pd
import numpy as np
import os 
import IDmapping

paritioned_data = f'/home/poojaparameswaran/Documents/SoderlingLab/SpatialProteomics/data/Annotated_partitionSigned.csv'
kinase_data = pd.read_csv('/home/poojaparameswaran/Documents/generalData/data/pkinfam.csv')
proteomicsdata = pd.read_excel('/home/poojaparameswaran/Documents/SoderlingLab/iltp/SelectedGenes_Before_FCOnly_10595_SupplementalData_082324subSeq.xlsx',
                                )
proteomicsdata.insert(proteomicsdata.columns.get_loc('Accession')+1, 'ID',
                      proteomicsdata['Accession']+'_'+proteomicsdata['SubstrateSeq'])

data = pd.read_csv(paritioned_data)
# data.rename(columns={data.columns[0]: 'ID'}, inplace=True)
# data = data.merge(proteomicsdata, on='ID', how='outer')
# ## get line of best fit
# def get_polynomialfit(row):
#     FCs = [x for x in row.index if 'foldchange' in str(x).lower()]
#     sl = 0
#     if FCs:
#         sl, inte = np.polyfit(list(range(0, len(FCs))), list(row[FCs].values), 1) ## simple linear fit to get trend
#     return sl

# data['regulation'] = data.apply(get_polynomialfit, axis=1)
# data.drop(columns=['Accession', 'Description', ''Genes', 'Stripped_Seq', 'Detected_Imputed',
                #    ] , inplace=True) # + [x for x in data.columns if 'foldchange' in str(x).lower()]
## Pushpa data has 230 peptides that have variance close to 0, so discarded.
accessions = data.iloc[:,0].str.split('_').str[0].values
d = IDmapping.main(accessions, inputdb='UniProtKB_AC-ID', outputdb='Gene_Name')
## annotate Kinases
print('####################################')
data.insert(1, 'Gene Names', '')
data['Gene Names'] = data.iloc[:,0].str.split(';').str[0].str.split('_').str[0].map(d)
data.insert(2, 'Kinase', 'KINASE')
ORGANISM = 'Mouse'
kinase_dictionary = {acc: (False if acc.split(';')[0] not in kinase_data[f'{ORGANISM}ID'].values
                    else 'KINASE') for acc in accessions }
data['Kinase'] = data.iloc[:, 0].str.split(';').str[0].str.split('_').str[0].map(kinase_dictionary)
data.rename(columns={data.columns[0]: 'Protein'}, inplace=True)

out = data.groupby('Module')
with pd.ExcelWriter(f'/home/poojaparameswaran/Documents/SoderlingLab/SpatialProteomics/tables/'
                    f'{os.path.basename(paritioned_data).split(".")[0]}.xlsx') as file:
    for x, d in out:
        d.to_excel(file, sheet_name=f'Module{x}', index=False)

print(f'written to /home/poojaparameswaran/Documents/SoderlingLab/SpatialProteomics/tables/'
    f'{os.path.basename(paritioned_data).split(".")[0]}.xlsx')

