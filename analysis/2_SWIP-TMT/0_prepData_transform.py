#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 16:12:53 2023

@author: poojap
"""

import pandas as pd
import os
import numpy as np

def main():
    ## SNCA == Normalized Data, LRRK2 == Table 4 Normalized Data
    ## Vps35 == Table 3 Normalized Data
    # filename = '/Users/poojap/Documents/SoderlingLab/SwipProteomics-master copy/data/swip_tmt copy.csv'
    filename = 'LOPIT_LRRK2_young.xlsx'
    worksheet = 'Table 4 Normalized Data'
    
    currdir = os.getcwd()
    parent = os.path.dirname(currdir)
    gparent = os.path.dirname(os.path.dirname(currdir))
    rawdatafile = f"{gparent}/rawdata/{filename}"
    # do conversion
    # dataframe = convert_XLSX_file(rawdatafile, worksheet)
    dataframe = transform_raw(rawdatafile, worksheet)

    # write to file
    filename = filename.replace(".xlsx", "")
    newfilepath =  os.path.join(gparent, "transformeddata", filename)
    writeToFile(newfilepath, dataframe)


def transform_raw(filename, sheetname):
    if "xlsx" in filename:
        df = pd.read_excel(filename, sheet_name= sheetname)
    if "csv" in filename:
        df = pd.read_csv(filename)
    pd.set_option('display.max_rows', None)
    list_of_all = []
    
    M1 = list(range(3, 17))
    M2 = list(range(19, 33))
    M3 = list(range(35, 49))
    
    for i in range(df.shape[1]): # i is cols
        if i in M1:
            if "wt" in df.columns[i]:
                newdict = {'Protein': list(df['Accession']),
                            'Gene': list(df["Description"]),
                            'Function': list(df["Description"]),
                            'Origin': list(df["Description"]),
                            'Mixture': (['M1'] * len(df['Accession'])),
                            'Genotype': (['WildType']*len(df['Accession'])),
                            'BioFraction': None,
                            'Intensity': list(df.iloc[:, i]),
                            'Abundance': np.log2(list(df.iloc[:, i])),
                            'Rel_Intensity' :list(df.groupby('Accession').apply(lambda x: x.iloc[:, M1 +M2+M3].sum(axis=1)))
                            }
                for ind, val in enumerate(newdict["Rel_Intensity"]):
                    newdict["Rel_Intensity"][ind] = newdict["Intensity"][ind]/newdict["Rel_Intensity"][ind]

                for ind, val in enumerate(newdict["Gene"]):
                    genestr = "GN"
                    endstr = " "
                    startind = val.find(genestr) + len(genestr)
                    endind = val.find(endstr, startind)
                    genename = val[startind+1: endind]
                    newdict["Gene"][ind] = genename
                    
                    osstr = "OS"
                    endstr = "OX"
                    startind = val.rfind(osstr) + len(osstr)
                    endind = val.rfind(endstr, startind)
                    organism_species = val[startind+1: endind]
                    newdict["Origin"][ind] = organism_species
                    if endstr not in val:
                        endstr = "GN"
                        startind = val.rfind(osstr) + len(osstr)
                        endind = val.rfind(endstr, startind)
                        organism_species = val[startind+1: endind]
                        newdict["Origin"][ind] = organism_species
                    endstr = "OS"
                    endind = val.rfind(endstr, 0)
                    function = val[0: endind]
                    newdict["Function"][ind] = function
                    
                string = df.columns[i].lower()
                index = string.find('wt_f')
                if (index != -1) and index + len("wt_f") < len(string):
                    fraction_num = string[index + len("wt_f"): index + len("wt_f") + 2]
                    fraction = "F" + str(fraction_num)
                    dict1 = {'BioFraction': [fraction] *len(df['Accession'])}
                    newdict.update(dict1)
                list_of_all.append(newdict)
            
            elif "tg" in df.columns[i]:
                newdict = {'Protein': list(df['Accession']),
                            'Gene': list(df["Description"]),
                            'Function': list(df["Description"]),
                            'Origin': list(df["Description"]),
                            'Mixture': (['M1'] * len(df['Accession'])),
                            'Genotype': (['Transgenic']*len(df['Accession'])),
                            'BioFraction': None,
                            'Intensity': list(df.iloc[:, i]),
                            'Abundance': np.log2(list(df.iloc[:, i])),
                            'Rel_Intensity' :list(df.groupby('Accession').apply(lambda x: x.iloc[:, M1 +M2+M3].sum(axis=1)))
                            }
                for ind, val in enumerate(newdict["Rel_Intensity"]):
                    newdict["Rel_Intensity"][ind] = newdict["Intensity"][ind]/newdict["Rel_Intensity"][ind]
                    
                for ind, val in enumerate(newdict["Gene"]):
                    genestr = "GN"
                    endstr = " "
                    startind = val.find(genestr) + len(genestr)
                    endind = val.find(endstr, startind)
                    genename = val[startind+1: endind]
                    newdict["Gene"][ind] = genename
                    
                    osstr = "OS"
                    endstr = "OX"
                    startind = val.rfind(osstr) + len(osstr)
                    endind = val.rfind(endstr, startind)
                    organism_species = val[startind+1: endind]
                    newdict["Origin"][ind] = organism_species
                    if endstr not in val:
                        endstr = "GN"
                        startind = val.rfind(osstr) + len(osstr)
                        endind = val.rfind(endstr, startind)
                        organism_species = val[startind+1: endind]
                        newdict["Origin"][ind] = organism_species
                    
                    endstr = "OS"
                    endind = val.rfind(endstr, 0)
                    function = val[0: endind]
                    newdict["Function"][ind] = function


                string = df.columns[i].lower()
                index = string.find('tg_f')
                if (index != -1) and index + len("tg_f") < len(string):
                    fraction_num = string[index + len("tg_f"): index + len("tg_f") + 2]
                    fraction = "F" + str(fraction_num)
                    dict1 = {'BioFraction': [fraction] *len(df['Accession'])}
                    newdict.update(dict1)
                list_of_all.append(newdict)
                
        if i in M2:
            if "wt" in df.columns[i]:
                newdict = {'Protein': list(df['Accession']),
                            'Gene': list(df["Description"]),
                            'Function': list(df["Description"]),
                            'Origin': list(df["Description"]),
                            'Mixture': (['M2'] * len(df['Accession'])),
                            'Genotype': (['WildType']*len(df['Accession'])),
                            'BioFraction': None,
                            'Intensity': list(df.iloc[:, i]),
                            'Abundance': np.log2(list(df.iloc[:, i])),
                            'Rel_Intensity' :list(df.groupby('Accession').apply(lambda x: x.iloc[:, M1 +M2+M3].sum(axis=1)))
                            }
                for ind, val in enumerate(newdict["Rel_Intensity"]):
                    newdict["Rel_Intensity"][ind] = newdict["Intensity"][ind]/newdict["Rel_Intensity"][ind]
                
                for ind, val in enumerate(newdict["Gene"]):
                    genestr = "GN"
                    endstr = " "
                    startind = val.find(genestr) + len(genestr)
                    endind = val.find(endstr, startind)
                    genename = val[startind+1: endind]
                    newdict["Gene"][ind] = genename
                    
                    osstr = "OS"
                    endstr = "OX"
                    startind = val.rfind(osstr) + len(osstr)
                    endind = val.rfind(endstr, startind)
                    organism_species = val[startind+1: endind]
                    newdict["Origin"][ind] = organism_species
                    if endstr not in val:
                        endstr = "GN"
                        startind = val.rfind(osstr) + len(osstr)
                        endind = val.rfind(endstr, startind)
                        organism_species = val[startind+1: endind]
                        newdict["Origin"][ind] = organism_species
                    
                    endstr = "OS"
                    endind = val.rfind(endstr, 0)
                    function = val[0: endind]
                    newdict["Function"][ind] = function
                string = df.columns[i].lower()
                index = string.find('wt_f')
                if (index != -1) and index + len("wt_f") < len(string):
                    fraction_num = string[index + len("wt_f"): index + len("wt_f") + 2]
                    fraction = "F" + str(fraction_num)
                    newdict["BioFraction"] = [fraction] *len(df['Accession'])
                list_of_all.append(newdict)
            elif "tg" in df.columns[i]:
                newdict = {'Protein': list(df['Accession']),
                            'Gene': list(df["Description"]),
                            'Function': list(df["Description"]),
                            'Origin': list(df["Description"]),
                            'Mixture': (['M2'] * len(df['Accession'])),
                            'Genotype': (['Transgenic']*len(df['Accession'])),
                            'BioFraction': None,
                            'Intensity': list(df.iloc[:, i]),
                            'Abundance': np.log2(list(df.iloc[:, i])),
                            'Rel_Intensity' :list(df.groupby('Accession').apply(lambda x: x.iloc[:, M1 +M2+M3].sum(axis=1)))
                            }
                for ind, val in enumerate(newdict["Rel_Intensity"]):
                    newdict["Rel_Intensity"][ind] = newdict["Intensity"][ind]/newdict["Rel_Intensity"][ind]

                for ind, val in enumerate(newdict["Gene"]):
                    genestr = "GN"
                    endstr = " "
                    startind = val.find(genestr) + len(genestr)
                    endind = val.find(endstr, startind)
                    genename = val[startind+1: endind]
                    newdict["Gene"][ind] = genename
                    
                    osstr = "OS"
                    endstr = "OX"
                    startind = val.rfind(osstr) + len(osstr)
                    endind = val.rfind(endstr, startind)
                    organism_species = val[startind+1: endind]
                    newdict["Origin"][ind] = organism_species
                    if endstr not in val:
                        endstr = "GN"
                        startind = val.rfind(osstr) + len(osstr)
                        endind = val.rfind(endstr, startind)
                        organism_species = val[startind+1: endind]
                        newdict["Origin"][ind] = organism_species
                    
                    endstr = "OS"
                    endind = val.rfind(endstr, 0)
                    function = val[0: endind]
                    newdict["Function"][ind] = function
    
                string = df.columns[i].lower()
                index = string.find('tg_f')
                if (index != -1) and index + len("tg_f") < len(string):
                    fraction_num = string[index + len("tg_f"): index + len("tg_f") + 2]
                    fraction = "F" + str(fraction_num)
                    newdict["BioFraction"] = [fraction] *len(df['Accession'])
                list_of_all.append(newdict)
    for i in range(df.shape[1]):
        if i in M3:
            if "wt" in df.columns[i]:
                newdict = {'Protein': list(df['Accession']),
                            'Gene': list(df["Description"]),
                            'Function': list(df["Description"]),
                            'Origin': list(df["Description"]),
                            'Mixture': (['M3'] * len(df['Accession'])),
                            'Genotype': (['WildType']*len(df['Accession'])),
                            'BioFraction': None,
                            'Intensity': list(df.iloc[:, i]),
                            'Abundance': np.log2(list(df.iloc[:, i])),
                            'Rel_Intensity' :list(df.groupby('Accession').apply(lambda x: x.iloc[:, M1 +M2+M3].sum(axis=1)))
                            }
                for ind, val in enumerate(newdict["Rel_Intensity"]):
                    newdict["Rel_Intensity"][ind] = newdict["Intensity"][ind]/newdict["Rel_Intensity"][ind]

                for ind, val in enumerate(newdict["Gene"]):
                    genestr = "GN"
                    endstr = " "
                    startind = val.find(genestr) + len(genestr)
                    endind = val.find(endstr, startind)
                    genename = val[startind+1: endind]
                    newdict["Gene"][ind] = genename
                    
                    osstr = "OS"
                    endstr = "OX"
                    startind = val.rfind(osstr) + len(osstr)
                    endind = val.rfind(endstr, startind)
                    organism_species = val[startind+1: endind]
                    newdict["Origin"][ind] = organism_species
                    if endstr not in val:
                        endstr = "GN"
                        startind = val.rfind(osstr) + len(osstr)
                        endind = val.rfind(endstr, startind)
                        organism_species = val[startind+1: endind]
                        newdict["Origin"][ind] = organism_species
                    
                    endstr = "OS"
                    endind = val.rfind(endstr, 0)
                    function = val[0: endind]
                    newdict["Function"][ind] = function
                            
                string = df.columns[i].lower()
                index = string.find('wt_f')
                if (index != -1) and index + len("wt_f") < len(string):
                    fraction_num = string[index + len("wt_f"): index + len("wt_f") + 2]
                    fraction = "F" + str(fraction_num)
                    dict1 = {'BioFraction': [fraction] *len(df['Accession'])}
                    newdict.update(dict1)
                list_of_all.append(newdict)
            
            elif "tg" in df.columns[i]:
                newdict = {'Protein': list(df['Accession']),
                            'Gene': list(df["Description"]),
                            'Function': list(df["Description"]),
                            'Origin': list(df["Description"]),
                            'Mixture': (['M3'] * len(df['Accession'])),
                            'Genotype': (['Transgenic']*len(df['Accession'])),
                            'BioFraction': None,
                            'Intensity': list(df.iloc[:, i]),
                            'Abundance': np.log2(list(df.iloc[:, i])),
                            'Rel_Intensity' :list(df.groupby('Accession').apply(lambda x: x.iloc[:, M1 +M2+M3].sum(axis=1)))
                            }
                for ind, val in enumerate(newdict["Rel_Intensity"]):
                    newdict["Rel_Intensity"][ind] = newdict["Intensity"][ind]/newdict["Rel_Intensity"][ind]

                for ind, val in enumerate(newdict["Gene"]):
                    genestr = "GN"
                    endstr = " "
                    startind = val.find(genestr) + len(genestr)
                    endind = val.find(endstr, startind)
                    genename = val[startind+1: endind]
                    newdict["Gene"][ind] = genename
                    
                    osstr = "OS"
                    endstr = "OX"
                    startind = val.rfind(osstr) + len(osstr)
                    endind = val.rfind(endstr, startind)
                    organism_species = val[startind+1: endind]
                    newdict["Origin"][ind] = organism_species
                    if endstr not in val:
                        endstr = "GN"
                        startind = val.rfind(osstr) + len(osstr)
                        endind = val.rfind(endstr, startind)
                        organism_species = val[startind+1: endind]
                        newdict["Origin"][ind] = organism_species
                    
                    endstr = "OS"
                    endind = val.rfind(endstr, 0)
                    function = val[0: endind]
                    newdict["Function"][ind] = function

                string = df.columns[i].lower()
                index = string.find('tg_f')
                if (index != -1) and index + len("tg_f") < len(string):
                    fraction_num = string[index + len("tg_f"): index + len("tg_f") + 2]
                    fraction = "F" + str(fraction_num)
                    dict1 = {'BioFraction': [fraction] *len(df['Accession'])}
                    newdict.update(dict1)
                list_of_all.append(newdict)
    newdf = pd.DataFrame(list_of_all)
    # newdf = newdf.set_index()
    # print(newdf.head())
    newdf = newdf.apply(lambda x: pd.Series(x).explode())
    newdf = newdf.reset_index(drop = True)
    newdf.to_csv('.csv', index=True)
    # result = newdf[newdf['Genotype'] == 'Transgenic']
    print(newdf.shape)
    print(newdf.iloc[: :-4].tail())
    return newdf


def writeToFile(filename, dataframe):
    write = filename + '_Transformed.csv'
    dataframe.to_csv(write, index=True)
    print("file saved to: {}".format(write))
    return

if __name__ == "__main__":
    main()
