#!/usr/bin/env python
import sys
import pandas as pd
import re
from pandas import DataFrame
import numpy

def lFactorial(val):
    returnValue = 0
    i = 2
    while i <= val:
        returnValue = returnValue + numpy.log(i)
        i += 1
    return returnValue
def lNchooseK(n, k):
    answer = 0
    if k > (n-k):
        k = n-k
    i = n
    while i > (n - k):
        answer = answer + numpy.log(i)
        i -= 1
    answer = answer - lFactorial(k)
    return answer
def hypergeometric(n, p, k, r):
    """
    extraido del script de perl de GeneMerge
    n = total poblacion
    p = parte de la poblacion con un item especifico
    k = total de la muestra
    r = parte de la muestra con el mismo item que p
    n, p, k, r = 6157, 222, 473, 81
    hypergeometric(n, p, k, r)
    """
    np = p
    nq = n - p
    p = p/n
    log_n_choose_k = lNchooseK(n, k)
    
    top = k
    if np < k:
        top = np
    lfoo = lNchooseK(np, top) + lNchooseK(n*(1-p), k-top)
    suma = 0
    i = top
    while i >= r:
        suma = suma + numpy.exp(lfoo - log_n_choose_k)
        if i > r:
            lfoo = lfoo + numpy.log(i / (np-i+1)) +  numpy.log((nq - k + i) / (k-i+1))
        i -= 1
    if suma > 1:
        suma = 1
    return suma



# open files
association=pd.read_csv('data/Association.txt',sep='\t', names = ['Entry', 'base'])
association['Entry'] = [str(i) for i in association.Entry]

description = pd.read_csv('../NeVOmics_DataBase/'+sys.argv[1],sep='\t', names = ['base' , 'Term'])

background = pd.read_csv('data/Background.txt',sep='\t', names = ['Entry'])
background['Entry'] = [str(i) for i in background.Entry]

List = pd.read_csv('data/List.txt',sep='\t', names = ['Entry'])
List['Entry'] = [str(i) for i in List.Entry]

# Total of proteins with terms in list (for hypergeometric dustribution)
total_protein_list = len(set(List.Entry.tolist()))
#total_protein_list

# Get all list terms involved in a category
list_cat = List.merge(association , on = 'Entry' , how = 'left').dropna().drop_duplicates()
yyy = list_cat.merge(description , on = 'base', how='left').dropna().drop_duplicates()
list_category = DataFrame(yyy.groupby('base').Entry.size()).reset_index()

# Total of proteins with terms in background
total_proteins_bg = len(set(background.Entry.tolist()))#.drop_duplicates().count()

# Get all background terms involved in a category
category = background.merge(association , on = 'Entry', how = 'left').reset_index(drop=True).dropna()
xxx = category.merge(description , on = 'base', how = 'left').dropna().drop_duplicates()
cat = DataFrame(xxx.groupby('base').Entry.count()).reset_index()

#######################################################################
#########################  Statistical test   #########################
#######################################################################
## list_count = proteins in list
## back_count = proteins in background by category (P or F or C or kegg)
statistics = list_category.merge(cat, on = 'base', how = 'left')
statistics.columns = ['base','list_count', 'back_count']



## following the GeneMerge1.4 approach we use non-singletons as multiple testing value
multiple_testing_value = statistics[statistics.back_count > 1].count()[0]

#print(multiple_testing_value)
statistics['tot_list']=total_protein_list
statistics['tot_back']=total_proteins_bg

## Hypergeometric Distribution
#hypergeom.sf(k, M, n, N, loc=0)
# k = number of genes/proteins associated to the process "cell cycle"
# M = total number of genes/proteins with some annotation
# n = total number of genes/proteins annotated for "cell cycle" inside M
# N = number of genes associated to at least one Biological Process in the Gene Ontology.
# https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.stats.hypergeom.html
# example =  hypergeom.sf(8-1, total_proteins_bg, 199, total_protein_list, loc=0)



## Loop for calculate hypergeometric distribution
p_val=[]
for index, row in statistics.iterrows():
    # from scipy.stats import hypergeom
    #b=hypergeom.sf(row['list_count']-1, total_proteins_bg, row['back_count'], total_protein_list, loc=0)
    b = hypergeometric(total_proteins_bg, row['back_count'], total_protein_list, row['list_count'])
    p_val.append(b)
statistics['P'] = p_val


## Loop for calculate Bonferroni correction
Bonf_cor=[]
for x in statistics.P:
    Bonf_cor.append(x*multiple_testing_value)
statistics['Bonf_corr'] = Bonf_cor


## Loop for calculate FDR, sorting P-value before this test
statistics = statistics[statistics.list_count > 1].reset_index(drop=True)

statistics = statistics.sort_values(by ='P',ascending=True).reset_index(drop=True)

statistics['Rank'] = statistics.index + 1


#
#sys.argv[2] = 0.1
#
FDR_val=[]
for x in statistics.Rank:
    FDR_val.append((x/statistics.count()[0])*float(sys.argv[2]))
statistics['FDR'] = FDR_val

## Loop to add boolean value to statistically significant
# T = True if P < FDR
# F = False if P > FDR

significant = []
for index, row in statistics.iterrows():
    if row.P <= row.FDR:
        significant.append('T')
    if row.P > row.FDR:
        significant.append('F')
statistics['Sig'] = significant

statistics = statistics.merge(description, on = 'base', how = 'left')

ff = []
for i in statistics.base.drop_duplicates():
    df = list_cat[list_cat.base == i]
    ff.append([i, ';'.join(df.Entry.tolist())])
enrichment = DataFrame(ff, columns = ['base','entry'])

statistics = statistics.merge(enrichment, on = 'base', how = 'left')
statistics.to_csv('data/Enrichment_analysis_'+sys.argv[1].split('.')[0]+'.tsv', index=None,sep='\t')
