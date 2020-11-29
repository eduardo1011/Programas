import re
import pandas as pd
import requests
import os
from pandas import DataFrame
import urllib.request
import subprocess
import numpy as np
import xlsxwriter
import warnings
warnings.filterwarnings("ignore")
import sys


ppparametros = open('NeVOmics_params.txt', 'r')
parametros = ppparametros.read()
ppparametros.close()

# parametros elegidos
method_P = 'FDR' # default
# definimos la localizacion y nombre del archivo elegido
file_path = re.search('filelocation.*', parametros).group().split('=')[1]
# Biological process
bpfdr = float(re.search('bpfdr.*', parametros).group().split('=')[1])  / 100
# Molucular Function
mffdr = float(re.search('mffdr.*', parametros).group().split('=')[1])  / 100
# Cellular Component
ccfdr = float(re.search('ccfdr.*', parametros).group().split('=')[1])  / 100

anotacion_uniprot = '1' # default

anotacion_goa = re.search('anotacion_goa.*', parametros).group().split('=')[1]

labelnode = 'Gene Name'

## read file
inp_file=pd.read_csv(file_path,sep='\t',header=None)   
    
## explore input file
if len(inp_file.columns) == 1:
    hayvalores = 'nohayvalores'
    inp_file['values'] = 1
    ## only gene list
    list_input=inp_file.rename(columns={0:'Entry'},index=str) 
if len(inp_file.columns) == 2:
    hayvalores = 'sihay'
    ## gene list and vales
    list_input=inp_file.rename(columns={0:'Entry',1:'values'},index=str) ## gene list and values
if len(inp_file.columns) == 3:
    hayvalores = 'sihay'
    ## gene list, values and background
    list_input=inp_file.rename(columns={0:'Entry',1:'values',2:'Background'},index=str)

## exttract id-organism
id_organism = requests.get("https://www.uniprot.org/uniprot/?query="+list_input.Entry[0]+"&sort=score&columns=organism-id&format=tab&limit=1").content.decode()
Prefix = id_organism.split('\n')[1]

if id_organism == '':
    print('Organism not found')
    #del_stop_process()
else:
    pass

os.makedirs('data',exist_ok=True)

# save all ontology from go-basic.obo file from GOC

import os, fnmatch
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

gobasic = open('../NeVOmics_img/go-basic.obo', 'r')
for line in gobasic:
    if re.search('data-version: .*', line):
        pat = re.search('data-version: .*', line).group()
        go_version = re.sub('data-version: releases.', '', pat)
        print('\nOntology version: ', go_version)
        print('Downloaded from: http://geneontology.org/docs/download-ontology/', '\n')
        break
        
with open('../NeVOmics_img/go-basic.obo', 'r') as g:
    go_obo = g.read()
g.close()
ontology_file = go_obo.split('[Term]')

aspect = {'biological_process':'P', 'molecular_function':'F', 'cellular_component':'C'}
items = []
for i in ontology_file[1:len(ontology_file)]:
    items.append([i.split('\n')[1].split(': ')[1],
                 i.split('\n')[2].split(': ')[1],
                 aspect[i.split('\n')[3].split(': ')[1]]])
ontologia = DataFrame(items, columns = ['GO', 'Term', 'Aspect'])


ontologia[ontologia['Aspect'].str.contains('P') == True][['GO','Term']].to_csv('data/GO_BP.txt', sep = '\t', index = None)
ontologia[ontologia['Aspect'].str.contains('F') == True][['GO','Term']].to_csv('data/GO_MF.txt', sep = '\t', index = None)
ontologia[ontologia['Aspect'].str.contains('C') == True][['GO','Term']].to_csv('data/GO_CC.txt', sep = '\t', index = None)

##############################################################
################           Uniprot         ###################
##############################################################

file_uniprot = find('annotation_'+Prefix, '../')
if file_uniprot == ('' or []):
    uni = urllib.request.urlretrieve('https://www.uniprot.org/uniprot/?query=organism:'+Prefix+'&format=tab&columns=id,genes,go-id', 'annotation_'+Prefix)
    prot_version = uni[1]['Last-Modified']
    go_uniptot_version = uni[1]['Last-Modified']
    print('UniProtKB version: ', prot_version)
    print('Entries: ', uni[1]['X-Total-Results'])
    with open(uni[0], 'a') as fq:
        fq.write('#'+prot_version)
        fq.close()
    
    acc_GOid=pd.read_csv('annotation_'+Prefix,sep='\t')#.dropna().reset_index(drop=True)
    acc_GOid.columns = ['Entry', 'Gene', 'GO']
else:
    file_uniprot = re.sub('\\\\', '/', file_uniprot[0])
    print('It already exists:', file_uniprot)
    acc_GOid=pd.read_csv(file_uniprot,sep='\t')#.dropna().reset_index(drop=True)
    acc_GOid.columns = ['Entry', 'Gene', 'GO']
    go_uniptot_version = re.sub('#', '', acc_GOid.Entry.tolist()[-1])
    print('UniProtKB version: ', re.sub('#', '', acc_GOid.Entry.tolist()[-1]))
    acc_GOid = acc_GOid[acc_GOid.Entry.str.contains('#') == False]
    print('Entries: ', acc_GOid.Entry.count())

# ## exploracion de la anotacion de Uniprot

# es un df con Entries sin anotación GO en Uniprot
con_nas = acc_GOid[pd.isna(acc_GOid['GO'])]

# es un df con Entries que tienen anotacion en Uniprot, pero algunos genes no tienen identificador
sin_nas = acc_GOid[pd.notna(acc_GOid['GO'])]
# del df anterior extraigo un df de genes que no tienen nombre, en otra celda edito estos genes
genes_sin_name = sin_nas[pd.isna(sin_nas['Gene'])]
sin_nas = sin_nas[pd.notna(sin_nas['Gene'])]

# asigno el Entry como nombre del gen
new_gene_name = []
for i, j in genes_sin_name.iterrows():
    new_gene_name.append([j.Entry, 'Entry:'+j.Entry, j.GO])

df_new_gene_name = DataFrame(new_gene_name, columns = ['Entry', 'Gene', 'GO'])


Entry_GOid = pd.concat([sin_nas, df_new_gene_name]).drop_duplicates()

uniprot_anotation_org = []
for index, row in Entry_GOid.iterrows():
    for j in row.GO.split('; '):
        if re.search('GO:\d+', j):
            uniprot_anotation_org.append([row.Entry, row.Gene, re.search('GO:\d+', j).group()])
# este df contiene toda la anotación GO en uniprot del organismo en estudio
Entry_GOid_annotated =DataFrame(uniprot_anotation_org, columns = ['Entry', 'Gene', 'GO'])

# ## df con toda la informacion funcional GOA capturada usando los ids de uniprot
uniprot_entry_go_term = Entry_GOid_annotated.merge(ontologia, on = 'GO', how = 'left').dropna()
uniprot_entry_go_term['Gene'] = [i.split(' ')[0] for i in uniprot_entry_go_term.Gene.tolist()]
uniprot_entry_go_term = uniprot_entry_go_term[['Entry', 'GO', 'Term', 'Aspect', 'Gene']]

total = len(uniprot_entry_go_term.Entry.drop_duplicates().tolist())
print('\nUniProt Entries with GO Terms:', total)


def hojas(dict_hoja = dict()):
    final = []
    for i in dict_hoja.keys():
        for j in dict_hoja[i]['results']:
            final.append([re.sub('UniProtKB:','', j['geneProductId']),
                           j['qualifier'],
                           j['goId'],
                           j['goName'],
                           j['goEvidence'],
                           j['goAspect'],
                           j['evidenceCode'],
                           j['reference'],
                           j['withFrom'],
                           j['taxonId'],
                           j['taxonName'],
                           j['assignedBy'],
                           j['extensions'],
                           j['targetSets'],
                           j['symbol'],
                           j['date'],
                           j['synonyms'],
                           j['name']])
    return final


if anotacion_goa == '1':

    file_goa = ''.join(find('Complete_Annotation_'+Prefix+'_goa', '../'))
    file_goa1 = re.sub('\\\\', '/', file_goa)
    file_goa2 = file_goa1.split('/')[-1]
    
    if file_goa2 == '':
        from bioservices import QuickGO
        qg = QuickGO()
        info_goa = requests.get('https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?geneProductId='+Entry_GOid_annotated.Entry.drop_duplicates().tolist()[0])
        goa_information = info_goa.headers['Date']
        print('GOA annotation:', goa_information)
        print('The full annotation will be downloaded, this may take some time, from 10 min to 2 h')
        # descarga de la anotacion completa
        ##############################################################
        ################       Uniprot GOA         ###################
        ##############################################################
        ## GOA Proteome Sets
        # descarga anotacion completa GOA usando bioservices module
        # los arhcivos .gaf están incompletos y no muestran la anotación completa, por eso se usa bioservices
        #

        from datetime import datetime

        print('\nTotal UniProt entries expected:', total)
        print('Mapping: (UniprotKB Entries | Found Entries | Missing Entries)')
        inicio = datetime.now()
        dl = 0
        resultado = []
        suma = 0
        for k in range(0, total,100):
            tim = datetime.now() - inicio
            serie_cien = ''.join(','.join(Entry_GOid_annotated.Entry.drop_duplicates().tolist()[k:k+100]))
            ##
            exploracion = qg.Annotation(page = 1, geneProductId  = serie_cien)
            total_pages = exploracion['pageInfo']['total']
            lista_paginas = list(range(1,total_pages+1))
            res = {}
            ids = []
            for i in list(range(1,total_pages+1)):
                res[i] = qg.Annotation(page = i, geneProductId  = serie_cien)
                if res[i] == 400:
                    res.pop(i)
                for f in res[i]['results']:
                    ids.append(re.sub('-.*', '',f['geneProductId'].split(':')[1]))
            suma += len(set(ids))
            ##
            dl += len(serie_cien.split(','))
            done = int(40 * dl / total)
            sys.stdout.write("\r"+'{}'.format(tim).split('.')[0]+" [%s%s] (%s | %s | %s)" % ('>' * done, ' ' * (40-done),
                                                                                             dl, suma, (dl-suma)))
            sys.stdout.flush()

    
            header = ['Entry','qualifier','goId','goName','goEvidence','goAspect','evidenceCode','reference',
                        'withFrom','taxonId','taxonName','assignedBy','extensions','targetSets','symbol',
                        'date','synonyms','name']
            resultado.append(DataFrame(hojas(dict_hoja = res), columns = header))
    
        print('\n')
        complete_annotation = pd.concat(resultado)
        complete_annotation.to_csv('Complete_Annotation_'+Prefix+'_goa', sep = '\t',index=None)
        with open('Complete_Annotation_'+Prefix+'_goa', 'a') as fq:
            fq.write('#'+goa_information)
            fq.close()
    else:
        print('\nIt already exists:', file_goa1)
        complete_annotation = pd.read_csv(file_goa1, sep='\t')
        goa_information = re.sub('#', '', complete_annotation.Entry.tolist()[-1])
        print('GOA annotation:', goa_information)
        complete_annotation = complete_annotation[complete_annotation.Entry.str.contains('#') == False]
        
    ### recuperacion de Entries no encontrados en Quickgo, se obtienen a partir de uniprot   
    # removemos Entries que tienen guiones
    complete_annotation = complete_annotation[complete_annotation['Entry'].str.contains('-') == False]
    goa_annotation = complete_annotation[['Entry', 'goId', 'symbol']]
    goa_annotation.columns = ['Entry', 'GO', 'Gene']
    goa_entry_go_term = goa_annotation.merge(ontologia, on = 'GO', how = 'left').dropna()
    # recuperadas no mapeadas = rnm
    rnm = uniprot_entry_go_term.merge(goa_entry_go_term, on = 'Entry', how = 'left')
    rnm = rnm[pd.isnull(rnm).any(axis=1)]
    rnm2 = rnm[['Entry']].drop_duplicates().merge(uniprot_entry_go_term, on = 'Entry', how = 'left')
    ## df con toda la informacion funcional GOA capturada usando los ids de uniprot
    goa_entry_go_term = pd.concat([goa_entry_go_term, rnm2])
    total = len(goa_entry_go_term.Entry.drop_duplicates().tolist())
    print('\nGOA Entries:', total)  

else:
    pass

# descarga el modulo para la estadistica
hd = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/HD.py', './HD.py')

#############  funciones

# funcion para crear archivos para enriquecimiento, ingresar diferentes data frames
def enrichment_files(df = DataFrame([])):
    if len(inp_file.columns) == 3:
        print('With background')
        provicional = list_input[['Background']].rename(columns={'Background':'Entry'})
        background_info = provicional.merge(df, on = 'Entry', how = 'inner')
    
        # guardar archivo background, con la columna de genes
        background_info[['Entry']].drop_duplicates().to_csv('data/Background.txt',index=None)
        print('data/Background.txt')

        # 2.- Preparation of list with pathways
        ## Protein list mapping against "background_info" and then save this list
        list_input_match = list_input[['Entry']].merge(background_info,how="left", on='Entry').dropna().drop_duplicates()
        list_input_match[['Entry']].drop_duplicates().to_csv('data/List.txt',index=None)
        print('data/List.txt')
    
        # 3.- background with: Entry	GO, for association file
        background_info[['Entry','GO']].to_csv('data/Association.txt',index=None,sep='\t')
        print('data/Association.txt')
        
    else:
        print('No background')
        background_info = uniprot_entry_go_term
        
        # guardar archivo background, con la columna de genes
        background_info[['Entry']].drop_duplicates().to_csv('data/Background.txt',index=None)
        print('data/Background.txt')
    
        # 2.- Preparation of list with pathways
        ## Protein list mapping against "background_info" and then save this list
        list_input_match = list_input[['Entry']].merge(background_info,how="left", on='Entry').dropna().drop_duplicates()
        list_input_match[['Entry']].drop_duplicates().to_csv('data/List.txt',index=None)
        print('data/List.txt')
    
        # 3.- background with: Entry	GO, for association file
        background_info[['Entry','GO']].to_csv('data/Association.txt',index=None,sep='\t')
        print('data/Association.txt')
    ######
    no_anotadas = []
    for i in list_input.Entry.drop_duplicates().dropna().tolist():
        if i in list_input_match.Entry.drop_duplicates().tolist():
            continue
        else:
            no_anotadas.append(i)
    return list_input_match, no_anotadas

# funcion para explorar si hay terminos enriquecidos, si hay crea un df,
# si no hay crea los archivos excel y termina el proceso

def filtro_significancia(df = DataFrame([]), info = '', asso_file = '', fdr_val = 0, no_annot = [], db = ''):
    if df[df.Sig == 'T']['FDR'].count() >= 1: # al menos un valor de FDR es significativo      
        results_sig = df[df.Sig == 'T']
        vertice = [(k, j) for i, k in zip(results_sig.entry.tolist(), results_sig.base.tolist()) for j in i.split(';')]
        proteins_count = len(set([i[1] for i in vertice]))
        GO_count = results_sig.base.count()
        
        return results_sig 
        
    else:
        # si no hay GO terms enriquecidos se generan los archivos excel sin los edges
        if len(df) < 1:
            input_background = 0
            go_background = 0
            go_lista = 0
            singletons_value = 0
        else:
            input_background = int(float(df.tot_back.iloc[0:1]))
            go_background = df.tot_back.iloc[0]
            go_lista = df.tot_list.iloc[0]
            singletons_value = int(round(float(df.Bonf_corr.iloc[-1:]) / float(df.P.iloc[-1:]), 0))
        
        results_sig = df[df.Sig == 'T']
        reporte = {'base':[np.nan,
                           'GO DB Last-Modified',
                           'Input file name',
                           'Association file name',
                           'Total number of background',
                           'Total number of list',
                           'Background with GO Terms',
                           'List input with GO Terms',
                           'Non-singletons value for Bonf_corr',
                           'Correction Method',
                           'Value',
                           np.nan,
                           'Proteins with no information in UniProtKB',
                           ';'.join(no_annot)],
                'list_count':[np.nan,
                              info,
                              file_path,
                              asso_file,
                              input_background,
                              list_input['Entry'].drop_duplicates().count(),
                              go_background,
                              go_lista,
                              singletons_value,
                              'FDR',
                              str(fdr_val)+' ('+str(np.round(fdr_val * 100,1))+'%)',
                              np.nan,
                              len(no_annot),
                              np.nan]}
        information = DataFrame(reporte)
        informe_final = pd.concat([results_sig, information], axis=0, sort=False).rename(columns={'base':'GO'})
    
        informe_final = informe_final[['GO', 'list_count', 'back_count', 'tot_list', 'tot_back', 'P', 'Bonf_corr',
           'Rank', 'FDR', 'Sig', 'Term', 'entry']]
        writer = pd.ExcelWriter(db+'_Enrichment_'+asso_file.split('.')[0]+'_FDR_'+str(fdr_val)+'.xlsx')

        informe_final.to_excel(writer,'Significant GO Terms',index=False)

        df.to_excel(writer,'Enrichment Results',index=False)
    
        writer.save()
        
        print('!!!!!!!!!!!!!!\nLess than 2 significant terms were identified in '+asso_file.split('.')[0]+' for the chosen FDR,'+              ' try another FDR.\nTo create networks it is necessary to obtain at least 2 terms.')

        return DataFrame([None])




def termino_corto(df = DataFrame([])):
    etiquetas = []
    for i in df.Term:
        i = i.rstrip()
        if len(i.split(' ')) == 1:
            etiquetas.append(i)
        if len(i.split(' ')) == 2:
            etiquetas.append(re.sub(' ', '\n', i))
        if len(i.split(' ')) == 3:
            etiquetas.append(re.sub(' ', '\n', i))
        if len(i.split(' ')) == 4:
            etiquetas.append(' '.join(i.split()[0:2])+'\n'+' '.join(i.split()[2:4]))
        if len(i.split(' ')) > 4:
            etiquetas.append(' '.join(i.split()[0:2])+'\n'+' '.join(i.split()[2:4])+'...')
    return etiquetas

#########################################################################
categorias = ['GO_BP.txt', 'GO_MF.txt', 'GO_CC.txt']
fdrs = [bpfdr, mffdr, ccfdr]

if anotacion_uniprot == '1':
    print('*****UniProtKB')
    no_anotadas_uniprot = enrichment_files(df = uniprot_entry_go_term)
    uniprot_enrich = {}
    uniprot_signif = {}
    for i, j in zip(categorias, fdrs):
        subprocess.call(["python", "HD.py", i, str(j)])
        enrich = pd.read_csv('data/Enrichment_analysis_'+i.split('.')[0]+'.tsv',sep='\t')
        uniprot_enrich[i.split('.')[0]] = enrich
        significantes = enrich[enrich.Sig == 'T']
        uniprot_signif[i.split('.')[0]] = significantes
        print('Finished (UniProt):', i.split('.txt')[0])
    ###
    # los que no tienen terminos significantes se descartarán
    uni_info = list(np.repeat(go_uniptot_version, len(categorias)))
    uni_noanno = [no_anotadas_uniprot[1], no_anotadas_uniprot[1], no_anotadas_uniprot[1]]
    uni_database = ['UniProt', 'UniProt', 'UniProt']
    aprobados_uniprot = {}
    for i, j, k, l, m in zip(categorias, fdrs, uni_info, uni_noanno, uni_database):
        #print(i.split('.')[0])
        aprobados_uniprot[i.split('.')[0]] = filtro_significancia(df = uniprot_enrich[i.split('.')[0]],
                                      asso_file = i,
                                      fdr_val = j,
                                      info = k,
                                      no_annot = l,
                                      db = m)
else:
    pass
if anotacion_goa == '1':
    print('\n*****UniProt GOA')
    no_anotadas_goa = enrichment_files(df = goa_entry_go_term)
    goa_enrich = {}
    goa_signif = {}
    for i, j in zip(categorias, fdrs):
        subprocess.call(["python", "HD.py", i, str(j)])
        enrich = pd.read_csv('data/Enrichment_analysis_'+i.split('.')[0]+'.tsv',sep='\t')
        goa_enrich[i.split('.')[0]] = enrich
        significantes = enrich[enrich.Sig == 'T']
        goa_signif[i.split('.')[0]] = significantes
        print('Finished (GOA):', i.split('.txt')[0])
    ###
    # los que no tienen terminos significantes se descartarán
    goa_info = list(np.repeat(goa_information, len(categorias)))
    goa_noanno = [no_anotadas_goa[1], no_anotadas_goa[1], no_anotadas_goa[1]]
    goa_database = ['UniProtGOA', 'UniProtGOA', 'UniProtGOA']
    aprobados_goa = {}
    for i, j, k, l, m in zip(categorias, fdrs, goa_info, goa_noanno, goa_database):
        #print(i.split('.')[0])
        aprobados_goa[i.split('.')[0]] = filtro_significancia(df = goa_enrich[i.split('.')[0]],
                                      asso_file = i,
                                      fdr_val = j,
                                      info = k,
                                      no_annot = l,
                                      db = m)
else:
    pass





cats = ['GO_BP', 'GO_MF', 'GO_CC']
###### UniProt ########################

if anotacion_uniprot == '1':
    go_tablas_uniprot = {}
    edges_frame_excel_uniprot = {}
    for z in cats:
        if aprobados_uniprot[z].count().iloc[0] > 1:
            df = aprobados_uniprot[z]
            df['Short_Term'] = termino_corto(df = aprobados_uniprot[z])
            
            significativos = []
            for x in df.base.drop_duplicates():
                dff = df[df.base == x]
                for index, row in dff.iterrows():
                    for i in row.entry.split(';'):
                        significativos.append([x, row.P, row.FDR, row.Term, row.Short_Term, i])
            gotabla = DataFrame(significativos, columns = ['GO', 'P', 'FDR', 'Term', 'Short_Term', 'Entry'])
            gotabla['LogminFDR'] = -np.log10(gotabla.FDR)
            gotabla['LogminP'] = -np.log10(gotabla.P)
            n = 0
            ranked = []
            for i in gotabla['Entry'].drop_duplicates():
                n+=1
                ranked.append([i, str(n)])
            rank = DataFrame(ranked, columns = ['Entry', 'label'])
            
            gotabla = gotabla.merge(rank, on = 'Entry', how = 'left')
            gotabla = gotabla.merge(no_anotadas_uniprot[0][['Entry', 'GO']], on = ['Entry', 'GO'], how = 'left')
            gotabla = gotabla.merge(uniprot_entry_go_term[['Entry', 'Gene']], on = 'Entry', how = 'left')
            gotabla = gotabla.merge(list_input[['Entry', 'values']], on = 'Entry', how = 'left')

            edges_frame_excel = gotabla[['GO','Entry', 'Gene', 'Term','values']]
            edges_frame_excel_uniprot[z] = edges_frame_excel
            if labelnode == 'Gene Name':
                gotabla = gotabla.rename({'Gene':'Entry', 'Entry':'Gene'}, axis='columns')
            if labelnode == 'UniProt ID':
                pass
            go_tablas_uniprot[z] = gotabla.drop_duplicates().reset_index(drop = True)
            del gotabla
            del edges_frame_excel
        else:
            if aprobados_uniprot[z].count().iloc[0] == 1:
                df = aprobados_uniprot[z]
                df['Short_Term'] = termino_corto(df = aprobados_uniprot[z])

                significativos = []
                for x in df.base.drop_duplicates():
                    dff = df[df.base == x]
                    for index, row in dff.iterrows():
                        for i in row.entry.split(';'):
                            significativos.append([x, row.P, row.FDR, row.Term, row.Short_Term, i])
                gotabla = DataFrame(significativos, columns = ['GO', 'P', 'FDR', 'Term', 'Short_Term', 'Entry'])
                gotabla['LogminFDR'] = -np.log10(gotabla.FDR)
                gotabla['LogminP'] = -np.log10(gotabla.P)
                n = 0
                ranked = []
                for i in gotabla['Entry'].drop_duplicates():
                    n+=1
                    ranked.append([i, str(n)])
                rank = DataFrame(ranked, columns = ['Entry', 'label'])

                gotabla = gotabla.merge(rank, on = 'Entry', how = 'left')
                gotabla = gotabla.merge(no_anotadas_uniprot[0][['Entry', 'GO']], on = ['Entry', 'GO'], how = 'left')
                gotabla = gotabla.merge(uniprot_entry_go_term[['Entry', 'Gene']], on = 'Entry', how = 'left')
                gotabla = gotabla.merge(list_input[['Entry', 'values']], on = 'Entry', how = 'left')

                edges_frame_excel = gotabla[['GO','Entry', 'Gene', 'Term','values']]
                edges_frame_excel_uniprot[z] = edges_frame_excel
                if labelnode == 'Gene Name':
                    gotabla = gotabla.rename({'Gene':'Entry', 'Entry':'Gene'}, axis='columns')
                if labelnode == 'UniProt ID':
                    pass
                go_tablas_uniprot[z] = gotabla.drop_duplicates().reset_index(drop = True)
                del gotabla
                del edges_frame_excel
        

###### GOA ########################
if anotacion_goa == '1':
    go_tablas_goa = {}
    edges_frame_excel_goa = {}
    for z in cats:
        if aprobados_goa[z].count().iloc[0] > 1:
            df = aprobados_goa[z]
            df['Short_Term'] = termino_corto(df = aprobados_goa[z])
            
            significativos = []
            for x in df.base.drop_duplicates():
                dff = df[df.base == x]
                for index, row in dff.iterrows():
                    for i in row.entry.split(';'):
                        significativos.append([x, row.P, row.FDR, row.Term, row.Short_Term, i])
            gotabla = DataFrame(significativos, columns = ['GO', 'P', 'FDR', 'Term', 'Short_Term', 'Entry'])
            gotabla['LogminFDR'] = -np.log10(gotabla.FDR)
            gotabla['LogminP'] = -np.log10(gotabla.P)
            n = 0
            ranked = []
            for i in gotabla['Entry'].drop_duplicates():
                n+=1
                ranked.append([i, str(n)])
            rank = DataFrame(ranked, columns = ['Entry', 'label'])
            
            gotabla = gotabla.merge(rank, on = 'Entry', how = 'left')
            gotabla = gotabla.merge(no_anotadas_goa[0][['Entry', 'GO']], on = ['Entry', 'GO'], how = 'left')
            gotabla = gotabla.merge(goa_entry_go_term[['Entry', 'Gene']], on = 'Entry', how = 'left')
            gotabla = gotabla.merge(list_input[['Entry', 'values']], on = 'Entry', how = 'left')

            edges_frame_excel = gotabla[['GO','Entry', 'Gene', 'Term','values']]
            edges_frame_excel_goa[z] = edges_frame_excel
            if labelnode == 'Gene Name':
                gotabla = gotabla.rename({'Gene':'Entry', 'Entry':'Gene'}, axis='columns')
            if labelnode == 'UniProt ID':
                pass
            go_tablas_goa[z] = gotabla.drop_duplicates().reset_index(drop = True)
        else:
            if aprobados_goa[z].count().iloc[0] == 1:
                df = aprobados_goa[z]
                df['Short_Term'] = termino_corto(df = aprobados_goa[z])

                significativos = []
                for x in df.base.drop_duplicates():
                    dff = df[df.base == x]
                    for index, row in dff.iterrows():
                        for i in row.entry.split(';'):
                            significativos.append([x, row.P, row.FDR, row.Term, row.Short_Term, i])
                gotabla = DataFrame(significativos, columns = ['GO', 'P', 'FDR', 'Term', 'Short_Term', 'Entry'])
                gotabla['LogminFDR'] = -np.log10(gotabla.FDR)
                gotabla['LogminP'] = -np.log10(gotabla.P)
                n = 0
                ranked = []
                for i in gotabla['Entry'].drop_duplicates():
                    n+=1
                    ranked.append([i, str(n)])
                rank = DataFrame(ranked, columns = ['Entry', 'label'])

                gotabla = gotabla.merge(rank, on = 'Entry', how = 'left')
                gotabla = gotabla.merge(no_anotadas_goa[0][['Entry', 'GO']], on = ['Entry', 'GO'], how = 'left')
                gotabla = gotabla.merge(goa_entry_go_term[['Entry', 'Gene']], on = 'Entry', how = 'left')
                gotabla = gotabla.merge(list_input[['Entry', 'values']], on = 'Entry', how = 'left')

                edges_frame_excel = gotabla[['GO','Entry', 'Gene', 'Term','values']]
                edges_frame_excel_goa[z] = edges_frame_excel
                if labelnode == 'Gene Name':
                    gotabla = gotabla.rename({'Gene':'Entry', 'Entry':'Gene'}, axis='columns')
                if labelnode == 'UniProt ID':
                    pass
                go_tablas_goa[z] = gotabla.drop_duplicates().reset_index(drop = True)

def crear_excel(df = DataFrame([]), df_edges = DataFrame([]), info = '',
                asso_file = '', fdr_val = 0, no_annot = [], db = ''):
    results_sig = df[df.Sig == 'T']
    #
    gos = results_sig.base.drop_duplicates().tolist()
    lista_sup = re.sub(':', '%3A', '%2C'.join(gos))
    quickgo_url = 'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/%7Bids%7D/chart?ids='
    quickgo = quickgo_url+lista_sup
    #
    reporte = {'base':[np.nan,
                           'GO DB Last-Modified',
                           'Input file name',
                           'Association file name',
                           'Total number of background',
                           'Total number of list',
                           'Background with GO Terms',
                           'List input with GO Terms',
                           'Non-singletons value for Bonf_corr',
                           'Correction Method',
                           'Value',
                           np.nan,
                           quickgo,
                           np.nan,
                           'Proteins with no information in UniProtKB',
                           ';'.join(no_annot)],
                'list_count':[np.nan,
                              info,
                              file_path,
                              asso_file,
                              int(float(df.tot_back.iloc[0:1])),
                              list_input['Entry'].drop_duplicates().count(),
                              df.tot_back.iloc[0],
                              df.tot_list.iloc[0],
                              int(round(float(df.Bonf_corr.iloc[-1:]) / float(df.P.iloc[-1:]), 0)),
                              'FDR',
                              str(fdr_val)+' ('+str(np.round(fdr_val * 100,1))+'%)',
                              np.nan,
                              'Copy and paste the url into your browser',
                              np.nan,
                              len(no_annot),
                              np.nan]}
    information = DataFrame(reporte)
    informe_final = pd.concat([results_sig, information], axis=0, sort=False).rename(columns={'base':'GO'})

    writer = pd.ExcelWriter(db+'_Enrichment_'+asso_file.split('.')[0]+'_FDR_'+str(fdr_val)+'.xlsx',
                           engine='xlsxwriter',options={'strings_to_urls': False})

    informe_final.to_excel(writer,'Significant GO Terms',index=False)

    df.to_excel(writer,'Enrichment Results',index=False)
    
    df_edges.drop_duplicates().to_excel(writer,'Edges Pathways',index=False)
    
    writer.save()

################ 
## UniProt
################
if anotacion_uniprot == '1':
    if ('GO_BP' in list(go_tablas_uniprot.keys())) == True:
        crear_excel(df = uniprot_enrich['GO_BP'],
                    df_edges = edges_frame_excel_uniprot['GO_BP'],
                    info = uni_info[0],
                    asso_file = categorias[0],
                    fdr_val = fdrs[0],
                    no_annot = uni_noanno[0],
                    db = uni_database[0])
    else:
        print('There are no enriched terms for BP in UniProtKB')
        pass
    if ('GO_MF' in list(go_tablas_uniprot.keys())) == True:
        crear_excel(df = uniprot_enrich['GO_MF'],
                    df_edges = edges_frame_excel_uniprot['GO_MF'],
                    info = uni_info[1],
                    asso_file = categorias[1],
                    fdr_val = fdrs[1],
                    no_annot = uni_noanno[1],
                    db = uni_database[1])
    else:
        print('There are no enriched terms for MF in UniProtKB')
        pass
    if ('GO_CC' in list(go_tablas_uniprot.keys())) == True:
        crear_excel(df = uniprot_enrich['GO_CC'],
                    df_edges = edges_frame_excel_uniprot['GO_CC'],
                    info = uni_info[2],
                    asso_file = categorias[2],
                    fdr_val = fdrs[2],
                    no_annot = uni_noanno[2],
                    db = uni_database[2])
    else:
        print('There are no enriched terms for CC in UniProtKB')
        pass
################ 
## GOA
################
if anotacion_goa == '1':
    if ('GO_BP' in list(go_tablas_goa.keys())) == True:
        crear_excel(df = goa_enrich['GO_BP'],
                    df_edges = edges_frame_excel_goa['GO_BP'],
                    info = goa_info[0],
                    asso_file = categorias[0],
                    fdr_val = fdrs[0],
                    no_annot = goa_noanno[0],
                    db = goa_database[0])
    else:
        print('There are no enriched terms for BP in GOA')
        pass
    if ('GO_MF' in list(go_tablas_goa.keys())) == True:
        crear_excel(df = goa_enrich['GO_MF'],
                    df_edges = edges_frame_excel_goa['GO_MF'],
                    info = goa_info[1],
                    asso_file = categorias[1],
                    fdr_val = fdrs[1],
                    no_annot = goa_noanno[1],
                    db = goa_database[1])
    else:
        print('There are no enriched terms for MF in GOA')
        pass
    if ('GO_CC' in list(go_tablas_goa.keys())) == True:
        crear_excel(df = goa_enrich['GO_CC'],
                    df_edges = edges_frame_excel_goa['GO_CC'],
                    info = goa_info[2],
                    asso_file = categorias[2],
                    fdr_val = fdrs[2],
                    no_annot = goa_noanno[2],
                    db = goa_database[2])
    else:
        print('There are no enriched terms for CC in GOA')
        pass
