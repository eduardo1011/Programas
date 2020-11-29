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


ppparametros = open('NeVOmics_params.txt', 'r')
parametros = ppparametros.read()
ppparametros.close()

han = open('../NeVOmics_img/KEGG_Organisms.txt', 'r')
dict_org = {}
for line in han:
    line = line.rstrip()
    if re.search('^#', line):
        infokegg = line
        pass
    else:
        separados = line.split('\t')
        dict_org[separados[0]] = [separados[1], separados[2]]
han.close()

# parametros elegidos
method_P = 'FDR' # default

# definimos la localizacion y nombre del archivo elegido
file_path = re.search('filelocation.*', parametros).group().split('=')[1]

# fdr elegido
FDR = float(re.search('keggfdr.*', parametros).group().split('=')[1]) / 100

# T number identificado a partir del organismo
t_number = re.search('keggorganism.*', parametros).group().split('=')[1]

# frefijo identificado a partir del organismo
pref = dict_org[t_number][0]

# organismo seleccionado
organism = dict_org[t_number][1]

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


"""
extrae el numero taxonomico y el nombre del organismo desde uniprot
para compararlo con el ingresado por el usuario,
si no son iguales el proceso se detiene
"""
id_organism = requests.get("https://www.uniprot.org/uniprot/?query="+list_input.Entry[0]+"&sort=score&columns=organism-id,organism&format=tab&limit=1").content.decode()
Prefix = id_organism.split('\t', 2)[1].split('\n')[1]


k = requests.get("https://www.kegg.jp/dbget-bin/www_bget?gn:"+t_number)
k = k.text.rstrip()
tax = re.findall('TAX:.*', k)[0]
tax2 = re.findall('Info&id=\d+..\d+', tax)[0]
Prefix_user = tax2.split('">')[1]

if Prefix_user == Prefix:
    pass
else:
    from tkinter import messagebox
    import tkinter as tk
    import ctypes

    ctypes.windll.shcore.SetProcessDpiAwareness(1)
    root = tk.Tk()
    root.withdraw()
    messagebox.showwarning('Status',
                        'Your selected organism ('+organism+') does not correspond to the \
 organism identified with the UniProt identifiers ('+id_organism.rstrip().split('\t')[-1]+').\n\n\
 !!!Choose the organism correctly!!!')

## Create a folder
os.makedirs('data',exist_ok=True)


# all kegg-id and pathway-description
ee=requests.get('http://rest.kegg.jp/list/pathway/'+pref+'').content.decode()
ee = re.sub('path:|- '+organism[0:5]+'.*','',ee)
kegg_pathways = DataFrame([i.split('\t') for i in ee.split('\n')], columns = ['Path','Term'])
kegg_pathways.to_csv('data/Pathways.txt',sep='\t',index=None)

print(re.sub('#', '', infokegg))

######## extraccion de informacion
if Prefix == '9606': # human
    inf1=requests.get('https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:'+Prefix+'&format=tab&columns=id,database(DisGeNET)').content.decode().rstrip()
    guardar = []
    for i in inf1.split('\n')[1:]:
        if i.split('\t')[1] == '':
            uno = np.nan
        else:
            uno = i.split('\t')[1].split(';')[0]
        guardar.append([i.split('\t')[0], uno])
    inf1 = DataFrame(guardar, columns = ['Entry', 'Entry_Kegg'])
    inf2=requests.get('https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:'+Prefix+'&format=tab&columns=id,database(GeneID)').content.decode().rstrip()
    guardar2 = []
    for i in inf2.split('\n')[1:]:
        if i.split('\t')[1] == '':
            uno = np.nan
        else:
            uno = i.split('\t')[1].split(';')[0]
        guardar2.append([i.split('\t')[0], uno])
    inf2 = DataFrame(guardar2, columns = ['Entry', 'Entry_Kegg'])
    Kegg_Uniprot=pd.concat([inf1,inf2], axis=0,sort=False).dropna().reset_index(drop=True).drop_duplicates()
    dd=requests.get('http://rest.kegg.jp/link/pathway/'+pref+'').content.decode().rstrip()
    guardar3 = [(re.sub('^'+pref+':', '', i.split('\t')[0]), re.sub('path:', '', i.split('\t')[1])) for i in dd.split('\n')]
    kegg_path_ID = DataFrame(guardar3, columns = ['Entry_Kegg','Path'])
else:
    # all kegg-id and pathway-id
    dd=requests.get('http://rest.kegg.jp/link/pathway/'+pref+'').content.decode().rstrip()
    guardar3 = [(re.sub('^'+pref+':', '', i.split('\t')[0]), re.sub('path:', '', i.split('\t')[1])) for i in dd.split('\n')]
    kegg_path_ID = DataFrame(guardar3, columns = ['Entry_Kegg','Path'])
    # all kegg-id and pathway-description
    # all kegg-id and uniprot
    ff=requests.get('http://rest.kegg.jp/conv/uniprot/'+pref+'').content.decode().rstrip()
    guardar4 = [(re.sub('^'+pref+':', '', i.split('\t')[0]), re.sub('up:', '', i.split('\t')[1])) for i in ff.split('\n')]
    Kegg_Uniprot = DataFrame(guardar4, columns = ['Entry_Kegg','Entry'])

allanotacion = Kegg_Uniprot.merge(kegg_path_ID, on = 'Entry_Kegg', how = 'left').dropna().reset_index(drop = True)

 ## definimos el tipo de etiqueta para los nodos
labelnode = 'Gene Name'

if len(inp_file.columns) == 3:
    print('Data with background column\n')
    provicional = list_input[['Background']].rename(columns={'Background':'Entry'})
    background_info = provicional.merge(allanotacion, on = 'Entry', how = 'inner')
    
    # guardar archivo background, con la columna de genes
    background_info[['Entry_Kegg']].drop_duplicates().to_csv('data/Background.txt',index=None)

    # 2.- Preparation of list with pathways
    ## Protein list mapping against "background_info" and then save this list
    list_input_match = list_input[['Entry']].merge(background_info,how="left", on='Entry').dropna().drop_duplicates()
    list_input_match[['Entry_Kegg']].drop_duplicates().to_csv('data/List.txt',index=None)
    
  
    # 3.- background with: Entry	GO, for association file
    background_info[['Entry_Kegg','Path']].to_csv('data/Association.txt',index=None,sep='\t')
else:
    print('Data without background column\n')
    background_info = allanotacion
    
    # guardar archivo background, con la columna de genes
    background_info[['Entry_Kegg']].drop_duplicates().to_csv('data/Background.txt',index=None)

    # 2.- Preparation of list with pathways
    ## Protein list mapping against "background_info" and then save this list
    list_input_match = list_input[['Entry']].merge(background_info,how="left", on='Entry').dropna().drop_duplicates()
    list_input_match[['Entry_Kegg']].drop_duplicates().to_csv('data/List.txt',index=None)
    
  
    # 3.- background with: Entry	GO, for association file
    background_info[['Entry_Kegg','Path']].to_csv('data/Association.txt',index=None,sep='\t')

# descarga el modulo para la estadistica
hd = urllib.request.urlretrieve('https://raw.githubusercontent.com/bioinfproject/bioinfo/master/Folder/HD.py', './HD.py')


analysis = 'Pathways.txt'
subprocess.call(["python", "HD.py", analysis,str(FDR)])

# abrimos los resultados
enrich_P = pd.read_csv('data/Enrichment_analysis_'+analysis.split('.')[0]+'.tsv',sep='\t')

# #En este paso si no hay mas de un termino terminar el proceso porque no se pueden realizar redes con un nodo/proteinas<font>

if enrich_P[enrich_P.Sig == 'T']['FDR'].count() >= 1: # al menos un valor de FDR es significativo      
    results_process_P = enrich_P[enrich_P.Sig == 'T']
    vertice = [(k, j) for i, k in zip(results_process_P.entry.tolist(), results_process_P.base.tolist()) for j in i.split(';')]
    proteins_count_P = len(set([i[1] for i in vertice]))
    GO_count_P = results_process_P.base.count()
else:
    # sin informacion en el kegg
    results_process_P = enrich_P
    no_anotadas = []
    for i in list_input.Entry.drop_duplicates().dropna().tolist():
        if i in list_input_match.Entry.drop_duplicates().tolist():
            continue
        else:
            no_anotadas.append(i)
    
    if len(results_process_P) == 0:
        singleton = 0
    else:
        singleton = int(round(float(results_process_P.Bonf_corr.iloc[-1:]) / float(results_process_P.P.iloc[-1:]), 0))
    
    
    reporte = {'base':[np.nan,
                       'KEGG DB Last-Modified',
                       'Input file name',
                       'Association file name',
                       'Total number of background',
                       'Total number of list',
                       'Background with Pathways',
                       'List input with Pathways',
                       'Non-singletons value for Bonf_corr',
                       'Correction Method',
                       'Value',
                       np.nan,
                       'Proteins with no information in KEGG Pathways',
                       ';'.join(no_anotadas)],
               'list_count':[np.nan,
                             infokegg,
                             file_path, analysis,
                             background_info['Entry'].drop_duplicates().count(),
                             list_input['Entry'].drop_duplicates().count(),
                             background_info['Entry'].drop_duplicates().count(),
                             list_input_match['Entry'].drop_duplicates().count(),
                             int(round(float(results_process_P.Bonf_corr.iloc[-1:]) / float(results_process_P.P.iloc[-1:]), 0)),
                             'FDR',
                             str(FDR)+' ('+str(FDR * 100)+'%)',
                             np.nan, 
                             len(no_anotadas),
                             np.nan]}
    information = DataFrame(reporte)
    informe_final = pd.concat([results_process_P, information], axis=0, sort=False).rename(columns={'base':'Path'})
    informe_final = informe_final[['Path', 'list_count', 'back_count', 'tot_list', 'tot_back', 'P', 'Bonf_corr',
           'Rank', 'FDR', 'Sig', 'Term', 'entry']]
    writer = pd.ExcelWriter('Enrichment_Pathways_Analysis_FDR_'+str(FDR)+'.xlsx')

    informe_final.to_excel(writer,'Significant KEGG Pathways',index=False)

    enrich_P.to_excel(writer,'Enrichment Results',index=False)
    
    writer.save()
    

etiquetas = []
for i in results_process_P.Term:
    i = i.rstrip()
    if len(i.split(' ')) == 1:
        etiquetas.append(i)
    if len(i.split(' ')) == 2:
        etiquetas.append(re.sub(' ', '\n', i))
    if len(i.split(' ')) == 3:
        etiquetas.append(re.sub(' ', '\n', i))
    if len(i.split(' ')) == 4:
        etiquetas.append(' '.join(i.split(' ')[0:2])+'\n'+' '.join(i.split(' ')[2:4]))
    if len(i.split(' ')) > 4:
        etiquetas.append(' '.join(i.split()[0:2])+'\n'+' '.join(i.split()[2:4])+'...')
results_process_P['Short_Term'] = etiquetas



significativos = []
for x in results_process_P.base.drop_duplicates():
    dff = results_process_P[results_process_P.base == x]
    for index, row in dff.iterrows():
        for i in row.entry.split(';'):
            significativos.append([x, row.P, row.FDR, row.Term, row.Short_Term, i])


keggtabla = DataFrame(significativos, columns = ['Path', 'P', 'FDR', 'Term', 'Short_Term', 'Entry_Kegg'])
keggtabla['LogminFDR'] = -np.log10(keggtabla.FDR)
keggtabla['LogminP'] = -np.log10(keggtabla.P)
n = 0
ranked = []
for i in keggtabla['Entry_Kegg'].drop_duplicates():
    n+=1
    ranked.append([i, str(n)])
rank = DataFrame(ranked, columns = ['Entry_Kegg', 'label'])


keggtabla = keggtabla.merge(rank, on = 'Entry_Kegg', how = 'left')
keggtabla = keggtabla.merge(list_input_match, on = ['Entry_Kegg', 'Path'], how = 'left')

keggtabla = keggtabla.merge(list_input[['Entry', 'values']], on = 'Entry', how = 'left')

edges_frame_excel = keggtabla[['Path','Entry_Kegg','Entry','Term','values']]


if labelnode == 'Gene Name':
    pass
if labelnode == 'UniProt ID':
    keggtabla = keggtabla.rename({'Entry_Kegg':'Entry', 'Entry':'Entry_Kegg'}, axis='columns')

no_anotadas = []
for i in list_input.Entry.drop_duplicates().dropna().tolist():
    if i in list_input_match.Entry.drop_duplicates().tolist():
        continue
    else:
        no_anotadas.append(i)


reporte = {'base':[np.nan,
                   'KEGG DB Last-Modified',
                   'Input file name',
                   'Association file name',
                   'Total number of background',
                   'Total number of list',
                   'Background with Pathways',
                   'List input with Pathways',
                   'Non-singletons value for Bonf_corr',
                   'Correction Method',
                   'Value',
                   np.nan,
                   'Proteins with no information in KEGG Pathways',
                   ';'.join(no_anotadas)],
        'list_count':[np.nan,
                      infokegg,
                      file_path, analysis,
                      background_info['Entry'].drop_duplicates().count(),
                      list_input['Entry'].drop_duplicates().count(),
                      background_info['Entry'].drop_duplicates().count(),
                      list_input_match['Entry'].drop_duplicates().count(),
                      int(round(float(results_process_P.Bonf_corr.iloc[-1:]) / float(results_process_P.P.iloc[-1:]), 0)),
                      'FDR',
                      str(FDR)+' ('+str(FDR * 100)+'%)',
                      np.nan,
                      len(no_anotadas),
                      np.nan]}
information = DataFrame(reporte)
informe_final = pd.concat([results_process_P, information], axis=0, sort=False).rename(columns={'base':'Path'})


writer = pd.ExcelWriter('Enrichment_Pathways_Analysis_FDR_'+str(FDR)+'.xlsx')

informe_final.to_excel(writer,'Significant KEGG Pathways',index=False)

enrich_P.to_excel(writer,'Enrichment Results',index=False)

edges_frame_excel.to_excel(writer,'Edges Pathways',index=False)
writer.save()

