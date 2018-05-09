
# coding: utf-8

# In[ ]:


print('=============== Loading data ===============\n..........')
## Importamos los paquetes requeridos
from pandas import Series, DataFrame 
import pandas as pd
pd.set_option('max_rows',100000)
pd.set_option('max_colwidth',100000)
import urllib.request
import webbrowser
import requests
from bs4 import BeautifulSoup
import re
import shutil, os
import numpy as np
from subprocess import Popen, PIPE, STDOUT
from subprocess import call
import shlex, subprocess
import subprocess


# In[ ]:


keyword=input("\nKeyword search: ")
url_keyword="https://pfam.xfam.org/search/keyword?query="+keyword+""
urllib.request.urlretrieve(url_keyword,'file_download_from_pfam')
###
#!runipy -q ./software/convert_html_to_txt_Pfam_keyword.ipynb
####
aa="sed 's/<\/p><\/td>/#####/g' file_download_from_pfam > 001"
subprocess.call(aa, shell=True)
bb="perl -p -e 's/\n//' 001 > 002"
subprocess.call(bb, shell=True)
cc="perl -p -e 's/#####/\n#####/g' 002 | sed 's/        <[^>]*>//g; s/<[^>]*>/!/g; s/#####.*[0-9]!!PF/PF/g' | grep '^PF.*' | sed 's/!!!/\t/g; s/!!/\t/g; s/!.*//g; /^$/d' > download_correctec"
subprocess.call(cc, shell=True)
####
PFAM_DOMAINS_keyword=pd.read_table('download_correctec',names=['Accession','ID','Description'])
PFAM_DOMAINS_keyword.to_csv('PFAM_DOMAINS_'+keyword,index=None)
###
if os.path.exists("001"): os.remove("001")
if os.path.exists("002"): os.remove("002")   
if os.path.exists("file_download_from_pfam"): os.remove("file_download_from_pfam")
if os.path.exists("download_correctec"): os.remove("download_correctec")
###
print('\n',PFAM_DOMAINS_keyword[0:5])
print('\n========= Filename created: PFAM_DOMAINS_'+keyword+'\n')
count_domains=DataFrame(PFAM_DOMAINS_keyword.count()).rename(columns={0:'Found domains'},inplace=False)
print(count_domains)
print('..........')

