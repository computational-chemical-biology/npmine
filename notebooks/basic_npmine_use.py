#!/usr/bin/env python
# coding: utf-8

# In[1]:


# NOTE: The way the functions are arranged 
# is a bit repetitive, you can group method in modules
# then you can from npmine.util import something
from npmine.retrieve_doi import retrieve_doi
from npmine.retrieve_chemical_entities import retrieve_chemical_entities
from npmine.retrieve_chemical_entities import download_pdf
from npmine.retrieve_chemical_entities_from_image import retrieve_chemical_entities_from_image

import json


# In[2]:


# Retrieve all DOis from this journal
# NOTE: The function should contain the parameter referring to the journal
doi = retrieve_doi()


# In[3]:


doi


# In[11]:


# write doi file
#with open('../data/doi_dnp.json', 'w+') as outfile:
#    json.dump(doi, outfile, indent=4)

# compress the json file
get_ipython().system('zip -j ../data/doi_dnp.zip ../data/doi_dnp.json')


# In[2]:


# uncompress the json file
get_ipython().system('unzip ../data/doi_dnp.zip -d ../data/')

# Load the file, so you don't need to generate it
# every time
with open('../data/doi_dnp.json', 'r') as inputfile:
    doi = json.load(inputfile)


# In[3]:


len(doi)


# In[9]:

#z = 'https://pubs.acs.org/doi/pdf/10.1021/acs.jnatprod.8b00292' 
chemical_entities = [retrieve_chemical_entities(d, save=True) for d in doi]
fls = os.listdir()
fls = [x.replace('.pdf', '') for x in fls]

ndoi = [x for x in doi if x.split('/')[-1] in fls]
with open('remaining.json', 'w+') as f:
    json.dump(ndoi, f, indent=4)

chemical_entities = [retrieve_chemical_entities(d, save=True) for d in ndoi]
chemical_entities = retrieve_chemical_entities(doi[-3:])

for url in alldoi:
    doi = link.split('/')[-1]
    download_pdf(url, '%s.pdf' % doi)


# In[8]:
os.chdir('jnatprod/')
chemical_entities = [retrieve_chemical_entities(d) for d in doi[:5]]

os.chdir('jnatprod')

chemical_entities_images = [retrieve_chemical_entities_from_image(d) for d in doi[:5]]


# In[10]:


# is that correct???
chemical_entities


# In[ ]:
fls = [os.path.join('jnatprod_pdfs/', x) for x in os.listdir('jnatprod_pdfs/')]

entities = []
for f in fls[2016:]:
    try:
        entities.append(retrieve_chemical_entities(f))
    except:
        entities.append({f:{'oscar': []}})
    with open("entities.json", "w+") as jsonFile:
        json.dump(entities, jsonFile)

dflist = []
for e in entities:
    for k,v in e.items():
        if len(v['oscar'])==0:
            continue
        elif len(v['oscar'])>1:
            print('Warning, list with more than one element!')
            break
        else:
            tmp = pd.DataFrame(v['oscar'][0]['chemicalData']).T
            tmp['doi'] = k
            dflist.append(tmp)

df = pd.concat(dflist, sort=False)

# stats
df.reset_index(drop=True, inplace=True)
df.shape
len(df['standardInChIKey'].unique())
df['doi'].value_counts()
df['standardInChIKey'].value_counts()

exactMol = []

for x in df.standardInChI:
    try:
        exactMol.append(rdMD.CalcExactMolWt(Chem.MolFromInchi(x)))
    except:
        exactMol.append(0)

df['ExactMolWt'] = exactMol
df['ExactMolWt'].median()
df.to_csv('entities_dataframe.tsv', sep='\t', index=None)


nms = []
fls = os.listdir('gn_txt/')
for fl in fls:
    with open(os.path.join('gn_txt', fl)) as f:
        gn = json.load(f)

    if gn['names']==None:
        continue

    for nm in gn['names']:
        matchType = nm['verification']['BestResult']['matchType']
        if matchType == 'NoMatch':
            nms.append([fl.replace('_gn.txt', ''),
                        nm['verbatim'], nm['odds'],
                        '', '', '', '', matchType])
        elif 'classificationRank' in nm['verification']['BestResult'].keys():
            nms.append([fl.replace('_gn.txt', ''),
                        nm['verbatim'], nm['odds'],
                        nm['verification']['BestResult']['dataSourceId'],
                        nm['verification']['BestResult']['taxonId'],
                        nm['verification']['BestResult']['classificationPath'],
                        nm['verification']['BestResult']['classificationRank'],
                        matchType])
        else:
            nms.append([fl.replace('_gn.txt', ''),
                        nm['verbatim'], nm['odds'],
                        nm['verification']['BestResult']['dataSourceId'],
                        nm['verification']['BestResult']['taxonId'],
                        '', '', matchType])



dfnms = pd.DataFrame(nms)
dfnms.columns = ['doi', 'verbatim', 'odds', 'dataSourceId',
                 'taxonId', 'classificationPath',
                 'classificationRank', 'matchType']

sum(dfnms['taxonId']=='')
sum(dfnms['matchType']=='NoMatch')
dfnms['matchType'].unique()
len(dfnms['doi'].unique())
len(dfnms.loc[dfnms['matchType']!='NoMatch','taxonId'].unique())
sum(dfnms['classificationPath']=='')

dfnms.to_csv('gn_dataframe.tsv', sep='\t', index=None)

with open('jnatprod_pdfs/image_entity.json') as f:
    imgs = json.load(f)

imglist = []
for img in imgs:
    for k,v in img.items():
        for s in v['osra']:
            s = s.strip()
            m = Chem.MolFromSmiles(s)
            if m==None:
                imglist.append([k, s, '', 0])
            else:
                imglist.append([k, s, Chem.MolToInchiKey(m),
                                rdMD.CalcExactMolWt(m)])

dfimg = pd.DataFrame(imglist)
dfimg.columns = ['doi', 'smiles', 'standardInChIKey', 'ExactMolWt']

sum(dfimg['standardInChIKey']!='')
len(dfimg.loc[dfimg['standardInChIKey']!='', 'standardInChIKey'].unique())
len(dfimg.loc[dfimg['standardInChIKey']!='', 'doi'].unique())

dfimg.to_csv('img_dataframe.tsv', sep='\t', index=None)

sall = pd.concat([df[['doi', 'standardInChIKey']],
                  dfimg.loc[dfimg['standardInChIKey']!='', ['doi', 'standardInChIKey']]])

sall['doi'] = sall['doi'].str.replace('jnatprod_pdfs/', '')
sall = sall[~sall.duplicated()]
sall.reset_index(inplace=True, drop=True)

keys = sall['standardInChIKey'].unique()

dbid = []
for k in keys:
    dbid.append([k, inchikey2cid(k)])

dbids = pd.DataFrame(dbid)
dbids.fillna(0, inplace=True)
dbids[1] = dbids[1].astype(int)
dbids[1] = dbids[1].astype(str)
dbids.replace('0', '', inplace=True)
sum(dbids[1]=='')

dbids.to_csv('all_ids_pubchem_chemspider.tsv', sep='\t', index=None)

idfimg = pd.merge(dfimg, dbids, left_on='standardInChIKey', right_on=0, how='left')
idfimg.drop(0, axis=1, inplace=True)
idfimg.fillna('', inplace=True)
idfimg.rename(columns={1: 'pubchem'}, inplace=True)
idfimg.to_csv('img_dataframe.tsv', sep='\t', index=None)

idf = pd.merge(df, dbids, left_on='standardInChIKey', right_on=0, how='left')
idf.drop(0, axis=1, inplace=True)
idf.fillna('', inplace=True)
idf.rename(columns={1: 'pubchem'}, inplace=True)
idf.to_csv('entities_dataframe.tsv', sep='\t', index=None)

idf['doi'] = idf['doi'].str.replace('jnatprod_pdfs/', '')
idf['smiles'] = ''
idf['source'] = 'oscar'

idfimg['standardInChI'] = ''
idfimg['source'] = 'osra'

dfall = pd.concat([idf[idfimg.columns], idfimg])
dfall.to_csv('entities_img_dataframe.tsv', sep='\t', index=None)
