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

