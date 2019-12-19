# npmine 
Documentation on the npmine package

# Installation
```
conda create -n nplibrary   
conda install -n nplibrary rdkit -c rdkit 
pip install bs4 
source activate nplibrary
python setup.py install 
```
# Getting started

Load python terminal and load packages

```
from npmine.retrieve_doi import retrieve_doi 
dois = retrieve_doi() 
```

