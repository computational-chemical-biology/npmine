#Import all necessary packages
import requests
from bs4 import BeautifulSoup
import os
import json
import subprocess
import itertools
from rdkit import Chem


def retrieve_chemical_entities(list_of_dois):
    """Obtain publication from DOI and retrieves chemical information
    Parameters
    ----------
    list_of_dois: list
        List of DOIs from JNP.
    Returns
        List of chemical entities.
    -------
    """

    #From each pdf link, the pdf is downloaded and then OSCAR is used to
    #transform the file into json and obtain the chemical entities
    json_list_doi = []
    for z in list_of_doi:
        subprocess.call(['wget', '-O', 'jnatprod.pdf', z])
        os.system('oscarpdf2json jnatprod.pdf > jnatprod.json')
        if os.path.isfile('jnatprod.json'):
            with open('jnatprod.json') as f:
                jnatprod = json.load(f)
                os.remove('jnatprod.pdf')
                os.remove('jnatprod.json')
            json_list_doi.append({z:{'oscar':jnatprod}})
    return json_list_doi

    #A list with all the InCHI keys is created so that the functions below
    #can be used to obtain either the PubChem ID or the SpiderChem ID of 
    #each chemical entity
    InChIkeys = []
    dictionary = json_list_doi[1]['https://pubs.acs.org/doi/pdf/10.1021/np50001a001']
    for i in dictionary.values():
        for j in i:
            for k in j["chemicalData"].values():
                InChIkeys.append(k['standardInChIKey'])



