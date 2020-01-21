#Import all necessary packages
import requests
from bs4 import BeautifulSoup
import os
import json
import subprocess
import itertools
from rdkit import Chem
import configparser
import npmine

fullname = os.path.join(npmine.__path__[0], '..', 'config.ini')

config = configparser.ConfigParser()
config.read(fullname)

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
    for z in list_of_dois:
        subprocess.call(['wget', '-O', 'jnatprod.pdf', z])
        os.system('%s jnatprod.pdf > jnatprod.json' % config['TOOLS']['OSCAR'])
        if os.path.isfile('jnatprod.json'):
            with open('jnatprod.json') as f:
                jnatprod = json.load(f)
                os.remove('jnatprod.pdf')
                os.remove('jnatprod.json')
            json_list_doi.append({z:{'oscar':jnatprod}})
    return json_list_doi
