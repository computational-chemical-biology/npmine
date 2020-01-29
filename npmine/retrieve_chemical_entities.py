#Import all necessary packages
import requests
from bs4 import BeautifulSoup
import os
import json
import subprocess
import itertools
from rdkit import Chem
import configparser
import requests
import npmine

fullname = os.path.join(npmine.__path__[0], '..', 'config.ini')

config = configparser.ConfigParser()
config.read(fullname)

def download_pdf(url, filename):
    """Download pdf
    Parameters
    ----------
    url: str
        Downloadable pdf link.
    Returns
        pdf file.
    -------
    """
    r = requests.get(url)
    if r.status_code!=200:
        with open(filename, 'wb') as fd:
            fd.write(r.content)
    else:
        r.raise_for_status()


def retrieve_chemical_entities(url):
    """Obtain chemical information
    Parameters
    ----------
    url: str
        url containing pdf file name after '/' or pdf file name.
    Returns
        Chemical entity dictionary.
    -------
    """

    #From each pdf link, the pdf is downloaded and then OSCAR is used to
    #transform the file into json and obtain the chemical entities
    doi = url.split('/')[-1] # for jnatprod
    os.system('%s %s.pdf > %s.json' % (config['TOOLS']['OSCAR'], doi, doi))
    if os.path.isfile('%s.json' % doi):
        with open('%s.json' % doi) as f:
            jnatprod = json.load(f)
            os.remove('%s.json' % doi)
        return {doi:{'oscar':jnatprod}}
    else:
        return {doi:{'oscar': []}}
