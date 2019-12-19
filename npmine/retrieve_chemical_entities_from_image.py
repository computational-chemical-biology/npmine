#Import all necessary packages
import requests
from bs4 import BeautifulSoup
import os
import json
import subprocess
import itertools
from rdkit import Chem


def retrieve_chemical_entities_from_image(modification_doi):
     """Obtains images from DOI and retrieves chemical information
    Parameters
    ----------
    modification_doi: list
        List of DOIs from JNP.
    Returns
        List of chemical entities.
    -------
    """
    #From each pdf link, the pdf is downloaded and then a txt file is
    #created that contains the simplified structure of molecules
    #represented by images on the paper
    for a in modification_doi:
        subprocess.call(['wget', '-O', 'jnatprod.pdf', a])
        if os.path.isfile('jnatprod.pdf'):
            with open('jnpd','a') as i:
                subprocess.call("docker run berlinguyinca/osra", shell=True, stdout=i, stderr=i)
                subprocess.call("docker run -v $PWD:/home berlinguyinca/osra osra /home/jnatprod.pdf", shell=True, stdout=i, stderr=i)
                with open('osra_list','a') as i:
                    subprocess.call("docker run -v $PWD:/home berlinguyinca/osra osra /home/jnatprod.pdf -v -w /home/jnatprod.tx", shell=True, stdout=i, stderr=i)
                    os.remove('jnatprod.pdf')
                    os.remove('jnatprod.tx')
                    os.remove('osra_list')
                    os.remove('osra_list')
