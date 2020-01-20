#Import all necessary packages
import requests
from bs4 import BeautifulSoup
import os
import json
import subprocess
import itertools
from rdkit import Chem

# Use the parameter journal_name to retrieve doi list from
# different journals, starting with ACS group
def retrieve_doi(journal_name='all'):
    """Performs web scraping to obtain DOIs from journal of ACS
    Parameters
    ----------
    journal_name: str
        Dictionary containing node attributes and edge list from GNPS.
    Returns
        List containing DOIs.
    -------
    """
    #Request html and use BeautifulSoup to parse table of contents page
    #from Journal of Natural Products website
    page = requests.get("https://pubs.acs.org/loi/jnprdf")
    soup = BeautifulSoup(page.text, "html.parser")

    soup.find_all(text=True)

    #From the html of the table of contents, obtain all the html links that
    #correspond to each issue
    toc_list = []
    for a in soup.find_all("a", href=True):
        if "/toc/jnprdf/" in a["href"]:
            toc_list.append(a["href"])

    #From the list of all issues, order them in sequence
    new_list = list(set(toc_list))
    new_list.sort()

    #From each html link, obtain the corresponding doi of each paper of
    #that issue 
    all_dois = {}
    for i in new_list[1:]:
        html1 = "https://pubs.acs.org" + i
        pag = requests.get(html1)
        parser = BeautifulSoup(pag.text, "html.parser")
        doi_list = []
        for a in parser.find_all("a", href=True):
            if "/doi/pdf/10.1021/" in a["href"]:
                doi_list.append(a["href"])
        all_dois[i] = doi_list

    #To be able to get the pdf file, it is necessary to transform into a 
    #list again and then do a string concatenation of each doi
    list_values = [v for v in dict.values(all_dois)]
    merged = list(itertools.chain(*list_values))
    modification_doi = ["https://pubs.acs.org" + i for i in merged]
    return modification_doi
