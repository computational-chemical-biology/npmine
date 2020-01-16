import requests
from bs4 import BeautifulSoup
import os
import json
import subprocess
import itertools

def retrieve_scietific_name():
    """Obtain publication from DOI and retrieves scientific names present
    Parameters
    ----------
    Returns
        Scientific Name.
    -------
    """
    page = requests.get("https://bmcchem.biomedcentral.com/articles?searchType=journalSearch&sort=PubDate&page=1")
    soup = BeautifulSoup(page.text, "html.parser")

    soup.find_all(text=True)

    toc_list = []
    for a in soup.find_all("a", href=True):
        if "/articles?searchType=journalSearch" in a["href"]:
            toc_list.append(a["href"])

    new_list = list(set(toc_list))
    new_list.sort()

    all_dois = {}
    for i in new_list[:3]:
        html1 = "https://bmcchem.biomedcentral.com" + i
        pag = requests.get(html1)
        parser = BeautifulSoup(pag.text, "html.parser")
        doi_list = []
    for a in parser.find_all("a", href=True):
        if "/articles/10.1186/" in a["href"]:
            doi_list.append(a["href"])
    all_dois[i] = doi_list

    list_values = [v for v in dict.values(all_dois)]
    merged = list(itertools.chain(*list_values))
    modification_doi = ["https://bmcchem.biomedcentral.com" + i for i in merged]
    list(set(modification_doi))

    r = requests.get('https://bmcchem.biomedcentral.com/articles/10.1186/s13065-018-0497-z')
    htmlstring = r.content
    gnfinder = htmlstring.decode('utf-8','ignore')

    f = open('gnfinder_input.txt', 'w' )
    f.write(gnfinder)
    f.close()
    subprocess.call("gnfinder find -c -l eng -s 4,12 gnfinder_input.txt > gnfinder_result.txt", shell=True)
