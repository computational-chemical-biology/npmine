#Import all necessary packages
import os
import json
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
    headers = {'User-Agent': 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:74.0) Gecko/20100101 Firefox/74.0'}
    r = requests.get(url, headers=headers)
    if r.status_code==200:
        with open(filename, 'wb') as fd:
            fd.write(r.content)
    else:
        r.raise_for_status()


def retrieve_chemical_entities(filename):
    """Obtain chemical information
    Parameters
    ----------
    filename: str
        pdf file name.
    Returns
        Chemical entity dictionary.
    -------
    """

    #From each pdf link, the pdf is downloaded and then OSCAR is used to
    #transform the file into json and obtain the chemical entities
    filename = filename.replace('.pdf', '')
    os.system('%s %s.pdf > %s.json' % (config['TOOLS']['OSCAR'], filename, filename))
    if os.path.isfile('%s.json' % filename):
        with open('%s.json' % filename) as f:
            jnatprod = json.load(f)
            os.remove('%s.json' % filename)
        return {filename:{'oscar':jnatprod}}
    else:
        return {filename:{'oscar': []}}
