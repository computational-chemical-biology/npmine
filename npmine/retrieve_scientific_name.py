import requests
import os
import json
import configparser
import npmine

fullname = os.path.join(npmine.__path__[0], '..', 'config.ini')

config = configparser.ConfigParser()
config.read(fullname)

def html2txt(url, filename):
    """ Downloads full text html and saves to text file
    Parameters
    ----------
    url: str
        Any url for open full text.
    Returns
        Text file.
    -------
    """
    # example: https://bmcchem.biomedcentral.com/articles/10.1186/s13065-018-0497-z 
    r = requests.get(url)
    htmlstring = r.content
    gnfinder = htmlstring.decode('utf-8','ignore')

    f = open(filename, 'w' )
    f.write(gnfinder)
    f.close()

def retrieve_scietific_name(input_file, out_file):
    """retrieves scientific names present in the document
    Parameters
    ----------
    input_file: str
        pdf (with extension .pdf) or text file (with extension .txt) name.
    out_file: str
        output file name.
    Returns
        json file.
    -------
    """
    if '.pdf' in input_file:
        cmd = '%s -layout %s'
        cmd = cmd % (config['TOOLS']['PDFTOTEXT'], input_file)
        os.system(cmd)
        input_file = input_file.replace('.pdf', '.txt')

    cmd = "%s find -c -l eng -s 4,12 %s > %s"
    cmd = cmd % (config['TOOLS']['GNFINDER'], input_file, out_file)
    os.system(cmd)
    os.remove(input_file)
