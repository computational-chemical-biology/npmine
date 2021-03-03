import json
import re
import os
import requests

from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import numpy as np
from IPython.display import SVG
import matplotlib.pyplot as plt

html = '''<!DOCTYPE html>
<html>
  <head>
    <meta http-equiv="content-type" content="text/html; charset=UTF-8">
    <title>NPMINE</title>
    <script src="https://code.jquery.com/jquery-3.3.1.js" integrity="sha256-2Kok7MbOyxpgUVvAk/HJ2jigOSYS2auK4Pfzbm7uH60=" crossorigin="anonymous"></script>
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js" integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>
    <link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.6.1/css/all.css" integrity="sha384-gfdkjb5BdAXd+lj+gudLWI+BXq4IuLW5IT+brZEZsLFm++aCMlF1V92rMkPaX4PP" crossorigin="anonymous">
</head>

  <body>
    <!-- A grey horizontal navbar that becomes vertical on small screens -->
    <nav class="navbar navbar-expand-sm bg-light navbar-light">
        <a class="navbar-brand" href="">
        </a>

      <!-- Links-->
      <ul class="navbar-nav">
        <li class="nav-item">
          <a class="nav-link" href="http://ccbl.fcfrp.usp.br">CCBL</a>
        </li>
      </ul>

    </nav>
    <div id="content">

    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js" integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>
    <link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.6.1/css/all.css" integrity="sha384-gfdkjb5BdAXd+lj+gudLWI+BXq4IuLW5IT+brZEZsLFm++aCMlF1V92rMkPaX4PP" crossorigin="anonymous">

    <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>

    <script src="https://cdn.datatables.net/1.10.4/js/jquery.dataTables.min.js"></script>
    <link rel="stylesheet" href="https://cdn.datatables.net/1.10.4/css/jquery.dataTables.min.css">
    <script>
        var dataSet = REPLACE;
        console.log(dataSet);
        $(document).ready(function() {
        $('#resTable').DataTable( {
            data : dataSet,
            // add column definitions to map your json to the table
            "columns": [
            {title: "doi"},
            {title: "pubchem"},
            {title: "ExactMolWt"},
            {title: "smiles"},
            {title: "source"},
            {title: "Structure"}
            ]
        } );
        });
    </script>
    <div class="m-5">
        <table id="resTable" class="table table-striped" style="width:100%" >
		<thead>
        <tr>
		<th>DOI</th>
		<th>PubChem</th>
		<th>ExactMolWt</th>
		<th>SMILES</th>
		<th>Source</th>
		<th>Structure</th>
        </tr>
        </thead>
        </table>
    </div>

    </div>
  </body>
</html>'''

def inchi2smiles(inchi):
    return Chem.MolToSmiles(Chem.MolFromInchi(inchi))

def format_source_table(oscar_list, osra_list):
    df_oscar = pd.concat(oscar_list)
    df_oscar['source'] = 'oscar'
    df_oscar['pubchem'] = [inchikey2cid(x) for x in df_oscar['standardInChIKey']]

    df_osra = pd.concat(osra_list)
    df_osra['source'] = 'osra'
    df_oscar['smiles'] = [inchi2smiles(x) for x in df_oscar['standardInChI']]

    pubchem = []
    for x in df_osra['standardInChIKey']:
        cid = inchikey2cid(x)
        if not cid is None:
            pubchem.append(cid)
        else:
            pubchem.append(0)

    df_osra['pubchem'] = pubchem

    report_print = pd.concat([df_oscar[['doi', 'pubchem', 'ExactMolWt', 'smiles','source']],
                              df_osra[['doi', 'pubchem', 'ExactMolWt', 'smiles','source']]]).reset_index(drop=True)
    return report_print

def create_report(report_print, dois, out_file='npmine_report.html',
                  useSVG=False):
    """Creates an html report from NPMINE's results
    Parameters
    ----------
    report_print: pd.DataFrame
        DataFrame containing columns 'doi', 'pubchem', 'ExactMolWt', 'smiles','source'.
    dois: str or list
        Directory of doi link files or link list.
    out_file: str
        HTML output file name.
    useSVG: bool
        If svg format should be used. Default is png.
    Returns
        Report html file.
    -------
    """
    if type(dois)==str:
        fls = [x for x in os.listdir(dois) if '.json' in x]
        dois_list = []
        for fl in fls:
            with open(os.path.join(dois, fl)) as f:
                dois_list.extend(json.load(f))
        dois = pd.DataFrame(dois_list)
    else:
        dois = pd.DataFrame(dois)

    if not os.path.exists('figs'):
        os.mkdir('figs')

    for i in report_print.index:
        ms = Chem.MolFromSmiles(report_print.loc[i, 'smiles'])
        if useSVG:
            try:
                img = Draw.MolsToGridImage(ms, subImgSize=(200,200),
                                           useSVG=True)
                with open('figs/%s.svg' % report_print.loc[i, 'pubchem'], 'w') as svgout:
                    svgout.write(img.data)
            except:
                fig1 = plt.gcf()
                fig1.savefig('figs/%s.png' % report_print.loc[i, 'pubchem'])
        else:
            try:
                Draw.MolToFile(ms, 'figs/%s.png' % report_print.loc[i, 'pubchem'], subImgSize=(200,200))
            except:
                fig1 = plt.gcf()
                fig1.savefig('figs/%s.png' % report_print.loc[i, 'pubchem'])

    report_print['Structure'] = ''
    for i in report_print.index:
        link = dois.loc[dois[0].str.contains(report_print.loc[i, 'doi']), 0].values[0]
        report_print.loc[i, 'doi'] = '<a href="{}" target="_blank">{}</a>'.format(*[link, report_print.loc[i, 'doi']])
        if useSVG:
            fig = "./figs/%s.svg" % report_print.loc[i, 'pubchem']
        else:
            fig = "./figs/%s.png" % report_print.loc[i, 'pubchem']
        report_print.loc[i, 'Structure'] = '<img src="%s" width="120" height="120">' % fig
        report_print.loc[i, 'pubchem'] = '<a href="https://pubchem.ncbi.nlm.nih.gov/compound/{}" target="_blank">{}</a>'.format(*[int(report_print.loc[i, 'pubchem']), int(report_print.loc[i, 'pubchem'])])
        smi = requests.utils.quote(report_print.loc[i, 'smiles'])
        report_print.loc[i, 'smiles'] = '<a href="https://pubchem.ncbi.nlm.nih.gov/edit3/index.html?smiles={}" target="_blank">{}</a>'.format(*[smi, 'PubChem Sketcher'])

    html_local = re.sub('REPLACE', json.dumps(report_print.apply(lambda a: a.tolist(), axis=1).tolist()), html)
    with open(out_file, 'w+') as f:
        f.write(html_local)

