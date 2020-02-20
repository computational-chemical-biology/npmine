import json
import requests

fullname = os.path.join(npmine.__path__[0], '..', 'config.ini')

config = configparser.ConfigParser()
config.read(fullname)


def inchikey2csid(inchi, consumer_key=''):
    """Obtains ID of chemical entity from InCHIKeys
    Parameters
    ----------
    inchi: string
        InCHIkey identifier.
    consumer_key: string
         ChemSpider consumer key.
    Returns
        ChemSpider ID.
    -------
    """
    if consumer_key=='':
        consumer_key = config['TOOLS']['CONSUMER_KEY']
        if consumer_key==None:
            raise Exception('Please set the consumer key!')
    data = {"inchikey": inchi}
    headers = {"apikey": consumer_key}
    purl = 'https://api.rsc.org/compounds/v1/filter/inchikey'
    r = requests.post(purl, data=json.dumps(data), headers=headers)
    if r.status_code==200:
        queryId = json.loads(r.text)['queryId']
        surl = 'https://api.rsc.org/compounds/v1/filter/%s/status' % queryId
        rg = requests.get(surl, headers=headers)
        while json.loads(rg.text)['status']!='Complete':
            rg = requests.get(surl, headers=headers)
        rurl = 'https://api.rsc.org/compounds/v1/filter/%s/results' % queryId
        rr = requests.get(rurl, headers=headers)
        res = json.loads(rr.text)['results']
        if len(res):
            return 'csid%s' % res[0]
        else:
            return None
    else:
        return None

def inchikey2cid(inchi, consumer_key=''):
    """Obtains ID of chemical entity from InCHIKeys
    Parameters
    ----------
    inchi: string
        InCHIkey identifier.
    consumer_key: string
         ChemSpider consumer key.
    Returns
        PubChem ID.
    -------
    """

    gurl = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/%s/cids/JSON' % inchi
    r = requests.get(gurl)
    if r.status_code==200:
        return json.loads(r.text)['IdentifierList']['CID'][0]
    else:
        return inchikey2csid(inchi)
