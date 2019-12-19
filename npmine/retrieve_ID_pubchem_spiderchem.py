def retrieve_ID_pubchem_spiderchem(inchi):
     """Obtains ID of chemical entity from InCHI Keys 
    Parameters
    ----------
    inchi: string
        InCHI key obtained from InChIkeys list.
    Returns
        PubChem ID or SpiderChem ID.
    -------
    """
    consumer_key = '3uCJAFueyvOnnLG1yS9JHpZzuIGDaglX'

    def inchikey2csid(inchi):
        if any(inchi in s for s in InChIkeys):
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

    def inchikey2cid(inchi):
        if any(inchi in s for s in InChIkeys):
            gurl = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/%s/cids/JSON' % inchi
            r = requests.get(gurl)
            if r.status_code==200:
                return json.loads(r.text)['IdentifierList']['CID'][0]
        else:
            return inchikey2csid(inchi)
