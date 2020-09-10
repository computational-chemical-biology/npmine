import requests

def get_bioactivity(queries, chunksize=50):
    """Accesses PubChem's API and queries for bioactivity
    Parameters
    ----------
    queries: list
        InChIKey list.
    chunksize: int
        Length of list chunks to submit.
    Returns
        List containing bioactivity results.
    ----------
    """
    jlist = []
    for i in range(0, len(queries), chunksize):
        to_query = queries[i:i+chunksize]
        ids = ','.join(map(str, to_query.tolist()))
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/assaysummary/JSON' % ids
        rq = requests.get(url)
        jlist.append(rq.json())

    return jlist
