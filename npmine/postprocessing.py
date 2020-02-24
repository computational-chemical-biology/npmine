import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rdMD

def get_exact_mass(inchi):
    """Download pdf
    Parameters
    ----------
    url: str
        InChI molecular representation.
    Returns
        Exact mass.
    -------
    """
    try:
        return rdMD.CalcExactMolWt(Chem.MolFromInchi(inchi))
    except:
        return 0

def entity_dict2dataframe(entities, exact_mass=True):
    """Converts entity dictionary obtained with oscar
       to pandas DataFrame
    Parameters
    ----------
    entities: dict
        Dictiorary containing chemicalData.
    exact_mass: bool
        Whether to calculate exact mass.
    Returns
        pandas DataFrame.
    -------
    """

    dflist = []
    for k,v in entities.items():
        if len(v['oscar'])==0:
            continue
        elif len(v['oscar'])>1:
            raise Exception('Warning, list with more than one element!')
        elif v['oscar'][0]['chemicalData']=={}:
            continue
        else:
            tmp = pd.DataFrame(v['oscar'][0]['chemicalData']).T
            tmp['doi'] = k
            dflist.append(tmp)

    if not len(dflist):
        return None

    df = pd.concat(dflist, sort=False)

    if exact_mass:
        exactMol = []

        for x in df.standardInChI:
            exactMol.append(get_exact_mass(x))

        df['ExactMolWt'] = exactMol

    return df

def sci_name_dict2dataframe(gn, identifier):
    """Converts scientific species dictionary obtained with gnfinder
       to pandas DataFrame
    Parameters
    ----------
    gn: dict
        Dictiorary containing scientific names.
    identifier: str
        Identifier for the dict data.
    Returns
        pandas DataFrame.
    -------
    """
    nms = []
    for nm in gn['names']:
        matchType = nm['verification']['BestResult']['matchType']
        if matchType == 'NoMatch':
            nms.append([identifier,
                        nm['verbatim'], nm['odds'],
                        '', '', '', '', matchType])
        elif 'classificationRank' in nm['verification']['BestResult'].keys():
            nms.append([identifier,
                        nm['verbatim'], nm['odds'],
                        nm['verification']['BestResult']['dataSourceId'],
                        nm['verification']['BestResult']['taxonId'],
                        nm['verification']['BestResult']['classificationPath'],
                        nm['verification']['BestResult']['classificationRank'],
                        matchType])
        else:
            nms.append([identifier,
                        nm['verbatim'], nm['odds'],
                        nm['verification']['BestResult']['dataSourceId'],
                        nm['verification']['BestResult']['taxonId'],
                        '', '', matchType])

    dfnms = pd.DataFrame(nms)
    dfnms.columns = ['doi', 'verbatim', 'odds', 'dataSourceId',
                     'taxonId', 'classificationPath',
                     'classificationRank', 'matchType']
    return dfnms

def image_dict2dataframe(img, exact_mass=True):
    """Converts image dictionary obtained with osra
       to pandas DataFrame
    Parameters
    ----------
    img: dict
        Dictiorary containing chemical data.
    exact_mass: bool
        Whether to calculate exact mass.
    Returns
        pandas DataFrame.
    -------
    """

    imglist = []
    for k,v in img.items():
        for s in v['osra']:
            s = s.strip()
            imglist.append([k, s])

    dfimg = pd.DataFrame(imglist)
    dfimg.columns = ['doi', 'smiles']

    if exact_mass:
        dfimg['standardInChIKey'] = ''
        dfimg['ExactMolWt'] = 0

        for i in dfimg.index:
            s = dfimg.loc[i, 'smiles']
            m = Chem.MolFromSmiles(s)
            if m!=None:
                dfimg.loc[i, 'standardInChIKey'] = Chem.InchiToInchiKey(Chem.MolToInchi(m))
                dfimg.loc[i, 'ExactMolWt'] = rdMD.CalcExactMolWt(m)

    return dfimg

