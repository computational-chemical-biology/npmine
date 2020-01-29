#Import all necessary packages
import os
import configparser
import npmine

fullname = os.path.join(npmine.__path__[0], '..', 'config.ini')

config = configparser.ConfigParser()
config.read(fullname)

def retrieve_chemical_entities_from_image(url):
    """Obtains images from DOI and retrieves chemical information
    Parameters
    ----------
    url: str
        url containing pdf file name after '/' or pdf file name.
    Returns
        List of chemical entities.
    -------
    """
    doi = url.split('/')[-1] # for jnatprod
    os.system(config['TOOLS']['OSRA'].format(*[doi]*2))

    if os.path.isfile('%s.txt' % doi):
        with open('%s.txt' % doi) as f:
            jnatprod = f.readlines()
            os.remove('%s.txt' % doi)
        return {doi:{'osra':jnatprod}}
    else:
        return {doi:{'osra': []}}
