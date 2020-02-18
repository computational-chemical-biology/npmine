#Import all necessary packages
import os
import configparser
import npmine

fullname = os.path.join(npmine.__path__[0], '..', 'config.ini')

config = configparser.ConfigParser()
config.read(fullname)

def retrieve_chemical_entities_from_image(filename):
    """Obtains images from DOI and retrieves chemical information
    Parameters
    ----------
    filename: str
        pdf file name.
    Returns
        List of chemical entities.
    -------
    """
    filename = filename.replace('.pdf', '') # for jnatprod
    os.system(config['TOOLS']['OSRA'].format(*[filename]*2))
    if os.path.isfile('%s.txt' % filename):
        with open('%s.txt' % filename) as f:
            jnatprod = f.readlines()
            os.remove('%s.txt' % filename)
        return {filename:{'osra':jnatprod}}
    else:
        return {filename:{'osra': []}}
