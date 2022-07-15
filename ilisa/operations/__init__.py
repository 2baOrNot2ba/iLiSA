import os
import time
import yaml
import logging

USER_CONF_DIR = os.path.expanduser('~/.config/ilisa/')
USER_DATA_DIR = os.path.expanduser('~/.local/share/ilisa/')
USER_CACHE_DIR = os.path.expanduser('~/.cache/ilisa/')

DATETIMESTRFMT = '%Y-%m-%dT%H:%M:%S'


def default_access_lclstn_conf(stnid=None):
    """
    Get local default access configuration to station

    Parameters
    ----------
    stnid : str, optional
        If stnid is given, try to get the access config from the file in the
        USER_CONF_DIR called 'access_{stnid}.conf'. If it not set then
        try to get file named 'access_lclstn.conf'.

    Returns
    -------
    accessconf : dict
        Access configuration dict, which contains info on LCU and DRU access.
    """
    if stnid==None:
        stnsuf = 'lclstn'
    else:
        stnsuf = stnid
    accessconffilename = 'access_'+stnsuf+'.conf'
    accessconffile = os.path.join(USER_CONF_DIR, accessconffilename)
    try:
        with open(accessconffile) as cfigfilep:
            accessconf = yaml.safe_load(cfigfilep)
    except IOError:
        accessconf = None

    return accessconf

logging.Formatter.converter = time.gmtime
logging.basicConfig(
    # filename='ilisa_operations.log',
    format='%(asctime)s.%(msecs)03d | %(levelname)s | %(message)s',
    datefmt=DATETIMESTRFMT,
    level=logging.INFO
)
