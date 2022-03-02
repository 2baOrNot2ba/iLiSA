import os
import yaml

USER_CONF_DIR = os.path.expanduser('~/.config/ilisa/')
USER_DATA_DIR = os.path.expanduser('~/.local/share/ilisa/')
USER_CACHE_DIR = os.path.expanduser('~/.cache/ilisa/')

LATESTDATAFILE = os.path.join(USER_CACHE_DIR, "latestdatafiles.txt")

def default_access_lclstn_conf():
    """Return the default access local station configuration."""
    accessconffile = os.path.join(USER_CONF_DIR, 'access_lclstn.conf')
    with open(accessconffile) as cfigfilep:
        accessconf = yaml.safe_load(cfigfilep)
    return accessconf
