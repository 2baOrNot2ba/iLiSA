import os
import yaml

user_conf_dir = os.path.expanduser('~/.config/ilisa/')
user_data_dir = os.path.expanduser('~/.local/share/ilisa/')

def default_access_lclstn_conf():
    """Return the default access local station configuration."""
    accessconffile = os.path.join(user_conf_dir, 'access_lclstn.conf')
    with open(accessconffile) as cfigfilep:
        accessconf = yaml.load(cfigfilep)
    return accessconf
