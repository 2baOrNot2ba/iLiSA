import threading
import plumbum

class DRUinterface():
    def __init__(self, accessconf_dru):
        hostname = accessconf_dru['hostname']
        user = accessconf_dru['user']
        self.dru = plumbum.SshMachine(hostname, user=user)
    
    def rec_bf_proxy(self, starttime, duration, lanes, band_, bf_data_dir,
                     port0, stnid_):
        recorders = ['ow', 'py']
        which_recorder = recorders[0]
        rec_cmd = self.dru['rec_bf_streams']
        cli_arg = "-s {} -d {} -l {} -b {} -p {} -w {}".format(
            starttime, duration, lanes, bf_data_dir, port0, which_recorder)
        rec_cmd(cli_arg)
    