import f90nml

from .WRFdefs import time_control#, domains, physics

class WRFnamelist(object):
    """Class to parse WRF inputs and convert to ERF inputs"""

    def __init__(self,nmlpath):
        with open(nmlpath,'r') as f:
            self.nml = f90nml.read(f)
        self.time_control = time_control(self.nml['time_control'])

    def __str__(self):
        s = ''
        s+= str(self.time_control)
        return s

