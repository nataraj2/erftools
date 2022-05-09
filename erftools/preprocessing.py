import f90nml

from .WRFdefs import time_control, domains#, physics
from .inputs import ERFinput

class WRFnamelist(object):
    """Class to parse WRF inputs and convert to ERF inputs"""

    def __init__(self,nmlpath):
        with open(nmlpath,'r') as f:
            self.nml = f90nml.read(f)
        self.time_control = time_control(self.nml['time_control'])
        self.domains = domains(self.nml['domains'])
        self.erf_input = ERFinput()
        self.calculate_inputs()

    def __str__(self):
        s = str(self.time_control) + '\n'
        s+= str(self.domains)+ '\n'
        return s

    def calculate_inputs(self):
        tdelta = self.time_control.end_datetime - self.time_control.start_datetime
        self.erf_input['stop_time'] = tdelta.total_seconds()
        n_cell = [self.domains.e_we[0], self.domains.e_sn[0], self.domains.e_vert[0]]
        self.erf_input['amr.n_cell'] = n_cell
        
