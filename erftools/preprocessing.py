import numpy as np
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

        # note: number vert pts is staggered
        # TODO: add vertical stretching
        n_cell = np.array([self.domains.e_we[0], self.domains.e_sn[0], self.domains.e_vert[0]-1])
        self.erf_input['geometry.prob_extent'] = n_cell * np.array([self.domains.dx[0], self.domains.dy[0], np.nan])
        # see https://github.com/a2e-mmc/mmctools/blob/dev/mmctools/wrf/preprocessing.py#L3336
        self.erf_input['geometry.prob_extent'][2] = 287.*300./9.81 * np.log(1e5/self.domains.p_top_requested)
        self.erf_input['amr.n_cell'] = n_cell

        # TODO: verify that refined regions will take finer time steps
        dt = np.array(self.domains.parent_time_step_ratio) * self.domains.time_step
        self.erf_input['erf.fixed_dt'] = dt[0]
        
        
