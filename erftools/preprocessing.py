import numpy as np
import xarray as xr
import f90nml

from .wrf_namelist_input import time_control, domains, physics
from .inputs import ERFInputFile

class WRFInputDeck(object):
    """Class to parse inputs from WRF and convert to inputs for ERF
    WRF inputs include:
    * namelist.input
    * wrfinput_d01[, wrfinput_d02, ...]
    """

    def __init__(self,nmlpath):
        with open(nmlpath,'r') as f:
            self.nml = f90nml.read(f)
        self.time_control = time_control(self.nml['time_control'])
        self.domains = domains(self.nml['domains'])
        self.physics = physics(self.nml['physics'])
        self.erf_input = ERFInputFile()
        self.calculate_inputs()

    def __str__(self):
        s = str(self.time_control) + '\n'
        s+= str(self.domains)+ '\n'
        s+= str(self.physics)+ '\n'
        return s

    def calculate_inputs(self):
        """Scrape inputs for ERF from a WRF namelist.input file

        Notes:
        * This does not import the WRF height levels (geopotential height); to
          get those coordinates, call get_heights() to extract from WRF initial
          conditions (`wrfinput_d01` from real.exe)
        * This does not import z0 (surface roughness) from the reanalysis data;
          to get that, call get_roughness_map()
        """
        tdelta = self.time_control.end_datetime - self.time_control.start_datetime
        self.erf_input['stop_time'] = tdelta.total_seconds()

        # note: number vert pts is staggered
        # TODO: add vertical stretching
        n_cell = np.array([self.domains.e_we[0], self.domains.e_sn[0], self.domains.e_vert[0]-1])
        ztop_est = 287.0 * 300.0 / 9.81 * np.log(1e5/self.domains.p_top_requested)
        self.erf_input['geometry.prob_extent'] = n_cell * np.array([self.domains.dx[0], self.domains.dy[0], np.nan])
        self.erf_input['geometry.prob_extent'][2] = ztop_est
        self.erf_input['amr.n_cell'] = n_cell

        # TODO: verify that refined regions will take finer time steps
        dt = np.array(self.domains.parent_time_step_ratio) * self.domains.time_step
        self.erf_input['erf.fixed_dt'] = dt[0]

        self.erf_input['amr.max_level'] = self.domains.max_dom - 1 # zero-based indexing
        grid_ratio = self.domains.parent_grid_ratio[1] # TODO: assume all nests have same ratio
        self.erf_input['amr.ref_ratio_vect'] = [grid_ratio, grid_ratio, 1]
        refine_names = ' '.join([f'box{idom:d}' for idom in range(1,self.domains.max_dom)])
        self.erf_input['amr.refinement_indicators'] = refine_names
        for idom in range(1,self.domains.max_dom):
            # zero-based indexing -- TODO: verify that these are relative to the parent box, not level0
            in_box_lo = [self.domains.i_parent_start[idom]-1     , self.domains.j_parent_start[idom]-1     ]
            # TODO: verify that these are incusive bounds
            in_box_hi = [in_box_lo[0] + self.domains.e_we[idom]-1, in_box_lo[1] + self.domains.e_sn[idom]-1]
            self.erf_input[f'amr.box{idom:d}.in_box_lo'] = in_box_lo
            self.erf_input[f'amr.box{idom:d}.in_box_hi'] = in_box_hi

        restart_period = self.time_control.restart_interval * 60.0 # [s]
        self.erf_input['amr.check_int'] = int(restart_period / dt[0])

        # TODO: specify PBL scheme per level
        self.erf_input['erf.pbl_type'] = self.physics.bl_pbl_physics[0]
        
        return self.erf_input
        
    def get_heights(self):
        # Note: This calculates the heights for the outer WRF domain only
        wrfinp = xr.open_dataset('wrfinput_d01')
        ph = wrfinp['PH'] # perturbation geopotential
        phb = wrfinp['PHB'] # base-state geopotential
        hgt = wrfinp['HGT'] # terrain height
        self.terrain = hgt
        geo = ph + phb # geopotential, dims=(Time: 1, bottom_top_stag, south_north, west_east)
        geo = geo/9.81 - hgt
        geo = geo.isel(Time=0).mean(['south_north','west_east']).values
        self.heights = (geo[1:] + geo[:-1]) / 2 # destaggered
        # TODO: update erf_input with specified height levels
        return self.heights

    def get_roughness_map(self):
        wrfinp = xr.open_dataset('wrfinput_d01')
        

