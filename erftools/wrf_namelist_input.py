"""
Processing for each namelist within a WRF namelist.input file
"""
from datetime import datetime

class TimeControl(object):
    """&time_control namelist"""

    def __init__(self,nmldict):
        self.nml = nmldict
        self.restart_interval = self.nml['restart_interval'] # [min]
        self.parse_datetime_range()

    def __str__(self):
        return f"""WRF `time_control` parameters
  sim range: {self.start_datetime} to {self.end_datetime}"""

    def parse_datetime_range(self):
        # assume all domains have the same start/end datetime for now; TODO: allow for variable start times per domain
        idom = 0
        start_year  = self.nml['start_year'][idom]
        start_month = self.nml['start_month'][idom]
        start_day   = self.nml['start_day'][idom]
        start_hour  = self.nml['start_hour'][idom]
        end_year  = self.nml['end_year'][idom]
        end_month = self.nml['end_month'][idom]
        end_day   = self.nml['end_day'][idom]
        end_hour  = self.nml['end_hour'][idom]
        self.start_datetime = datetime(start_year, start_month, start_day, start_hour)
        self.end_datetime = datetime(end_year, end_month, end_day, end_hour)
        
    
class Domains(object):
    """&domains namelist"""

    def __init__(self,nmldict):
        self.nml = nmldict
        self.parse_time_integration()
        self.parse_grid()

    def __str__(self):
        s = 'WRF `domains` parameters\n'
        for dom in range(self.max_dom):
            s += f'  d{dom+1:02d}:' \
                 f' [(1,{self.e_we[dom]:d}),(1,{self.e_sn[dom]:d}),(1,{self.e_vert[dom]:d})]' \
                 f' ds=[{self.dx[dom]},{self.dy[dom]}]\n'
        return s.rstrip()

    def parse_time_integration(self):
        self.time_step = self.nml['time_step'] # seconds
        self.parent_time_step_ratio = self.nml['parent_time_step_ratio']

    def parse_grid(self):
        self.max_dom = self.nml['max_dom']
        # - assume s_we == s_sn == s_vert == [1, 1, ...]
        self.e_we = self.nml['e_we'][:self.max_dom] # west--east, unstaggered
        self.e_sn = self.nml['e_sn'][:self.max_dom] # south--north, unstaggered
        self.e_vert = self.nml['e_vert'][:self.max_dom] # bottom--top, STAGGERED
        self.dx = self.nml['dx'][:self.max_dom]
        self.dy = self.nml['dy'][:self.max_dom]
        self.p_top_requested = self.nml['p_top_requested']
        self.i_parent_start = self.nml['i_parent_start'][:self.max_dom]
        self.j_parent_start = self.nml['j_parent_start'][:self.max_dom]
        self.parent_grid_ratio = self.nml['parent_grid_ratio'][:self.max_dom]
        for dom in range(1,self.max_dom):
            assert (self.dx[dom-1]/self.dx[dom] == self.parent_grid_ratio[dom])
            assert (self.dy[dom-1]/self.dy[dom] == self.parent_grid_ratio[dom])
            assert self.i_parent_start[0] == 1
            assert self.j_parent_start[0] == 1


pbl_mapping = {
    0:  'none',
    1:  'YSU',
    2:  'MYJ',
    3:  'GFS',
    4:  'QNSE',
    5:  'MYNN2.5',
    6:  'MYNN3',
    7:  'ACM2',
    8:  'BouLac',
    9:  'UW',
    10: 'TEMF',
    11: 'Shin-Hong',
    12: 'GBM',
    99: 'MRF',
}
sfclay_mapping = {
    1: 'MOST', # w/ Carslon-Boland viscous sub-layer and standard similarity function lookup
    2: 'Eta', # Janjic Eta similarity, w/ Zilitinkevich thermal roughness length, standard similarity function lookup
    5: 'MYNN',
}
valid_sfclay = {
    # for each PBL scheme
    1:  [1],
    2:  [2],
    5:  [1,2,5,91],
    6:  [1,2,5,91],
    7:  [1,7,91],
    8:  [1,2,91],
    9:  [1,2,91],
    11: [1,91],
    12: [1,91],
    91: [1,91],
}


class Physics(object):
    """&physics namelist"""

    def __init__(self,nmldict):
        self.nml = nmldict
        self.parse_all()

    def __str__(self):
        s = 'WRF `physics` parameters\n'
        for dom in range(len(self.bl_pbl_physics)):
            s+=f'  d0{dom+1:d}: {self.bl_pbl_physics[dom]} PBL, {self.sf_sfclay_physics[dom]} surface layer\n'
        return s.rstrip()

    def parse_all(self):
        pbl_idx_list = self.nml['bl_pbl_physics']
        sfclay_idx_list = self.nml['sf_sfclay_physics']
        for pbl_idx,sfclay_idx in zip(pbl_idx_list, sfclay_idx_list):
            if sfclay_idx not in valid_sfclay[pbl_idx]:
                print(f'WARNING: unexpected pairing of bl_pbl_physics={pbl_idx} with sf_sfclay_idx={sfclay_idx}')
        self.bl_pbl_physics = [pbl_mapping.get(idx,'UNKNOWN') for idx in pbl_idx_list]
        self.sf_sfclay_physics = [sfclay_mapping.get(idx,'UNKNOWN') for idx in sfclay_idx_list]
        self.num_land_cat = self.nml['num_land_cat']


