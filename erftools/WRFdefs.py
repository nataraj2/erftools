from datetime import datetime

class time_control(object):
    def __init__(self,nmldict):
        self.nml = nmldict
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
        
    
class domains(object):
    def __init__(self,nmldict):
        self.nml = nmldict
        self.parse_grid()
        self.parse_time_integration()

    def __str__(self):
        s = 'WRF `domains` parameters\n'
        for dom in range(self.max_dom):
            s += f'  domain d{dom+1:02d}:' \
                 f' [(1,{self.e_we[dom]:d}),(1,{self.e_sn[dom]:d}),(1,{self.e_vert[dom]:d})]' \
                 f' ds=[{self.dx[dom]},{self.dy[dom]}]\n'
        return s

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
        parent_grid_ratio = self.nml['parent_grid_ratio'][:self.max_dom]
        for dom in range(1,self.max_dom):
            assert (self.dx[dom-1]/self.dx[dom] == parent_grid_ratio[dom])
            assert (self.dy[dom-1]/self.dy[dom] == parent_grid_ratio[dom])

    def parse_time_integration(self):
        self.time_step = self.nml['time_step'] # seconds
        self.parent_time_step_ratio = self.nml['parent_time_step_ratio']


