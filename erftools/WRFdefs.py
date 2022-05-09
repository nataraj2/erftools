from datetime import datetime

class time_control(object):
    def __init__(self,nmldict):
        self.legacy = nmldict
        self.parse_datetime_range()

    def __str__(self):
        return f"""WRF `time_control` parameters
  sim range: {self.start_datetime} to {self.end_datetime}"""

    def parse_datetime_range(self):
        # assume all domains have the same start/end datetime for now; TODO: allow for variable start times per domain
        idom = 0
        start_year  = self.legacy['start_year'][idom]
        start_month = self.legacy['start_month'][idom]
        start_day   = self.legacy['start_day'][idom]
        start_hour  = self.legacy['start_hour'][idom]
        end_year  = self.legacy['end_year'][idom]
        end_month = self.legacy['end_month'][idom]
        end_day   = self.legacy['end_day'][idom]
        end_hour  = self.legacy['end_hour'][idom]
        self.start_datetime = datetime(start_year, start_month, start_day, start_hour)
        self.end_datetime = datetime(end_year, end_month, end_day, end_hour)
        
    
class domains(object):
    def __init__(self,nmldict):
        self.legacy = nmldict
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
        self.max_dom = self.legacy['max_dom']
        # UNstaggered grid points == cells
        # - assume s_we == s_sn == s_vert == [1, 1, ...]
        self.e_we = self.legacy['e_we'][:self.max_dom] # west--east
        self.e_sn = self.legacy['e_sn'][:self.max_dom] # south--north
        self.e_vert = self.legacy['e_vert'][:self.max_dom] # bottom--top
        self.dx = self.legacy['dx'][:self.max_dom]
        self.dy = self.legacy['dy'][:self.max_dom]
        self.i_parent_start = self.legacy['i_parent_start'][:self.max_dom]
        self.j_parent_start = self.legacy['j_parent_start'][:self.max_dom]
        parent_grid_ratio = self.legacy['parent_grid_ratio'][:self.max_dom]
        for dom in range(1,self.max_dom):
            assert (self.dx[dom-1]/self.dx[dom] == parent_grid_ratio[dom])
            assert (self.dy[dom-1]/self.dy[dom] == parent_grid_ratio[dom])

    def parse_time_integration(self):
        self.time_step = self.legacy['time_step'] # seconds
        self.parent_time_step_ratio = self.legacy['parent_time_step_ratio']


