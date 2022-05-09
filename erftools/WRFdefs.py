from datetime import datetime

class time_control(object):
    def __init__(self,nmldict):
        self.legacy = nmldict
        self.parse_datetime_range()

    def __str__(self):
        return f"""WRF `time_control` parameters
  sim range: {self.start_datetime} to {self.end_datetime} ({self.stop_time} s)"""

    def parse_datetime_range(self):
        # assume all domains have the same start/end datetime for now
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
        self.stop_time = (self.end_datetime - self.start_datetime).total_seconds()
        
    
