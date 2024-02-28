import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class SurfaceHistory(object):
    """Process text diagnostic output that was written out with:
        ```
        erf.v           = 1
        erf.data_log    = hist.txt
        ```
    """
    timename = 't' # 'time'
    heightname = 'z' # 'height'
    surfvars = ['t','ustar','tstar','L']

    def __init__(self, histfile, t0=0.0, dt=None):
        """Load diagnostic profile data from 3 datafiles

        Parameters
        ----------
        *args : list or glob string
            Averaged profile datafiles to load
        t0 : float, optional
            With `dt`, used to overwrite the timestamps in the text data
            to address issues with insufficient precision
        dt : float, optional
            Overwrite the time dimension coordinate with
            t0 + np.arange(Ntimes)*dt
        """
        self.df = pd.read_csv(histfile, delim_whitespace=True,
                              names=self.surfvars)
        self.df = self.df.drop_duplicates().set_index('t')
        if dt is not None:
            dt0 = surf['time'].iloc[0]
            if not dt==dt0:
                print('Warning: inconsistent specified dt and first step',dt,dt0)
            self.df = self.df.reset_index()
            self.df['t'] = t0 + np.arange(1,len(self.df)+1)*dt
            self.df = self.df.set_index('t')

    def plot(self,*args):
        """Quick plotting function"""
        varns = self.surfvars[1:] if (len(args) == 0) else args
        if isinstance(varns,str):
            varns = [varns]
        fig,axs = plt.subplots(nrows=len(varns),figsize=(5,2.5*len(varns)))
        if len(varns) == 1:
            axs = [axs]
        for varn,ax in zip(varns,axs):
            ax.plot(self.df.index, self.df[varn])
            ax.set_ylabel(varn)
            ax.grid()
        axs[-1].set_xlabel('simulation time [s]')
        return fig,axs
