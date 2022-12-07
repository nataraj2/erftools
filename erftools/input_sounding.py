import numpy as np
import matplotlib.pyplot as plt

class InputSounding(object):

    def __init__(self,fpath='input_sounding'):
        with open(fpath,'r') as f:
            # The first line includes the surface pressure (hPa), potential
            # temperature (K) and moisture mixing ratio (g/kg).
            p_surf, T_surf, qv_surf = [float(val) for val in f.readline().split()]

            # Each subsequent line has five input values: height (meters above
            # sea-level), dry potential temperature (K), vapor mixing ratio
            # (g/kg), x-direction wind component (m/s), and y-direction wind
            # component (m/s)
            init = np.loadtxt(f)

        self.p_surf = p_surf * 100 # [Pa]
        self.T_surf = T_surf # [K]
        self.qv_surf = qv_surf * 0.001 # [g/g]

        self.z  = init[:,0]
        self.T  = init[:,1]
        self.qv = init[:,2] * 0.001
        self.u  = init[:,3]
        self.v  = init[:,4]

    def plot(self):
        fig,ax = plt.subplots(ncols=3,sharey=True,figsize=(10,8))
        ax[0].plot(self.T, self.z)
        ax[0].plot(self.T_surf, 0, 'ko', markerfacecolor='none')
        ax[1].plot(self.qv, self.z)
        ax[1].plot(self.qv_surf, 0, 'ko', markerfacecolor='none')
        ax[2].plot(self.u, self.z, label='u')
        ax[2].plot(self.v, self.z, label='v')
        ax[0].set_ylabel('altitude AGL [m]')
        ax[0].set_xlabel('dry potential temperature [K]')
        ax[1].set_xlabel('water vapor mixing ratio [g/g]')
        ax[2].set_xlabel('horizontal wind component [m/s]')
        ax[2].legend(loc='best')
        for axi in ax:
            axi.grid()
        return fig,ax

