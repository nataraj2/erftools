import numpy as np
from .constants import Gamma, p_0, R_d, R_v, Cp_d, Cp_v

def getdPdRgivenConstantTheta(rho, theta):
    return Gamma * p_0 * (R_d * theta / p_0)**Gamma * rho**(Gamma-1.0)

def getPgivenRTh(rhotheta, qv=0, precip=True):
    """Function to return pressure given density times theta

    Parameters
    ----------
    rhotheta: float
        Density times potential temperature
    qv: float, optional
        Water vapor mixing ratio
    precip: bool
        Precipitation is modeled by a single-moment microphysics model,
        otherwise a warm moisture model is used.
    """
    if qv==0:
        return p_0 * (R_d * rhotheta / p_0)**Gamma
    elif not precip:
        rhotheta_t = rhotheta*(1.0+qv);
        return p_0 * (R_d * rhotheta_t / p_0)**Gamma
    else:
        R_t        =  R_d + qv*R_v;
        Cp_t       = Cp_d + qv*Cp_v;
        Gamma_t    = Cp_t / (Cp_t-R_t);
        rhotheta_t = rhotheta*(1.0+qv);
        return p_0 * (R_t * rhotheta_t / p_0)**Gamma_t

def getRhoThetagivenP(p):
    return (p * p_0**(Gamma-1))**(1.0/Gamma) / R_d

def getRhogivenThetaPress(theta, p, rdOcp=R_d/Cp_d, qv=0.0):
    return p_0**rdOcp * p**(1.0/Gamma) / (R_d * theta * (1.0 + R_v/R_d*qv))
