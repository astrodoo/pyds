"""
filename:
    astroeq.py

PURPOSE:
    collection of formulae in astronomy

Routines:
    cs(gamma=5./3., **keywords): sound speed 
    vkep(pointmass,r): Keplerian velocity in the gravitional potential due to point mass
    rsh(mbh): Schwarschild radius
    rBondi(mbh,cs): Bondi radius
    LEdd(mbh): Eddington Luminosity
    MdotEdd(mbh): Eddington BH mass accretion rate

Written by:
    Doosoo Yoon
    Shanghai Astronomical Observatory

History:
    Written, 22 November 2017
"""
from pyds.astro import astrounit as unit
import numpy as np
import sys

def cs(gamma=5./3., **keywords):
    """ sound speed
        Either of T (temperature) or P (pressure) & rho (density) should be entered
        keywords:
             gamma: adiabatic index (default: 5./3.)
        **keywords:
             T: temperature in K
             P: pressure in cgs
             rho: density in cgs
             mmw: mean molecular weight (default: 0.62 for fully ionized)
    """
    if 'T' in keywords.keys():
        T = keywords['T']

        if 'mmw' in keywords.keys():
            mmw = keywords['mmw']
        else:
            mmw = 0.62   # fully ionized

        result = np.sqrt(gamma*unit.kb*T / mmw/unit.mh)

    else:
        try:
            P = keywords['P']
            rho = keywords['rho']
        except:
            sys.exit('should input either of T or (P & rho)')

        result = np.sqrt(gamma* P / rho)
    return result

def vkep(pointmass,r):
    """ Keplerian velocity in the gravitional potential due to point mass
    """
    result = np.sqrt(unit.g*pointmass/r)
    return result

def rsh(mbh):
    """ Schwarschild radius
        mbh shoud be in cgs unit
    """
    result = 2.*unit.g*mbh / unit.c /unit.c
    return result

def rBondi(mbh,cs):
    """ Bondi radius
        mbh and cs(sound speed) should be in cgs unit
    """
    result = 2.*unit.g*mbh / cs/cs
    return result

def LEdd(mbh):
    """ Eddington Luminosity
        mbh should be in cgs unit
    """
    result = 4.*np.pi* unit.g * mbh * unit.mp*unit.c / unit.sigmaT
    return result

def MdotEdd(mbh,radeff=0.1):
    """ Eddington BH mass accretion rate 
        mbh should be in cgs unit
    """
    result = 4.*np.pi* unit.g * mbh*unit.mp/radeff/unit.sigmaT/unit.c
    return result
