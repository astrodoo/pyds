"""
filename:
    astroeq.py

PURPOSE:
    collection of formulae in astronomy

Routines:
    cs(gamma=5./3., **keywords): sound speed 
    vkep(pointmass,r=1.): Keplerian velocity in the gravitional potential due to point mass
    rsh(mbh): Schwarschild radius
    rbondi(mbh=1.,cs=1.): Bondi radius
    Mdotbondi(mbh=1.,rho=1.cs=1.): Bondi accretion rate
    Ledd(mbh): Eddington Luminosity
    Mdotedd(mbh,radeff=0.1): Eddington BH mass accretion rate
    Mjeans(T=0., rho=1., mmw=1.3): Jean's Mass
    Ljeans(T=0., rho=1., mmw=1.3): Jean's length
    tff(rho=1.): free-fall time scale

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

def vkep(pointmass,r=1.):
    """ Keplerian velocity in the gravitional potential due to point mass.
        returned value would be in cm/s

        args:
            pointmass: mass of the central point in g
        keywords:
            r: radius 
    """
    result = np.sqrt(unit.g*pointmass/r)
    return result

def rsh(mbh):
    """ Schwarschild radius
        mbh shoud be in cgs unit.
        returned value would be in cm unit.
    """
    result = 2.*unit.g*mbh / unit.c /unit.c
    return result

def rbondi(mbh=1.,cs=1.):
    """ Bondi radius
        mbh and cs(sound speed) should be in cgs unit
        returned value would be in cm unit.

        keywords:
            mbh: black hole mass in g unit
            cs: sound speed in cm/s unit
    """
    result = 2.*unit.g*mbh / cs/cs
    return result 

def Mdotbondi(mbh=1.,rho=1.,cs=1.):
    """ Bondi Accretion rate
        mbh, rho, cs(sound speed) should be in cgs unit
        returned value would be in cm unit.

        keywords:
            mbh: black hole mass in g unit
            rho: density in g/cm3 unit
            cs: sound speed in cm/s unit
    """
    result = 4.*np.pi * unit.g**2. * mbh**2. * rho / cs**3.
    return result 


def Ledd(mbh):
    """ Eddington Luminosity
        mbh should be in cgs unit.
        returned value would be in erg/s unit.
    """
    result = 4.*np.pi* unit.g * mbh * unit.mp*unit.c / unit.sigmaT
    return result

def Mdotedd(mbh,radeff=0.1):
    """ Eddington BH mass accretion rate 
        mbh should be in cgs unit.
        returned value would be in g/s unit

        args:
            mbh: black hole mass in g unit
        keywords:
            radeff: radiative efficiency (default=0.1)
    """
    result = 4.*np.pi* unit.g * mbh*unit.mp/radeff/unit.sigmaT/unit.c
    return result

def Mjeans(T=0., rho=1., mmw=1.3):
    """ Jean's mass (assume the uniform density at the spherical shape) in g unit
        (eq.(5.26) in astropedia)

        keywords:
           T: temperature in K
           rho: density in g/cm3
           mmw: mean molecular weight (default=1.3 for neutral solar abundance)
    """
    result = np.power(5.*unit.kb*T / (unit.g*mmw*unit.mh), 3./2.) \
            * np.sqrt(3./ (4.*np.pi*rho))
    return result

def Ljeans(T=0., rho=1., mmw=1.3):
    """ Jean's length (assume the uniform density at the spherical shape) in cm unit
        Ljeans = 2 x Rjeans, where Rjeans = (Mjeans/ (4/3 pi rho))^(1/3)
        (using eq.(5.27) for Rjeans in astropedia)

        keywords:
           T: temperature in K
           rho: density in g/cm3
           mmw: mean molecular weight (default=1.3 for neutral solar abundance)
    """
    result = np.sqrt(15.*unit.kb*T / (unit.g*mmw*unit.mh*np.pi*rho) )
    return result

def tff(rho=1.):
    """ free-fall time scale (assume the uniform density at the spherical shape) in s unit
        (eq.(5.28) in astropedia)

        keywords:
            rho: density in g/cm3
    """
    result = np.sqrt(3.*np.pi/32./unit.g/rho)
    return result
