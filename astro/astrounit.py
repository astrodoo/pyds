from __future__ import print_function
from future.utils import iteritems
from astropy import constants as cons
from astropy import units

consmatch = {'g'      : 'G',       \
             'lsun'   : 'L_sun',   \
             'msun'   : 'M_sun',   \
             'rsun'   : 'R_sun',   \
             'au'     : 'au',      \
             'c'      : 'c',       \
#             'e'      : 'e',       \
             'alpha'  : 'alpha',   \
             'h'      : 'h',       \
             'kb'     : 'k_B',     \
             'pc'     : 'pc',      \
             'kpc'    : 'kpc',     \
             'me'     : 'm_e',     \
             'mn'     : 'm_n',     \
             'mp'     : 'm_p',     \
#             'mh'     : 'm_p',     \
             'amu'    : 'u',       \
             'sigmaT' : 'sigma_T', \
             'sigmaSB': 'sigma_sb'}

mh = cons.m_p.cgs.value + cons.m_e.cgs.value
esu = cons.e.esu.value
year = 60.*60.*24.*365.
lyr = cons.c.cgs.value *year
eV  = units.eV.to(units.erg)
Jy = (1*units.Jy).cgs.value
re  = 2.81794092e-13

class astrounit: 
    """
    assign the astrophysical constants.
    The available constants will be shown as follows:
    >>> from pyds.astrounit import astrounit
    >>> astrounit.info()
    >>> astrounit.g
    """
    for (key,val) in iteritems(consmatch):
        exec('%s = cons.%s.cgs.value'%(key,val))

    mh   = mh
    esu  = esu
    year = year
    lyr  = lyr
    eV   = eV
    Jy   = Jy
    re   = re

    @staticmethod
    def info():

        info()

        pass


def info():
    print('%8s \t %15s \t %s \t %s'%('name','name in astropy','value in cgs','unit'))
    print('----------------------------------------------------------------------')
    for (key,val) in iteritems(consmatch):
        execscope = {}           # dictionary for exec scope     
        exec('astroval  = cons.%s.cgs.value'%val, globals(), execscope)
        astroval = execscope['astroval']
        exec('astrounit = cons.%s.cgs.unit'%val, globals(), execscope)
        astrounit = execscope['astrounit']
        print('%8s \t %15s \t %e \t %s'%(key,val,astroval,astrounit ))

    print('%8s \t %15s \t %e \t %s'%('esu','e.esu',esu,'cgs' ))
    print('----------------------------------------------------------------------')
    print('not in astropy:')
    print('%8s \t %15s \t %e \t %s'%('mh','not in astropy',mh,'g' ) )
    print('%8s \t %15s \t %e \t %s'%('year','not in astropy',year,'s' ) )
    print('%8s \t %15s \t %e \t %s'%('lyr','not in astropy',lyr,'cm' ) )
    print('%8s \t %15s \t %e \t %s'%('eV','not in astropy',eV,'erg' ) )
    print('%8s \t  %15s \t %e \t %s'%('Jy','not in astropy',Jy,'erg / (cm2 s Hz)' ) )
    print('%8s \t %15s \t %e \t %s'%('re','not in astropy',re,'cm' ) )

