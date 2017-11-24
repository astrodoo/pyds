from astropy import constants as cons
from astropy import units

consmatch = {'g'      : 'G',       \
             'lsun'   : 'L_sun',   \
             'msun'   : 'M_sun',   \
             'rsun'   : 'R_sun',   \
             'au'     : 'au',      \
             'c'      : 'c',       \
#             'e'      : 'e',       \
             'h'      : 'h',       \
             'kb'     : 'k_B',     \
             'pc'     : 'pc',      \
             'kpc'    : 'kpc',     \
             'me'     : 'm_e',     \
             'mn'     : 'm_n',     \
             'mp'     : 'm_p',     \
             'mh'     : 'm_p',     \
             'amu'    : 'u',       \
             'sigmaT' : 'sigma_T', \
             'sigmaSB': 'sigma_sb'}

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
    >>> astrounit().info()
    """
    for key,val in consmatch.iteritems():
        exec('%s = cons.%s.cgs.value'%(key,val))

    year = year
    lyr  = lyr
    eV   = eV
    Jy   = Jy
    re   = re

    def info(self):

        info()

        return None


def info():
    print '%8s \t %15s \t %s \t %s'%('name','name in astropy','value in cgs','unit')
    print '----------------------------------------------------------------------'
    for key,val in consmatch.iteritems():
        exec('astroval  = cons.%s.cgs.value'%val)
        exec('astrounit = cons.%s.cgs.unit'%val)
        print '%8s \t %15s \t %e \t %s'%(key,val,astroval,astrounit )

    print '----------------------------------------------------------------------'
    print 'not in astropy:'
    print '%8s \t %15s \t %e \t %s'%('year','not in astropy',year,'s' )
    print '%8s \t %15s \t %e \t %s'%('lyr','not in astropy',lyr,'cm' )
    print '%8s \t %15s \t %e \t %s'%('eV','not in astropy',eV,'erg' )
    print '%8s \t %15s \t %e \t %s'%('Jy','not in astropy',Jy,'erg / (cm2 s Hz)' )
    print '%8s \t %15s \t %e \t %s'%('re','not in astropy',re,'cm' )

