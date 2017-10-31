from astropy import constants as cons

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
             'me'     : 'm_e',     \
             'mn'     : 'm_n',     \
             'mp'     : 'm_p',     \
             'mh'     : 'm_p',     \
             'amu'    : 'u',       \
             'sigmaT' : 'sigma_T', \
             'sigmaSB': 'sigma_sb'}

lyr = cons.c.cgs.value *60.*60.*24.*365

class astrounit:

    for key,val in consmatch.iteritems():
        exec('%s = cons.%s.cgs.value'%(key,val))

    lyr = lyr


def info():
    print '%8s \t %15s \t %s \t %s'%('name','name in astropy','value in cgs','unit')
    print '----------------------------------------------------------------------'
    for key,val in consmatch.iteritems():
        exec('astroval  = cons.%s.cgs.value'%val)
        exec('astrounit = cons.%s.cgs.unit'%val)
        print '%8s \t %15s \t %e \t %s'%(key,val,astroval,astrounit )

    print '----------------------------------------------------------------------'
    print 'not in astropy:'
    print '%8s \t %15s \t %e \t %s'%('lyr','not in astropy',lyr,'cm' )


