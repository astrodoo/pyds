"""
filename: 
     mimicAlpha_to_rgb.py
 
PURPOSE:
     get the tuple of rgb, mimicking the rgb with alpha value

DEPENDENCY:
     matplotlib.colors
     numpy

Written by:
     Doosoo Yoon
     Shanghai Astronomical Observatory
   
History:
     Written, 18 January 2018
"""
import matplotlib.colors as cls
import numpy as np

def mimicAlpha_to_rgb(color,alpha=0.5,bg='w'):
    """
    get the tuple of rgb, mimicking the rgb with alpha value

    args:
        color: color name or color tuple (r,g,b)

    keywords:
        alpha: alpha value, representing the transparency (0 ~ 1)
        bg: back ground color name or tuple (r,g,b)
    return:
        tuple of returned (r,g,b) value
    """

    if type(color).__name__ == 'str':
        clr_rgb = cls.to_rgb(color)
    elif ((type(color).__name__ != 'tuple') | (len(color) != 3)):
        raise ValueError('color should be in string or (r,g,b) tuple')
    else:
        clr_rgb = color

    if type(bg).__name__ == 'str':
        clr_bg  = cls.to_rgb(bg)
    elif ((type(bg).__name__ != 'tuple') | (len(bg) != 3)):
        raise ValueError('bg should be in string or (r,g,b) tuple')
    else:
        clr_bg = bg
        
    clr_rgb = np.asarray(clr_rgb,dtype='float')
    clr_bg  = np.asarray(clr_bg,dtype='float')
    
    clrA_rgb = (1.-alpha)*clr_bg + alpha*clr_rgb

    return tuple(clrA_rgb)
