"""
filename: 
     mtickexp.py
 
PURPOSE:
     transform the tick format to exponential math type (e.g. 3 x 10^{10})

Written by:
     Doosoo Yoon 
     Shanghai Astronomical Observatory
   
History:
     Written, 26 January 2018
"""
import matplotlib.ticker as mtick

def mtickexp(ax,axis='x',pointnum=1):
    """
    transform the tick format to exponential math type (e.g. 3 x 10^{10})

    args:
        ax: axis class

    keywords:
        axis: indicator of x/y axis (string, default: 'x')
        pointnum: the number of digits below the point (integer, default:1)
    """


    fmt = r'%.'+str(pointnum)+'fx$10^{%i}$'

    def format_tick(x,pos=None):
	xstr_e = '%e'%x
    	xsplt = xstr_e.split('e')
    	xfact = float(xsplt[0])
    	xexp  = int(xsplt[1])
    
        if x != 0:
            xstr = fmt%(xfact,xexp)
        else:
            xstr = '0'
        return xstr 

    if (axis=='x'):
        ax.xaxis.set_major_formatter(mtick.FuncFormatter(format_tick))
    else:
        ax.yaxis.set_major_formatter(mtick.FuncFormatter(format_tick))
