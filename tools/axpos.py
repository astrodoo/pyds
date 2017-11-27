"""
filename: 
     axpos.py
 
PURPOSE:
     set the position of multiplots

Routines:
     axpos: set the configuration of multiplots

Written by:
     Doosoo Yoon 
     Shanghai Astronomical Observatory
   
History:
     Written, 26 November 2017
"""
import numpy as np
import sys

class axpos:
    def __init__(self,multplt,pltx0=[50.,10.],plty0=[30.,10.],pltxs=200.,pltys=100.,pltxw=0.,pltyw=5.):
	""" 
	set the position of multiplots

        args:
            multplt: two digits integer, indicating the number of row & column 
                    (ex, 23 -> 2 rows and 3 columns)
                    
        keywords:
            pltx0/plty0: two elements array for [left margin, right margin]
            pltxs/pltys: (x/y) size of the plots
                         If single number is given, all plots have same (x/y) size.
                         If list of numbers is given, the length of the list should be same with
                         (x/y) size.
            pltxw/pltyw: gap between arrays

        return:
            self.winxs/self.winys: total (x,y) size of the figure.
                                   These are useful in setting the ratio of the figure. 
                                   (ex.  >>> plt.figure(figsize=(8.,8.*self.winys/self.winxs)))

	"""
        multplt = str(multplt)
        
        ny = int(multplt[0])
        nx = int(multplt[1])
        
        if (type(pltxs).__name__=='list'):
            if (len(pltxs) != nx):
                sys.exit('length of pltxs should be equal to nx.')
            self.pltxs = map(float,pltxs)
        else:
            self.pltxs = [float(pltxs) for i in range(nx)]
      
        if (type(pltys).__name__=='list'):
            if (len(pltys) != ny):
                sys.exit('length of pltys should be equal to ny.')
            self.pltys = map(float,pltys)
        else:
            self.pltys = [float(pltys) for i in range(ny)]
        
        self.winxs=np.sum(pltx0)+np.sum(self.pltxs)+(nx-1)*float(pltxw)
        self.winys=np.sum(plty0)+np.sum(self.pltys)+(ny-1)*float(pltyw)
    
        self.pltx0 = map(float,pltx0); self.plty0 = map(float,plty0)
        self.pltxw = float(pltxw); self.pltyw = float(pltyw)
        
            
    def pos(self,numplt):
        """
        set the position of multiplots with given identifier

        args:
            numplt: three digits integer, indicating the number of (row/column/plot identifier)
                    --> same rule with matplotlib subplots

        return:
            list of [left,bottom,width,height]
            --> This can be used as
                >>> ax1 = plt.axes(self.plot(322))
        """
        
        numpltstr = str(numplt)
        
        row    = int(numpltstr[0])
        col    = int(numpltstr[1])
        posplt = int(numpltstr[2])
        
        posrow = row - 1 - int((posplt-1) // col)
        poscol = int(np.mod((posplt-1),col))
        
        if (row*col < posplt):
            sys.exit('plot number is out of range.')
        
        left = (self.pltx0[0]+self.pltxw*poscol+np.sum(self.pltxs[:poscol]))/self.winxs
        width = (self.pltxs[poscol])/self.winxs 
        bottom = (self.plty0[0]+self.pltyw*posrow+np.sum(self.pltys[:posrow]))/self.winys
        height = (self.pltys[posrow])/self.winys
        
        return [left,bottom,width,height]
