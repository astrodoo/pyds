"""
filename: 
     img2data.py
 
PURPOSE:
     Collect the data from plot in image files (png, jpg ...)
     ### Note that the magic function should be the type that allows interactive plotting
         such as '%matplotlib notebook'.

Routines:
     readimg: read the image data
        args
           imgplot: 2D image (e.g., imgplot=ax.imshow(image))
           ax     : axis of the plot
        keywords  :
           wflag  : if given, widget will show the information interactively (boolean; default= True)
           xlim/ylim: range of the xaxis/yaxis -> type: list array [min, max]
           xlog/ylog: if given, the scale is considered as logarithmic scale --> type: boolean
           noline : if given, line will not be shown after clicking the data point (booine; default= False)
        Return variables
           self.xdata/self.ydata: the data coordinate chosen by mouse clikc
 
Written by:
     Doosoo Yoon
     Shanghai Astronomical Observatory
   
History:
     Written, 30 March 2017
"""
from __future__ import print_function
import ipywidgets as widgets
from IPython.display import display
import numpy as np

class readimg:
    """
     readimg: read the image data
        args
           imgplot: 2D image (e.g., imgplot=ax.imshow(image))
           ax     : axis of the plot
        keywords  :
           wflag  : if given, widget will show the information interactively (boolean; default= True)
           xlim/ylim: range of the xaxis/yaxis corresponidng to the point you select with first four clicks-> ([min, max])
           xlog/ylog: if given, the scale is considered as logarithmic scale --> (boolean; default= False)
        Return variables
           self.xdata/self.ydata: the data coordinate chosen by mouse click
    """

    def __init__(self, *args, **keywords):
        
        self.wflag = True
        if len(args) == 2:    
            self.imgplot = args[0]; self.ax = args[1]
        else:
            print('Number of Args should be 2.')
            return None  

        self.w = widgets.HTML()
        
        if 'xlim' in keywords.keys():
            self.xlim = keywords['xlim']
        else:
            print('xlim should be entered.')
            return None
        if 'ylim' in keywords.keys():
            self.ylim = keywords['ylim']
        else:
            print('ylim should be entered.')
            return None
        
        if 'xlog' in keywords.keys():
            self.xlog = keywords['xlog']
        else:
            self.xlog = False
        if 'ylog' in keywords.keys():
            self.ylog = keywords['ylog']
        else:
            self.ylog = False
        if 'wflag' in keywords.keys():
            self.wflag = keywords['wflag']
        else:
            self.wflag = True
        if 'noline' in keywords.keys():
            self.noline = keywords['noline']
        else:
            self.noline = False
            
        self.x0 = []; self.y0 = []
        self.xx = []; self.yy = []
        self.xdata = []; self.ydata = []
        self.ind = 0
        self.cid = self.imgplot.figure.canvas.mpl_connect('button_press_event', self)
        if self.wflag: display(self.w)

    def __call__(self, event):
        if self.wflag:
            self.w.value = 'xdata=%f, ydata=%f'%(event.xdata, event.ydata)
        if self.ind < 2:
            self.ax.plot(event.xdata,event.ydata,'o')
            self.x0.append(event.xdata)
            if self.wflag:
                self.w.value = 'xdata=%f, ydata=%f'%(event.xdata, event.ydata)
        elif self.ind < 4:
            self.ax.plot(event.xdata,event.ydata,'o')
            self.y0.append(event.ydata)
            if self.wflag:
                self.w.value = 'xdata=%f, ydata=%f'%(event.xdata, event.ydata)
        else:
            if not(self.noline):
                self.ax.axhline(y=event.ydata, linewidth=1)
                self.ax.axvline(x=event.xdata, linewidth=1)

            self.xx.append(event.xdata); self.yy.append(event.ydata)
            
            if self.xlog:
                xinplt = self.calcdataxlg(event.xdata)
            else:
                xinplt = self.calcdatax(event.xdata)

            self.xdata.append(xinplt)
            
            if self.ylog:
                yinplt = self.calcdataylg(event.ydata)
            else:
                yinplt = self.calcdatay(event.ydata)

            self.ydata.append(yinplt)
            
            if self.wflag:
                self.w.value = 'xpixel=%f, ypixel=%f ,xdata=%e, ydata=%e'%(event.xdata, event.ydata, xinplt, yinplt)

        self.ind = self.ind + 1
        
    # Calculate Data position in the plot coordinate
    def calcdatax(self, xpos):
        xdat = (xpos-self.x0[0])/(self.x0[1]-self.x0[0])*(self.xlim[1]-self.xlim[0])+self.xlim[0]
        return xdat
    
    def calcdatay(self, ypos):
        ydat = (ypos-self.y0[0])/(self.y0[1]-self.y0[0])*(self.ylim[1]-self.ylim[0])+self.ylim[0]
        return ydat
    
    def calcdataxlg(self, xpos):
        xlimlg = np.log10(self.xlim)
        xdat = np.power(10.,(xpos-self.x0[0])/(self.x0[1]-self.x0[0])*(xlimlg[1]-xlimlg[0])+xlimlg[0])
        return xdat
    
    def calcdataylg(self, ypos):
        ylimlg = np.log10(self.ylim)
        ydat = np.power(10.,(ypos-self.y0[0])/(self.y0[1]-self.y0[0])*(ylimlg[1]-ylimlg[0])+ylimlg[0])
        return ydat
