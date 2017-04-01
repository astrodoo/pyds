#!/usr/bin/env python
# filename: 
#     img2data.py
# 
# PURPOSE:
#     Collect the data from plot in image files (png, jpg ...)
#     ### Note that the magic function should be the type that allows interactive plotting
#         such as '%matplotlib notebook'.
#
# Routines:
#     readimg: read the image data
#        args
#           imgplot: imgplot=ax.imshow(image)
#           ax     : axis
#           w      : (optional) widget for showing the position and value interactively
#        keywords  :
#           xlim/ylim: range of the xaxis/yaxis -> type: list array [min, max]
#           xlog/ylog: if given, the scale is considered as logarithmic scale --> type: boolean
#        Return variables
#           self.xdata/self.ydata: the data coordinate chosen by mouse clikc
#       
# 
# Written by:
#     Doosoo Yoon
#     Shanghai Astronomical Observatory
#   
# History:
#     Written, 30 March 2017
###############################################################
import ipywidgets as widgets
from IPython.display import display
import numpy as np

class readimg:
    def __init__(self, *args, **keywords):
        
        self.wflag = True
        if len(args) == 2:    
            self.imgplot = args[0]; self.ax = args[1]
            self.wflag = False
        elif len(args) == 3:
            self.imgplot = args[0]; self.ax = args[1]
            self.w = args[2]
        else:
            print 'Number of Args should be 2 or 3.'
            return None  
        
        if 'xlim' in keywords.keys():
            self.xlim = keywords['xlim']
        else:
            print 'xlim should be entered.'
            return None
        if 'ylim' in keywords.keys():
            self.ylim = keywords['ylim']
        else:
            print 'ylim should be entered.'
            return None
        
        if 'xlog' in keywords.keys():
            self.xlog = keywords['xlog']
        else:
            self.xlog = False
        if 'ylog' in keywords.keys():
            self.ylog = keywords['ylog']
        else:
            self.ylog = False
            
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
            self.ax.axhline(y=event.ydata)
            self.ax.axvline(x=event.xdata)
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
                self.w.value = 'xpixel=%f, ypixel=%f ,xdata=%f, ydata=%f'%(event.xdata, event.ydata, xinplt, yinplt)

            
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
