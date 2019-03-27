"""
filename: 
     cosmonom.py
 
PURPOSE:
     Calculate cosmological variables under the given parameters:
     Hubble constant(H0), Omega_Lambda, Omega_m

     All equations and coding were referred to Paper-and-pencil cosmological calculator in Pilipenko (2013):
     http://adsabs.harvard.edu/abs/2013arXiv1303.5961P

Written by:
     Doosoo Yoon
     Shanghai Astronomical Observatory
   
History:
     Written, 2 April 2018
"""
from __future__ import print_function
import numpy as np

def solve(f,val,a0,b0,eps):
    a = a0
    b = b0
    epsilon = b - a
    df = f(b) - f(a)
    while (epsilon > eps):
        c = (a+b)/2.
        fc = f(c)
        if ((((fc-val)>=0)&(df>0)) | (((fc-val)<0)&(df<0))):
            b = c
        else:
            a = c
            
        epsilon = b - a
    solve = c
    return solve


def rombint(f,a,b,tol):
    MAXITER = 30; MAXJ = 5
    
    g = np.zeros(MAXJ+1, dtype=np.float64)
    h = 0.5*(b-a)
    gmax = h*( f(a) + f(b) )
    g[0] = gmax
    nint = 1
    error = 1e20
    i = 0
    
    while (i<=MAXITER):
        if ((i>5) & (np.abs(error)<tol)): 
            break
# Calculate next trapezoidal rule approximation to integral.
        g0 = 0.
    
        for k in np.arange(nint)+1:            
            g0 = g0 + f(a+(k+k-1)*h)
            
        g0 = 0.5 * g[0] + h*g0                                                      
        h  = 0.5 * h
        nint = nint+nint
        jmax = np.min([i,MAXJ])  
        fourj=1.0
        
        for j in range(jmax):
# Use Richardson extrapolation.
            fourj = 4.*fourj
            g1    = g0+(g0-g[j])/(fourj-1.)
            g[j]  = g0
            g0    = g1
        if (np.abs(g0)>tol):
            error = 1. - gmax/g0
        else:
            error = gmax

        gmax      = g0
        g[jmax] = g0
        i         = i+1

    rombint = g0
    if ((i>MAXITER) & (np.abs(error)>tol)):
        print('Rombint failed to converge; integral, error=%e,%e'%(rombint,error) )
    
    return rombint


class cosmonom:
    """
    Calculate cosmological variables under the given parameters:
    Hubble constant(H0), Omega_Lambda, Omega_m

    All equations and coding were referred to Paper-and-pencil cosmological calculator in Pilipenko (2013):
    http://adsabs.harvard.edu/abs/2013arXiv1303.5961P

    keywords -- H0: Hubble constants in km/s/Mpc (default= 67.15)
	        OL: Omega Lambda (default= 0.683)
                Om: Omega M (default=0.317)
                zmin: minimum value of red-shift (default=0.)
                zmax: maximum value of red-shift (default=20.)

    outputs   -- self.zext: extremum value of z
                 self.age0: self.age_z(0.)   in Gyr
    functions -- self.Hz: H(z) in km/s/Mpc    <->  self.zH: z(H)   
                 self.rz: r(z) in Mpc   <->  self.zr: z(r)   - comoving radius  
                 self.dmz: dm(z)  <->  self.zdm: z(dm)
                 self.size_z: 1" size(z) in kpc <->  self.zSize1, self.zSize2 - z(size) for z<zext, z>zext, respectively.   
                 self.angle_z: 1 kpc angle(z) in arcsec <-> self.zAngle1, self.zAngle2 - z(angle) for z<zext, z>zext, respectively.
                 self.age_z: age(z)  in Gyr <-> self.zage: z(age)
                 self.zt : self.zage(self.age0-t) - z for the corresponding lookback time in Gyr

    """
    def __init__(self,H0=67.15,OL=0.683,Om=0.317,zmin=0.,zmax=20.):
        self.H0 = H0
        self.OL = OL
        self.Om = Om

        self.zmin = zmin
        self.zmax = zmax
        if self.zmin==0.:
            self.zmin = np.min([self.zmax/10.,0.1])
	
        self.Hz = lambda z: self.H0 * np.sqrt( self.OL + (1.+z)**3 * self.Om)
        self.zH = lambda H: solve(self.Hz, H, self.zmin-1., self.zmax*2., 1e-4)
        rf_f         = lambda x: 1./np.sqrt((1.+x)**3. + self.OL/self.Om)
        self.rz      = lambda z: rombint(rf_f, 0., z, 0.1) * 299792.0 / self.H0 / np.sqrt(self.Om)
        self.zr      = lambda r: solve(self.rz, r, self.zmin, self.zmax*2, 1e-4)
        self.dmz     = lambda z: 5. * np.log10(self.rz(z) * (1.+z)) + 25.
        self.zdm     = lambda dm: solve(self.dmz, dm, self.zmin, self.zmax*2., 1e-4)

        eq_zext      = lambda z: self.rz(z)/(1.+z) - rf_f(z) * 299792. / self.H0 / np.sqrt(self.Om)
        self.zext    = solve(eq_zext, 0., self.zmin, self.zmax*2, 1e-4)
        
        self.size_z  = lambda z: self.rz(z)/(1.+z) * 1./206265. * 1e3  # in kpc
        self.zSize1  = lambda s: solve(self.size_z, s, self.zmin, self.zext, 1e-4)
        self.zSize2  = lambda s: solve(self.size_z, s, self.zext, self.zmax*2, 1e-4)
        self.angle_z = lambda z: 1. / self.size_z(z)
        self.zAngle1 = lambda a: solve(self.angle_z, a, self.zmin+1e-4, self.zext, 1e-4)
        self.zAngle2 = lambda a: solve(self.angle_z, a, self.zext, self.zmax*2, 1e-4)
        
        agef         = lambda z: 1./(1.+z) / np.sqrt(self.OL + self.Om*(1.+z)**3)
        self.age_z   = lambda z: rombint(agef, z, 1000.0, 0.001) * 977.8 / self.H0
        self.zage    = lambda age: solve(self.age_z, age, self.zmin-0.5, self.zmax*2, 1e-4)
        self.age0    = self.age_z(0.)
        self.zt      = lambda t: self.zage(self.age0 - t)
        
        self.Hz.__doc__      = "function: H(z) in km/s/Mpc"
        self.zH.__doc__      = "function: z(H)"
        self.rz.__doc__      = "function: comoving radius r(z) in Mpc" 
        self.zr.__doc__      = "function: z(r) where r is comving radius in unit of kpc" 
        self.dmz.__doc__     = "function: dm(z)"
        self.zdm.__doc__     = "function: z(dm)"
        self.size_z.__doc__  = 'function: 1" size(z) in kpc'
        self.zSize1.__doc__  = "function: z(size) where z < z_ext"
        self.zSize2.__doc__  = "function: z(size) where z > z_ext"
        self.angle_z.__doc__ = "function: 1kpc angle(z) in arcsec"
        self.zAngle1.__doc__ = "function: z(angle) where z < z_ext"
        self.zAngle2.__doc__ = "function: z(angle) where z > z_ext"
        self.age_z.__doc__   = "function: age(z) in Gyr"
        self.zage.__doc__    = "function: z(age)"
        self.zt.__doc__      = "function: zage(age0 - t)"


    def draw(self, **keywords):
        """
        Draw the scaled axes for the cosmological variables

        **keywords -- out: if given, the plot will be saved to the designated output file (type = string) 
                        z: if given, it will draw the horizontal line for the value of z so that it 
                           helps catching the corresponding variables. (type = float)
        """
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker

        print('z_extremum: %f'%self.zext)
        print('size(zext): %f, angle(zext): %f'%(self.size_z(self.zext),self.angle_z(self.zext)))

        ParamStr = r'$H_{0}$=%5.2f, $\Omega_{\Lambda}$=%5.3f, $\Omega_{M}$=%5.3f'%(self.H0,self.OL,self.Om)

        # Setup a plot such that only the bottom spine is shown
        def setup(ax):
    	    ax.spines['top'].set_color('none')
    	    ax.spines['bottom'].set_color('none')
    	    ax.xaxis.set_major_locator(ticker.NullLocator())
    	    ax.spines['left'].set_color('none')
    	    ax.yaxis.set_ticks_position('right')
    	    ax.tick_params(axis='y',which='major',length=10,direction='inout')
    	    #ax.tick_params(axis='y',which='minor',right='False')

        #figure setting
        x0 = [30.,50.]; y0 = [20.,70.]
        ys = 500.
        xw = 50.; xs = xw/10.
        naxes = 9
        winxs = np.sum(x0)+(naxes-1)*xw
        winys = np.sum(y0) + ys
        
        x0 = x0/winxs; y0 = y0/winys
        xs = xs/winxs; ys = ys/winys
        xw = xw/winxs

        fig = plt.figure(figsize=(15.*winxs/winys,15))
        
        # draw the 1st axis (z)
        ax1 = fig.add_axes([x0[0],y0[0],xs,ys])
        setup(ax1)
        ax1.set_ylim(self.zmax,self.zmin)
        ax1.set_yscale('log')
        ticklocz = np.array([0.2,0.4,0.6,0.8,1.,2.,3.,4.,5.,6.,7.,8.,10.,12.,14.,16.,18.,20.])
        ticklabz_str = ["%.1f"%z for z in ticklocz]
        ax1.yaxis.set_major_locator(ticker.FixedLocator(ticklocz))
        ax1.yaxis.set_major_formatter(ticker.FixedFormatter(ticklabz_str))
        ax1.annotate('z',(1.,1.02),xycoords='axes fraction',ha='center')
        
        ax1.annotate(ParamStr,(1.,1.07),xycoords='axes fraction',ha='left',fontsize='large')
        
        # draw the 2nd axis (H)
        ax2 = fig.add_axes([x0[0]+xw,y0[0],xs,ys])
        setup(ax2)
        ax2.set_ylim(ax1.get_ylim())
        ax2.set_yscale(ax1.get_yscale())
        ticklab = np.array([80.,100.,200.,300.,400.,600.,800.,1000.,2000.,3000.])
        tickloc = [self.zH(H) for H in ticklab]
        ticklab_str = ["%.0f"%H for H in ticklab]
        ax2.yaxis.set_major_locator(ticker.FixedLocator(tickloc))
        ax2.yaxis.set_major_formatter(ticker.FixedFormatter(ticklab_str))
        ax2.annotate('H',(1.,1.02),xycoords='axes fraction',ha='center')
        
        # draw the 3rd axis (r_comov)
        ax3 = fig.add_axes([x0[0]+2*xw,y0[0],xs,ys])
        setup(ax3)
        ax3.set_ylim(ax1.get_ylim())
        ax3.set_yscale(ax1.get_yscale())
        ticklab = np.array([i for i in np.arange(11)*1000.])
        tickloc = [self.zr(r) for r in ticklab]
        ticklab_str = ["%.0f"%r for r in ticklab]
        ax3.yaxis.set_major_locator(ticker.FixedLocator(tickloc))
        ax3.yaxis.set_major_formatter(ticker.FixedFormatter(ticklab_str))
        ax3.annotate('r_comov',(1.,1.02),xycoords='axes fraction',ha='center')
        
        # draw the 4th axis (dm)
        ax4 = fig.add_axes([x0[0]+3*xw,y0[0],xs,ys])
        setup(ax4)
        ax4.set_ylim(ax1.get_ylim())
        ax4.set_yscale(ax1.get_yscale())
        ticklab = np.array([35.,40.,42.,44.,46.,48.,50.,52.])
        tickloc = [self.zdm(dm) for dm in ticklab]
        ticklab_str = ["%.0f"%dm for dm in ticklab]
        ax4.yaxis.set_major_locator(ticker.FixedLocator(tickloc))
        ax4.yaxis.set_major_formatter(ticker.FixedFormatter(ticklab_str))
        ax4.annotate('dm',(1.,1.02),xycoords='axes fraction',ha='center')
        
        # draw the 5th axis (age)
        ax5 = fig.add_axes([x0[0]+4*xw,y0[0],xs,ys])
        setup(ax5)
        ax5.set_ylim(ax1.get_ylim())
        ax5.set_yscale(ax1.get_yscale())
        ticklab = np.array([13.,12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1.,0.8,0.6,0.4,0.2,0.1])
        tickloc = [self.zage(age) for age in ticklab]
        ticklab_str = ["%.1f"%age for age in ticklab]
        ax5.yaxis.set_major_locator(ticker.FixedLocator(tickloc))
        ax5.yaxis.set_major_formatter(ticker.FixedFormatter(ticklab_str))
        ax5.annotate('age',(1.,1.02),xycoords='axes fraction',ha='center')
        
        # draw the 6th axis (time)
        ax6 = fig.add_axes([x0[0]+5*xw,y0[0],xs,ys])
        setup(ax6)
        ax6.set_ylim(ax1.get_ylim())
        ax6.set_yscale(ax1.get_yscale())
        ticklab = np.array([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,12.2,12.4,12.6,12.8,13.,13.2,13.4,13.6,13.8])
        tickloc = [self.zt(t) for t in ticklab]
        ticklab_str = ["%.1f"%t for t in ticklab]
        ax6.yaxis.set_major_locator(ticker.FixedLocator(tickloc))
        ax6.yaxis.set_major_formatter(ticker.FixedFormatter(ticklab_str))
        ax6.annotate('time',(1.,1.02),xycoords='axes fraction',ha='center')
        
        # draw the 7th axis (size 1")
        ax7 = fig.add_axes([x0[0]+6*xw,y0[0],xs,ys])
        setup(ax7)
        ax7.set_ylim(ax1.get_ylim())
        ax7.set_yscale(ax1.get_yscale())
        ticklab1 = np.array([0.,2.,4.,6.,8.,8.7])
        ticklab1 = ticklab1[ticklab1<self.size_z(self.zext)]
        tickloc1 = [self.zSize1(s) for s in ticklab1]
        ticklab2 = ticklab1[::-1]
        tickloc2 = [self.zSize2(s) for s in ticklab2]
        ticklab = np.concatenate((ticklab1,ticklab2))
        tickloc = np.concatenate((tickloc1,tickloc2))
        ticklab_str = ["%.1f"%t for t in ticklab]
        ax7.yaxis.set_major_locator(ticker.FixedLocator(tickloc))
        ax7.yaxis.set_major_formatter(ticker.FixedFormatter(ticklab_str))
        ax7.annotate('size 1"',(1.,1.02),xycoords='axes fraction',ha='center')
        
        # draw the 8th axis (angle 1kpc)
        ax8 = fig.add_axes([x0[0]+7*xw,y0[0],xs,ys])
        setup(ax8)
        ax8.set_ylim(ax1.get_ylim())
        ax8.set_yscale(ax1.get_yscale())
        ticklab1 = np.array([1.,0.5,0.3,0.2,0.15,0.12])
        ticklab1 = ticklab1[ticklab1>self.angle_z(self.zext)]
        tickloc1 = [self.zAngle1(a) for a in ticklab1]
        ticklab2 = ticklab1[::-1]
        tickloc2 = [self.zAngle2(a) for a in ticklab2]
        ticklab = np.concatenate((ticklab1,ticklab2))
        tickloc = np.concatenate((tickloc1,tickloc2))
        ticklab_str = ["%.2f"%t for t in ticklab]
        ax8.yaxis.set_major_locator(ticker.FixedLocator(tickloc))
        ax8.yaxis.set_major_formatter(ticker.FixedFormatter(ticklab_str))
        ax8.annotate('angle 1kpc',(1.,1.02),xycoords='axes fraction',ha='center')
        
        # draw the 9th axis (z)
        ax9 = fig.add_axes([x0[0]+8*xw,y0[0],xs,ys])
        setup(ax9)
        ax9.set_ylim(ax1.get_ylim())
        ax9.set_yscale(ax1.get_yscale())
        ax9.yaxis.set_major_locator(ticker.FixedLocator(ticklocz))
        ax9.yaxis.set_major_formatter(ticker.FixedFormatter(ticklabz_str))
        ax9.annotate('z',(1.,1.02),xycoords='axes fraction',ha='center')	


        if 'z' in keywords.keys():
            zz = keywords['z']
            print("z: %f"%zz)
            print("H(z) = %f km/s/Mpc"%self.Hz(zz))
            print("comoving radius r(z) = %f Mpc"%self.rz(zz))
            print("1 arcsec size(z) = %f kpc"%self.size_z(zz))
            print("angle 1kpc(z) = %f arcsec"%self.angle_z(zz))
            print("age(z) = %f Gyr"%self.age_z(zz))
            print("lookback time(z) = %f Gyr"%(self.age0 - self.age_z(zz)))

            ax9.annotate("", xy=(1,zz),xytext=(-79.,zz),xycoords='data',textcoords='data' \
                ,arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=0.",color='magenta'))

        plt.show()
        # save the image
        if 'out' in keywords.keys():
            print('saved to %s'%keywords['out'])
            fig.savefig(keywords['out'])
