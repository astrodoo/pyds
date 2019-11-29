"""
filename: 
     bhscale.py
 
PURPOSE:
     compute the distance scale given the black hole properties

Written by:
     Doosoo Yoon
     University of Amsterdam
   
History:
     Written, 25 January 2019
"""
from __future__ import print_function
import numpy as np
from . import astrounit as unit

class bhscale:
    """
    compute the distance scale given the black hole properties

    keywords -- mbh: black hole mass in solar mass unit (default= 4.1e6)
	        distance: distance to the target in kpc unit (default= 8.1)
                spin: spin of black hole (+:prograde, -:retrograde) (default= 0.9375)

    outputs   -- self.rg: the gravitational radius
                 self.rh: the radius of event horizon in the unit of rg
                 self.risco: the radius of inner stable circular orbit in the unit of rg
    """

    def __init__(self, mbh=4.1e6, distance=8.1, spin=0.9375, rad_eff=0.1, cgs=False, silent=True):
        self.mbh = mbh
        self.dist_kpc = distance
        self.spin = spin

        self.rg = self.calc_rg()
        self.rh = self.calc_rh()

        dist_cgs = self.dist_kpc*1e3*unit.pc
        unit_rad2mias = 180./np.pi*60.*60.*1e6    # rad to micro-arcsec
        self.rh_mias = self.rg*self.rh / dist_cgs * unit_rad2mias
        self.rg_mias = self.rg / dist_cgs * unit_rad2mias
 
        self.risco = self.calc_risco()
        self.tg = self.calc_tg()

        self.Ledd = self.calc_Ledd(cgs=cgs)
        self.Mdotedd = self.calc_Mdotedd(rad_eff=rad_eff,cgs=cgs)

        self.cgs = cgs

        if (not silent):
            print('BH Params: Mbh=%e M_sun, spin=%f, Distance=%e kpc'%(self.mbh,self.spin,self.dist_kpc) )
            print('rg: %e cm'%self.rg)
            print('rg_mias: %f micro-arcsec'%self.rg_mias)
            print('rh: %f rg'%self.rh)
            print('rh_mias: %f micro-arcsec'%self.rh_mias)
            print('risco: %f rg'%self.risco)

            if (cgs):
                print('L_edd = %e erg s^-1'%self.Ledd)
                print('Mdot_edd = %e g s^-1'%self.Mdotedd)
            else:
                print('L_edd = %e L_sun'%self.Ledd)
                print('Mdot_edd = %e M_sun/year'%self.Mdotedd)


    def calc_rg(self):
        """ 
        calculate of gravitational radius in cgs
        """
        Mbh_cgs = self.mbh*unit.msun

        return unit.g*Mbh_cgs/unit.c/unit.c


    def calc_tg(self):
        """
        calculate of gravitational time unit in second 
        """
        Mbh_cgs = self.mbh*unit.msun
        return unit.g*Mbh_cgs / unit.c/unit.c/unit.c
    

    def calc_rh(self):
        """
        calculate the event horizon radius in rg unit
        """
        return 1.+np.sqrt(1.-self.spin*self.spin)
    

    def calc_risco(self):
        """
        calculate the radius of inner most stable circular orbit in rg unit
        """
        z1 = 1. + np.power(1.-self.spin*self.spin,1./3.) \
            * ( np.power(1.+self.spin,1./3.) + np.power(1.-self.spin,1./3.) ) 
        z2 = np.sqrt(3.*self.spin*self.spin + z1*z1)
        
        if self.spin < 0.:
            risco = 3.+z2 + np.sqrt( (3.-z1) * (3. + z1 + 2.*z2) )
        else:
            risco = 3.+z2 - np.sqrt( (3.-z1) * (3. + z1 + 2.*z2) )
        
        return risco


    def calc_torb(self, rr): 
        """
        calculate the orbital time in seconds with given radius in cm
        """

        vkep = np.sqrt(unit.g*self.mbh*unit.msun/rr)

        return 2.*np.pi*rr/vkep

    def calc_torb2r(self, tt): 
        """
        calculate the radius in cm with given orbital time in second
        """

        return np.power(tt*0.5/np.pi, 2./3.) * np.power(unit.g*self.mbh*unit.msun,1./3.)
 

    def calc_torbg(self, rr): 
        """
        calculate the orbital time in seconds with considering the BH spin with given radius in cm
        Formular is given in Hamaus+09
        """

        return 2.*np.pi*self.rg/unit.c * ( np.power(rr/self.rg, 1.5) + self.spin )


    def calc_torbg2r(self, tt): 
        """
        calculate the radius in cm with given orbital time in second with the consideration of BH spin
        Formular is given in Hamaus+09
        """

        return self.rg * np.power( unit.c*tt / (2.*np.pi*self.rg) - self.spin, 2./3.)
 

    def calc_Ledd(self,cgs=False): 
        """
        calculate the eddington luminosity
        
        keywords:
            cgs - if given, output is in cgs unit, otherwise the unit is L_sun (by default)
        """

        Ledd = 4.*np.pi*unit.g*self.mbh*unit.msun*unit.mp*unit.c/unit.sigmaT

        if (cgs):
            return Ledd
        else:
            return Ledd/unit.lsun


    def calc_Mdotedd(self,rad_eff=0.1,cgs=False): 
        """
        calculate the eddington mass accretion rate 

        keywords:
            rad_eff - radiative efficiency. default is 0.1
            cgs - if given, output is in cgs unit, otherwise the unit is M_sun/year (by default)
        """

        Mdotedd =  4.*np.pi*unit.g*self.mbh*unit.msun*unit.mp/rad_eff/unit.c/unit.sigmaT

        if (cgs):
            return Mdotedd
        else:
            return Mdotedd/unit.msun*unit.year

        
    def draw(self,rgmin=1., rgmax=1e3, log=True, hline=None, **keywords):
        """
        Draw the scaled axes with given black hole parameters

        keywords -- rgmin: minimum of gravitational radius (default= 1.)
                    rgmax: maximum of gravitational radius (default= 1e3)
                    log: logarithmic scale of axes (default= True)
                    hline: draw horizontal line in rg unit (default= None)

        **keywords -- out: if given, the plot will be saved to the designated output file (type = string) 
        """
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker

        dist_cgs = self.dist_kpc*1e3*unit.pc
        unit_rad2mias = 180./np.pi*60.*60.*1e6    # rad to micro-arcsec

        # Setup a plot such that only the bottom spine is shown
        def setup(ax):
            ax.spines['top'].set_color('none')
            ax.spines['bottom'].set_color('none')
            ax.xaxis.set_major_locator(ticker.NullLocator())
            ax.spines['left'].set_color('none')
            ax.yaxis.set_ticks_position('right')
            ax.tick_params(axis='y',which='major',length=10,direction='inout')
            ax.tick_params(axis='y',which='minor',length=5, direction='inout')
        
        #figure setting
        x0 = [30.,50.]; y0 = [20.,70.]
        ys = 300.
        xw = 50.; xs = xw/10.
        naxes = 6
        winxs = np.sum(x0)+(naxes-1)*xw
        winys = np.sum(y0) + ys
        
        x0 = x0/winxs; y0 = y0/winys
        xs = xs/winxs; ys = ys/winys
        xw = xw/winxs
        
        fig = plt.figure(figsize=(15.*winxs/winys,15))
        
        # draw the 1st axis (rg)
        ax1 = fig.add_axes([x0[0],y0[0],xs,ys])
        setup(ax1)
        ax1.set_ylim(rgmin,rgmax)
        if (log):
            ax1.set_yscale('log')
        else:
            ax1.set_yscale('linear')
        ax1.annotate(r'$r_{g}$',(1.,1.02),xycoords='axes fraction',ha='center')
        
        # draw the 2nd axis (pc)
        pcmin = self.rg*rgmin / unit.pc
        pcmax = self.rg*rgmax / unit.pc
        ax2 = fig.add_axes([x0[0]+xw,y0[0],xs,ys])
        setup(ax2)
        ax2.set_ylim(pcmin,pcmax)
        ax2.set_yscale(ax1.get_yscale())
        ax2.annotate('pc',(1.,1.02),xycoords='axes fraction',ha='center')
        
        # draw the 3rd axis (au)
        aumin = self.rg*rgmin / unit.au
        aumax = self.rg*rgmax / unit.au
        ax2 = fig.add_axes([x0[0]+2*xw,y0[0],xs,ys])
        setup(ax2)
        ax2.set_ylim(aumin,aumax)
        ax2.set_yscale(ax1.get_yscale())
        ax2.annotate('au',(1.,1.02),xycoords='axes fraction',ha='center')
        
        # draw the 4th axis (micro-arcsecond)
        dist_cgs = self.dist_kpc*1e3*unit.pc
        unit_rad2mias = 180./np.pi*60.*60.*1e6    # rad to micro-arcsec
        angmin = self.rg*rgmin / dist_cgs * unit_rad2mias
        angmax = self.rg*rgmax / dist_cgs * unit_rad2mias
        ax3 = fig.add_axes([x0[0]+3*xw,y0[0],xs,ys])
        setup(ax3)
        ax3.set_ylim(angmin,angmax)
        ax3.set_yscale(ax1.get_yscale())
        ax3.annotate(r'$\mu as$',(1.,1.02),xycoords='axes fraction',ha='center')
        
        # draw the 5th axis (v_esc escape velocity)
        v_esc_lmin = np.sqrt(2.*unit.g*self.mbh*unit.msun / (self.rg*rgmin)) / 1e5   # km/s
        v_esc_lmax = np.sqrt(2.*unit.g*self.mbh*unit.msun / (self.rg*rgmax)) / 1e5   # km/s
        ax5 = fig.add_axes([x0[0]+4*xw,y0[0],xs,ys])
        setup(ax5)
        ax5.set_ylim(v_esc_lmin,v_esc_lmax)
        ax5.set_yscale(ax1.get_yscale())
        ax5.annotate(r'$v_{esc}$ (km/s)',(1.,1.02),xycoords='axes fraction',ha='center')

        # draw the 6th axis (rg)
        ax6 = fig.add_axes([x0[0]+5*xw,y0[0],xs,ys])
        setup(ax6)
        ax6.set_ylim(ax1.get_ylim())
        ax6.set_yscale(ax1.get_yscale())
        ax6.annotate(r'$r_{g}$',(1.,1.02),xycoords='axes fraction',ha='center')
        
        # mark the line for rh and risco
        ax6.annotate("", xy=(1,self.rh),xytext=(-49.,self.rh),xycoords='data',textcoords='data' \
                    ,arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=0.",color='magenta'))
        ax1.annotate(r"r$_h$",(0.,self.rh*1.05),xycoords='data',color='magenta')
        
        ax6.annotate("", xy=(1,self.risco),xytext=(-49.,self.risco),xycoords='data',textcoords='data' \
                    ,arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=0.",color='purple'))
        ax1.annotate(r"r$_{\rm isco}$",(0.,self.risco*1.05),xycoords='data',color='purple')

        # mark the request line and print the informations
        print('horizontal line: %f (rg) = %f (mu-as) = %f (pc) = %f (au)'% \
                (hline, hline*self.rg/dist_cgs*unit_rad2mias, \
                 hline*self.rg/unit.pc, hline*self.rg/unit.au))
        if not(hline is None):
            ax6.annotate("", xy=(1,hline),xytext=(-49.,hline),xycoords='data',textcoords='data' \
                        ,arrowprops=dict(arrowstyle="-",connectionstyle="arc3,rad=0.",color='blue'))
 

        
        # make an annotation for the black hole parameters
        ParamStr = r'M$_{\rm BH}$=%5.2e M$_{\odot}$, Distance=%5.2e kpc, spin=%4.2f'%(self.mbh, self.dist_kpc,self.spin)
        ax1.annotate(ParamStr,(1.,1.07),xycoords='axes fraction',ha='left',fontsize='large')
        
        plt.show()
	# save the image
        if 'out' in keywords.keys():
            print('saved to %s'%keywords['out'])
            fig.savefig(keywords['out'])
