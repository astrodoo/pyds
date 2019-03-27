import numpy as np

def func(pr,ppre,pl,vr,vl,csr,csl,gamma):

    tmp1=np.sqrt((gamma+1.)/(2.*gamma)*(ppre/pr-1.)+1.)
    tmp2=1.+(gamma-1.)/(2.*csl)*(vl-vr-csr/gamma*(ppre/pr-1.)/tmp1)
    re_func=pl*np.power(tmp2,2.*gamma/(gamma-1.))-ppre

    return re_func


def sec(ppre0,p,v,cs,gamma):
    MAXIT = 30
    eps = 1e-8

    f0 = func(p[1],ppre0[0],p[0],v[1],v[0],cs[1],cs[0],gamma)
    f1 = func(p[1],ppre0[1],p[0],v[1],v[0],cs[1],cs[0],gamma)

    if (np.abs(f0) < np.abs(f1)):
        re_sec = ppre0[0] 
        ppl    = ppre0[1]
        swap   = f0
        f0     = f1
        f1     = swap
    else:
        ppl    = ppre0[0]
        re_sec = ppre0[1]

    for j in range(MAXIT):
        dp = (ppl - re_sec)*f1 / (f1-f0)
        ppl = re_sec
        f0 = f1
        re_sec = re_sec + dp
        f1 = func(p[1],re_sec,p[0],v[1],v[0],cs[1],cs[0],gamma)

        if ((np.abs(dp) < eps) | (f1 == 0.)): return re_sec

    raise ValueError('ERROR: secant exceed maximum iterations')


def sodshockSolve(d=[1.,0.125], p=[1.,0.1], v=[0.,0.], gamma=1.66667, \
        time=0.245, resol=200, xscale=[0.,1.]):

    d = np.asarray(d); p = np.asarray(p); v = np.asarray(v); xscale = np.asarray(xscale)

    dd = np.zeros(resol); vv = np.zeros(resol)
    pp = np.zeros(resol); cs = np.zeros(resol); ee = np.zeros(resol)

    x = np.linspace(xscale[0],xscale[1],resol)

    xhf = (xscale[1]+xscale[0])/2.    # half distance
    
    cs0 = np.sqrt(gamma*p/d)
    
    #--------------------------------------------------------------------
    # pre-shock region
    #--------------------------------------------------------------------
    # initial guess
    ppre0 = [0.7,0.05]*p
    ppre  = sec(ppre0,p,v,cs0,gamma) 

    # compute pre-shock quantities
    cspre = cs0[1]*np.sqrt(ppre/p[1]*((gamma+1.)/(gamma-1.)+ppre/p[1]) \
            / (1.+(gamma+1.)/(gamma-1.) * ppre/p[1]))
    temp = np.sqrt((gamma+1.)/(2.*gamma)*(ppre/p[1]-1.)+1.)
    vpre = v[1]+cs0[1]/gamma*(ppre/p[1]-1.)/temp
    dpre = gamma * ppre / (cspre**2.)
    vsh  = v[1] + cs0[1]*temp   # shock speed

    #--------------------------------------------------------------------
    # contact discontinuity
    #--------------------------------------------------------------------
    pcd = ppre
    vcd = vpre
    dcd = d[0] * np.power(pcd/p[0] , 1./gamma)
    cscd = np.sqrt(gamma * pcd / dcd)

    #--------------------------------------------------------------------
    # position of waves
    #--------------------------------------------------------------------
    xsh   = xhf + vsh*time
    xcd   = xhf + vcd*time
    xtail = xhf + (vcd-cscd)*time
    xhead = xhf + (v[0]-cs0[0])*time

    #--------------------------------------------------------------------
    # computing solution as the position separately
    #--------------------------------------------------------------------
    
    reg1 = (x < xhead)
    reg2 = (x >= xhead) & (x < xtail)
    reg3 = (x >= xtail) & (x < xcd)
    reg4 = (x >= xcd)   & (x < xsh)
    reg5 = (x >= xsh)

    pp[reg1] = p[0]
    dd[reg1] = d[0]
    vv[reg1] = v[0]
    cs[reg1] = cs0[0]
    ee[reg1] = 1./(gamma-1.)*p[0]/d[0]

    if (time != 0.):
        vv[reg2] = 2./(gamma+1.)*((x[reg2]-xhf)/time+(gamma-1.)/2.*v[0]+cs0[0])
        cs[reg2] = vv[reg2]-(x[reg2]-xhf)/time
        pp[reg2] = p[0]*np.power(cs[reg2]/cs[0], 2.*gamma/(gamma-1.))
        dd[reg2] = gamma*pp[reg2]/(cs[reg2]**2.)
        ee[reg2] = 1./(gamma-1.)*pp[reg2]/dd[reg2]

    pp[reg3] = pcd
    dd[reg3] = dcd
    vv[reg3] = vcd
    cs[reg3] = cscd
    ee[reg3] = 1./(gamma-1.)*pcd/dcd

    pp[reg4] = ppre
    dd[reg4] = dpre
    vv[reg4] = vpre
    cs[reg4] = cspre
    ee[reg4] = 1./(gamma-1.)*ppre/dpre

    pp[reg5] = p[1]
    dd[reg5] = d[1]
    vv[reg5] = v[1]
    cs[reg5] = cs0[1]
    ee[reg5] = 1./(gamma-1.)*p[1]/d[1]

    results = {'x':x,'p':pp,'d':dd,'v':vv,'cs':cs,'e':ee \
            ,'xsh':xsh,'xcd':xcd,'xhead':xhead,'xtail':xtail,'vsh':vsh}

    return results


class sodshock:
    """
    solving sodshocktube problem in adiabatic equation of state
        
    keywords -- d: initial densities (float array [d0,d1]; default = [1.,0.125])
                p: initial pressures (float array [p0,p1]; default = [1.,0.1])
                v: initial velocities (float array [v0,v1]: default = [0.,0.])
                xscale: xrange (float array[x0,x1]: default = [0.,1.])
                gamma: adiabatic index (float; defaut = 1.66667)
                time: the desired time to calcuate the profile (float; defaut = 0.245)
                resol: the resolution of x-grid (int; default = 200)

    output -- self.d / self.p / self.e / self.v / self.cs / self.x / 
              self.xtail / self.xhead / self.xsh / self.xcd / self.vsh
    """
    def __init__(self,d=[1.,0.125], p=[1.,0.1], v=[0.,0.], gamma=1.66667, \
            time=0.245, resol=200, xscale=[0.,1.]):

        results = sodshockSolve(d=d,p=p,v=v,gamma=gamma,time=time,resol=resol,xscale=xscale)

        self.x     = results['x']
        self.d     = results['d']
        self.p     = results['p']
        self.v     = results['v']
        self.e     = results['e']
        self.cs    = results['cs']
        self.xsh   = results['xsh']
        self.xcd   = results['xcd']
        self.xhead = results['xhead']
        self.xtail = results['xtail']
        self.vsh   = results['vsh']
        self.time  = time 

