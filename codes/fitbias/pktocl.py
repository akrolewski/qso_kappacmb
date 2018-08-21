import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as interpolate
from scipy.optimize import minimize
import matplotlib.pyplot as plt
plt.switch_backend('agg')

#local codes
from cmbcross_lensingZ import LensingZ
import cosmology
from pk_pregen import PkZCLEFT


class setup():
    dpath = '../../data/'
    dgenpath = '../../data/data-generated/'
    cosmo = cosmology.Cosmology(M=0.295, pfile=dpath+'pklin_ds14.dat')

    dndz = np.loadtxt(dgenpath + 'cls_sampled/dndz_qso_pdb_lf.txt').T
    idndz = interpolate(dndz[0], dndz[1])

    
stp = setup()


def chi2(p, lptz, lenz, cldata, error=None, lmin=10, lmax=1000):
    """   
    The chi^2 given the parameters, p.
    The behavior is dictated by the length of p.
    """
    # Evaluate the theory.
    ell = cldata[0]
    pk = lptz(p, auto=True)
    clqq = lenz.clzeff(pk, auto=True)    
    clqq = interpolate(lenz.l, clqq)(ell)
    
    # Set the weights
    wt  = np.exp(-(ell/lmax)**2)
    wt  = 1.0-ell/lmax
    wt[ell<lmin]=0
    wt[ell>lmax]=0

    # and compute chi^2.
    if error is None: error = np.ones_like(ell)
    c2 = np.sum((clqq-cldata[1])**2 *wt / error)
    
    # and add a prior on sn1 and sn2 to drive them to zero.
    # ADD prior on bs2 as well
    #p2 = 1.*p[4]**2 + 0.5*(p[6]/G.snerr)**2 #+ 100.*p[2]**2.
    p2 = 0 
    return(c2+p2)
    #



def fit_example(zmin, zmax, sigma=0.5, verbose=False):
    """    
    An example fitter.
    """
    # Clear out the past fits.
    pname = ["b1","b2","bs","bn","alpha","sn"]

    z0 = (zmin + zmax)/2.

    #data
    clqqfile = stp.dgenpath + \
               'cls_sampled/cl_qq_sampled_z%.2f_%.2f/cl_qq_z%.2f_%.2f_darksky_sigma%.2f.txt'%(zmin, zmax, zmin, zmax, sigma)
    cldata = np.loadtxt(clqqfile).T
    error = cldata[2]
    lmin, lmax = 10, 2000
    
    #theory
    lenz = LensingZ(z0, dndz=stp.idndz, cosmo=stp.cosmo)
    lptz = PkZCLEFT(zmin, zmax, z0=z0)

    #init
    pars = lptz.pars(b1=2.0, b2=0, alpha=0, sn=0)
    print(pars)
    print(chi2(pars, lptz, lenz, cldata, cldata[2]))

    #minimize
    solvopt = {'disp':True,'maxiter':5000,'maxfev':15000,'xtol':1e-4}
    tomin = lambda p : chi2(p, lptz, lenz, cldata, error=cldata[2], lmin=lmin, lmax=lmax)
    res = minimize(tomin,pars,method='nelder-mead',options=solvopt)
    #res = minimize(chi2,pars, lptz, lenz, cldata, cldata[2], method='nelder-mead',options=solvopt)
    print(res.x)
    print(chi2(res.x, lptz, lenz, cldata, cldata[2]))
    
    #plot cl
    clgginit = lenz.clzeff(lptz(pars, auto=True), auto=True)
    clgg = lenz.clzeff(lptz(res.x, auto=True), auto=True)
    plt.figure()
    plt.plot(lenz.l, clgginit)
    plt.plot(lenz.l, clgg)
    plt.plot(cldata[0], cldata[1], '--' , marker='.')
    plt.axvline(lmin, lw=0.5, color='gray')
    plt.axvline(lmax, lw=0.5, color='gray')
    plt.loglog()
    plt.savefig('tmp.png')

    #plot ratio
    plt.figure()
    label = ''
    for i in range(len(res.x)): label = label+pname[i] + ' = %0.2f\n'%res.x[i]
    print(label)
    plt.plot(cldata[0], cldata[1]/interpolate(lenz.l, clgg)(cldata[0]), label=label, marker='.')
    plt.axhline(1, lw=0.5, color='k')
    plt.axvline(lmin, lw=0.5, color='gray')
    plt.axvline(lmax, lw=0.5, color='gray')
    plt.xscale('log')
    plt.ylim(0.75, 1.25)
    plt.legend()
    plt.savefig('tmp2.png')
    

    
if __name__=="__main__":
    fit_example(zmin=1.5, zmax=2.0, verbose=False)
    #






##
##def fit_example(verbose=False):
##    """    
##    An example fitter.
##    """
##    # Clear out the past fits.
##    pname = ["b1","b2","bs","alpha","sn"]
##    fname = "./fit_results.log"
##    ff= open(fname,"w")
##    ff.write("# Results of fit to Clqq.\n")
##    ff.write("# "+str(datetime.datetime.now()).split('.')[0]+'\n')
##    st= "# iz"
##    for p in pname:
##        st += " %12s"%p
##    ff.write(st+" %12s\n"%"chi2")
##    ff.close()
##    # Some starting values for the fitter
##    b1 = {100: 0.42, 150: 0.85, 200: 1.60, 250:2.20, 300: 3.10}
##    kmax_x = [0.30,0.30,0.40,0.40,0.40]
##    kmax_a = [0.35,0.25,0.30,0.25,0.25]
##    if lo:
##        b1 = {100: 0.3, 150: 0.5, 200: 1.2, 250:1.80, 300: 2.50}
##        kmax_x = [0.45,0.30,0.55,0.40,0.40]
##        kmax_a = [0.40,0.25,0.50,0.25,0.40]
##    for ii,iz in enumerate([100,150,200,250,300]):
##        G.lpt    = PK.PkCLEFT(db+"PTdata/ps00_hh_RunPB_46_z%03d.dat"%iz)
##        G.kfit   = [1.0,1.0]
##        G.chi2wt = [1.0,1.0]
##        G.snerr  = 0.02*sn[ii]
##        print(ii,iz,G.snerr)
##        #
##        nbody = np.loadtxt(db+"hm_z%03d.pkr"%iz)
##        if lo:
##            nbody = np.loadtxt(db+"hm_z%03d_lo.pkr"%iz)
##        G.dat_kk = nbody[5:,0]
##        G.dat_dx = nbody[5:,1]
##        G.kfit[0]= kmax_x[ii]
##        nbody = np.loadtxt(db+"hh_z%03d.pkr"%iz)
##        if lo:
##            nbody = np.loadtxt(db+"hh_z%03d_lo.pkr"%iz)
##        G.dat_da = nbody[5:,1] - sn[ii]*nbody[5:,0]**3/2/np.pi**2
##        G.kfit[1]= kmax_a[ii]
##        # Use the Nelder-Mead method, because the other methods start
##        # using really odd values of the parameters.
##        pars    = [b1[iz],0.0,0.0,0.0,1.0,0.0,1.0]
##        solvopt = {'disp':True,'maxiter':5000,'maxfev':15000,'xtol':1e-4}
##        res = opt.minimize(chi2,pars,method='nelder-mead',options=solvopt)
##        # Optionally print the best fit and fit information.
##        #print(res)
##        p = res.x
##        p[4],p[6]=0.0,0.0	# Zero the shot noise terms.
##        ff= open(fname,"a")
##        st= " %03d"%iz
##        for x in p:
##            st += " %12.4e"%x
##        ff.write(st+" %12.4e\n"%res.fun)
##        ff.close()
##        if verbose:
##            print("\nMatch @ iz=",iz,"\n")
##            kk,px=G.lpt(p[0],p[1],p[2],0,p[3],p[4],auto=False)
##            kk,pa=G.lpt(p[0],p[1],p[2],0,p[5],p[6],auto=True)
##            delx =kk**3*px/(2*np.pi**2)
##            dela =kk**3*pa/(2*np.pi**2)
##            # Now interpolate the theory onto the N-body samples.
##            thy_x = np.interp(G.dat_kk,kk,delx)
##            thy_a = np.interp(G.dat_kk,kk,dela)
##            for i in range(G.dat_kk.shape[0]):
##                print("%3d %8.5f %12.4e %12.4e %8.3f %12.4e %12.4e %8.3f"%\
##                     (i,G.dat_kk[i],G.dat_dx[i],thy_x[i],\
##                     100*(G.dat_dx[i]-thy_x[i])/G.dat_dx[i],\
##                     G.dat_da[i],thy_a[i],\
##                     100*(G.dat_da[i]-thy_a[i])/G.dat_da[i]))
##    #
##
##
