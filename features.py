#code to compute Fisher matrix and errors on feature parameters (and all other params)
#takes about a minute or two per Fisher matrix, so if you have many fruquencies you want to try, let it run on a cluster

#other models and parameters can be run also. Already existing parameters can be added to the list below
#They are divided into cosmo params and initial params. Initial params are those that are defined in Power_tilt.f90
#If you want to run a different mode, you will have to modify the pycamb file initialpower.py

import sys, platform, os
from matplotlib import pyplot as plt
import numpy as np
import math
import camb
from camb import model, initialpower
from scipy.misc import derivative as deriv


ellmin = 2
ellmax = 2500

#get params from CAMB:
pars = camb.CAMBparams()
pars.set_for_lmax(ellmax, lens_potential_accuracy=0);

#just to be sure:
pars.Reion.Reionization = True
pars.Reion.delta_redshift = 0.2
pars.AccurateReionization = True
pars.AccuratePolarization = True
pars.Reion.use_optical_depth = True
pars.Reion.fraction = -1
pars.Reion.redshift = 10


#parameters to vary: (you can simply add a parameter here, it should work)

cosmoparsLaTeX = ["H_0","\Omega_bh^2","\Omega_ch^2","\tau"]
cosmopars =  ['H0','ombh2','omch2','tau']
initialparsLaTeX = ['n_s','\log_{10} f','\delta_{n_s}','p_f','\delta_{\phi}','A_s']
initialpars = ['ns','log10_f','delta_ns','paxion_f','delta_Phi','As']

#observational paramaters for Planck:
fsky = 0.75
w0T = 30./np.sqrt(fsky);
w0P = 60. #sqrt(2.d0)*w0T
beam = 7.


print "Generating default power spectra (TT, TE and EE)"
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06)
pars.InitPower.set_params(ns=0.965, r=0, log10_f = -3, delta_ns = 0, paxion_f = 0.75, parameterization = 2 )
results = camb.get_results(pars)
powersdefault =results.get_cmb_power_spectra(pars)
Cls=powersdefault['total']

#label creator 
def createlatexlabel(x):
    si = "$\sigma_{"
    fi = "}$"
    finv1 = "$F^{-1/2}_{"
    fend1 = "}$"
    finv2 = "$\sqrt{F^{-1}_{"
    fend2 = "}}$"
    pltlabel = "Error on $"
    pltlabelend ="$"
    return si+x+fi, finv1+x+fend1, finv2+x+fend2, pltlabel+x+pltlabelend

#generating parameter reset tool. Make sure you add the parameter here as well. 
def resettemppars():
    tmp_initial_pars = {}
    tmp_initial_pars['ns'] = 0.965
    tmp_initial_pars['r'] = 0
    tmp_initial_pars['log10_f'] = -2
    tmp_initial_pars['delta_ns'] = 0.01 #this sets overall modulation of signal
    tmp_initial_pars['paxion_f'] = 0.75
    tmp_initial_pars['delta_Phi'] = math.pi
    tmp_initial_pars['As']  = 2.1e-09
   
    tmp_cosmo_pars = {}    
    tmp_cosmo_pars['H0'] = 67.5
    tmp_cosmo_pars['ombh2'] = 0.022
    tmp_cosmo_pars['omch2'] = 0.12
    tmp_cosmo_pars['mnu'] = 0.06
    tmp_cosmo_pars['omk'] = 0
    tmp_cosmo_pars['tau'] = 0.06
    
    return tmp_cosmo_pars, tmp_initial_pars


#wrap to obtain power spectra after changing initial conditions
def get_cmb_power_spectra_initial(x, dummy_tuple):
    #print dummy_tuple, type(dummy_tuple)
    if type(dummy_tuple)==tuple: tmp_initial_pars = dummy_tuple[0]
    else: tmp_initial_pars = dummy_tuple
    # Calculates the power spectrum , given some value x of some variable i, where i is in tmp_args
    tmp_initial_pars[tmp_initial_pars['which_param']] = x
    set_cosmology_tmp_initial_pars = {}
    for k in tmp_initial_pars:
        if k != 'which_param': set_cosmology_tmp_initial_pars[k] = tmp_initial_pars[k]
    #print "initial parameter  values: ",        set_cosmology_tmp_initial_pars
    pars.InitPower.set_params(**set_cosmology_tmp_initial_pars)
    #print "cosmological parameter values: ", tmp_cosmo_pars
    pars.set_cosmology(**tmp_cosmo_pars)
    results= camb.get_results(pars)
    powers = results.get_cmb_power_spectra(pars)

    return powers['total']

#wrap to obtain power spectra after changing cosmological parameters
def get_cmb_power_spectra_cosmo(x, dummy_tuple):
    #print dummy_tuple, type(dummy_tuple)
    if type(dummy_tuple)==tuple: tmp_cosmo_pars = dummy_tuple[0]
    else: tmp_cosmo_pars = dummy_tuple
    # Calculates the power spectrum , given some value x of some variable i, where i is in tmp_args
    tmp_cosmo_pars[tmp_cosmo_pars['which_param']] = x

    pars.InitPower.set_params(**tmp_initial_pars)
    #pars.InitPower.set_params(ns=0.965, r=0, log10_f = -2.5, delta_ns = 0.1, paxion_f = 0.75, delta_Phi = 0.1)
    #pars.set_cosmology(H0=H0)
    set_cosmology_tmp_cosmo_pars = {}
    for k in tmp_cosmo_pars:
        if k != 'which_param': set_cosmology_tmp_cosmo_pars[k] = tmp_cosmo_pars[k]
    #print "cosmological parameter  values: ",        set_cosmology_tmp_cosmo_pars
    pars.set_cosmology(**set_cosmology_tmp_cosmo_pars)
    results= camb.get_results(pars)
    powers = results.get_cmb_power_spectra(pars)
    #
    return powers['total']

#defining derivative function 
def simple_deriv(x0, delta, func=None, which_param='', arg_dict={}):
        arg_dict['which_param'] = which_param
        print "derivative w.r.t.", arg_dict['which_param']  
        ttt = (arg_dict,)    
        return deriv(func, x0, dx=delta, args=ttt)

#noise T and E in l*(l+1)N_l/2/pi/T_cmb^2

def CMBNoise(sigmab,w0,l):
    arcmin = 0.000290888
    CMBtemp = (2.7255E06)**2 #in muK^2 --> to make noise dimensionless
    fac = arcmin / 2.*np.sqrt(2.*np.log(2.))    
    return  l*(l+1)*w0**2*arcmin**2*np.exp(l**2*fac**2*sigmab**2)/CMBtemp/2./math.pi



#assume covariance to be diagnoal
def inversecovariance(l,Cls,sb,n0T,n0P):
    TTTT = (Cls[l,0]+CMBNoise(sb,n0T,l))**2
    TTTE = (Cls[l,0]+CMBNoise(sb,n0T,l))*Cls[l,2]
    TTEE = (Cls[l,2])**2
    TETE = (Cls[l,2]**2+(Cls[l,0]+CMBNoise(sb,n0T,l))*(Cls[l,0]+CMBNoise(sb,n0P,l)))/2.
    TEEE = (Cls[l,1]+CMBNoise(sb,n0P,l))*Cls[l,2] 
    EEEE = (Cls[l,1]+CMBNoise(sb,n0P,l))**2            
    CovMatrix =np.array([[TTTT,TTTE,TTEE],[TTTE,TETE,TEEE],[TTEE,TEEE,EEEE]])
    return  (2.*l+1)*np.linalg.inv(CovMatrix)/2.

#test
invmatrix = inversecovariance(10,Cls,beam,w0T,w0T)
print "inverse matrix for l = 10", invmatrix

#how many frequencies and over what range
Nfreq = 100
freq = np.linspace(-3.0,-1,Nfreq)

totparams = np.shape(cosmopars)[0] + np.shape(initialpars)[0]
print "Total # of paramaters", totparams, "of which cosmo", np.shape(cosmopars)[0], "and initial", np.shape(initialpars)[0] 

sigma1D = np.zeros((totparams,Nfreq))
sigmaM = np.zeros((totparams,Nfreq))

for n in range(0, Nfreq): 


    initialderivs = np.zeros((np.shape(initialpars)[0]),dtype=object)
    cosmoderivs = np.zeros((np.shape(cosmopars)[0]),dtype=object)

    #computing all derivatives
    #Note: default values parameters all have to be non-zero.. Will fix this. 
    for i in range(0,np.shape(cosmopars)[0]):
        tmp_cosmo_pars, tmp_initial_pars = resettemppars()
        #set Freq
        tmp_initial_pars[initialpars[1]] = freq[n]
        #Set param you want to vary
        tmp_cosmo_pars['which_param'] = cosmopars[i]
        par = tmp_cosmo_pars[cosmopars[i]]
        cosmoderivs[i] = simple_deriv( par, 1e-1*par, func=get_cmb_power_spectra_cosmo, which_param=cosmopars[i], arg_dict=tmp_cosmo_pars)    

    for i in range(0,np.shape(initialpars)[0]):
        tmp_cosmo_pars, tmp_initial_pars = resettemppars()
        #set Freq
        tmp_initial_pars[initialpars[1]] = freq[n]
        #Set param you want to vary
        tmp_initial_pars['which_param'] = initialpars[i]
        par = tmp_initial_pars[initialpars[i]]
        initialderivs[i] = simple_deriv( par, 1e-2*par, func=get_cmb_power_spectra_initial, which_param=initialpars[i], arg_dict=tmp_initial_pars)    

    #generate Fisher:
    MyFisher = np.zeros((totparams,totparams))

    for i in range(0,np.shape(cosmopars)[0]):
        for j in range(i,np.shape(cosmopars)[0]):
        
            MyFisher[i,j] = sum(cosmoderivs[i][L,0:3].dot(inversecovariance(L,Cls,beam,w0T,w0T).dot(cosmoderivs[j][L,0:3])) for L in range(ellmin,ellmax))
            MyFisher[j,i] = MyFisher[i,j]
        
    for i in range(np.shape(cosmopars)[0],totparams):
        for j in range(i,totparams):
        
            MyFisher[i,j] = sum(initialderivs[i-int(np.shape(cosmopars)[0])][L,0:3].dot(inversecovariance(L,Cls,beam,w0T,w0T).dot(initialderivs[j-int(np.shape(cosmopars)[0])][L,0:3])) for L in range(ellmin,ellmax))
            MyFisher[j,i] = MyFisher[i,j]
        
    for i in range(np.shape(cosmopars)[0],totparams):
        for j in range(0,np.shape(cosmopars)[0]):
        
            MyFisher[i,j] = sum(initialderivs[i-int(np.shape(cosmopars)[0])][L,0:3].dot(inversecovariance(L,Cls,beam,w0T,w0T).dot(cosmoderivs[j][L,0:3])) for L in range(ellmin,ellmax))
            MyFisher[j,i] = MyFisher[i,j]

    InverseFisher = np.linalg.inv(MyFisher)

    #compute errors and print them on screen:

    for i in range(0,totparams):
        sigmaM[i,n] = (InverseFisher[i,i])**(1./2.)
        sigma1D[i,n] = (MyFisher[i,i])**(-1./2.)
        if i < np.shape(cosmopars)[0]:
            print "Marginalized Error on", cosmopars[i],":", sigma1D[i,n]
            print "UnMarginalized Error on", cosmopars[i],":", sigmaM[i,n]
        else:
            print "Marginalized Error on", initialpars[i-int(np.shape(cosmopars)[0])],":", sigma1D[i,n]
            print "UnMarginalized Error on", initialpars[i-int(np.shape(cosmopars)[0])],":", sigmaM[i,n]


#plotting results
for i in range(0,totparams):
    if i < np.shape(cosmopars)[0]:
        yl,l1,l2, pltitle = createlatexlabel(cosmoparsLaTeX[i])
        storeas = cosmopars[i]+'_X_error.pdf'
    else:
        yl,l1,l2, pltitle = createlatexlabel(initialparsLaTeX[i-int(np.shape(cosmopars)[0])])
        storeas = initialpars[i-int(np.shape(cosmopars)[0])]+'_X_error.pdf' 
    plt.figure("fig1")
    plt.semilogy(freq,sigma1D[i], label=l1)
    plt.semilogy(freq,sigmaM[i], label=l2)
    plt.ylabel(yl,fontsize=20)
    plt.xlabel('$\log_{10} f$',fontsize=20)
    plt.legend(loc='best')
    plt.title(pltitle)
    
    plt.savefig(storeas)
    plt.clf()

