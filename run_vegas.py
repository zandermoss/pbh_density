#! /usr/bin/python

import sys
import numpy as np
from scipy import special
import integrate
from prettytable import PrettyTable
from tqdm import tqdm
import math
import scipy
from scipy import special
pi=math.pi



npoints_arg=sys.argv[1]
nfields_arg=sys.argv[2]

npoints = float(npoints_arg)
nfields = int(nfields_arg)

print "NPOINTS: ",npoints
print "NFIELDS: ",nfields

"""
An implementation of the analytic signed number density calculated in the paper. The integrand is almost identical to those of each of the stationary point integrands (minimum, saddles, maximum), minus an absolute value and heaviside functions. The numerical estimate of the signed number density is used to refine the VEGAS grid, and comparisons between this analytic result and the numerical estimate provide sanity checks.
"""

def truint(myN,mynu,mys0,mys1):
	ret = ((mys1**3)/(mys0**3))*(mynu**(myN-4.0))*np.exp((-1.0/2.0)*mynu**2) * ((-6.0 +11.0*myN -6.0*myN**2 + myN**3 - 3.0*((myN-1.0)**2)*(mynu**2) + 3.0*myN*(mynu**4) - (mynu**6))/((2**((myN+1.0)/2.0))*3.0*(3.0**(1.0/2.0))*(np.pi**(3.0/2.0)) * special.gamma(myN/2.0)))
	return ret

"""
Define parameters s0 and s1 (from field theory power spectrum).
"""
gamma=0.01
s0=1.0
s1=1.0


"""
Define the grid on which to sample nu-gamma parameter space.
"""
nsamples_nu = 20
#nsamples_gamma = 1

nu=np.linspace(0.05,6.0, nsamples_nu)
#gamma=np.linspace(0.00,1.0, nsamples_gamma)


"""
Initialize output vectors.
"""

#The signed number density is gamma-independent
signed_analytic = np.zeros(nsamples_nu) 
signed_mean = np.zeros(nsamples_nu)
signed_mean = np.zeros(nsamples_nu)
signed_sdev = np.zeros(nsamples_nu)
#Four entries for these densities: one for each type of stationary point
means=np.zeros((nsamples_nu,4))
sdevs=np.zeros((nsamples_nu,4))

"""
Run the vegas algorithm over the specified grid in parameter space.
"""

for nu_samp in tqdm(range(nsamples_nu)):
	signed_analytic[nu_samp] =truint(nfields,nu[nu_samp],s0,s1)

	#Calculating sigma_2 for convenience in integration script.
	s2 = (s1**2)/(gamma*s0)

	#Calling the integration script. This script is the interface to VEGAS
	ret = integrate.integrate(nfields,nu[nu_samp],gamma,s1,s2,npoints)

	#The first entry in the returned vector is the signed number density mean estimate.
	#The VEGAS grid is refined using the signed number density.	
	signed_mean[nu_samp] = ret[0,0]		
	signed_sdev[nu_samp] = ret[0,1]		

	#The next four values are the mean estimates of the four types of stationary points.
	#They are ranked in order of the number of negative Hessian eigenvalues.
	#The first entry has zero negative eigenvalues: it is a density of minima.
	#The second and third entries are densities of saddle points with one and two
	#negative hessian eigenvalues, respectively.
	#The fourth entry has three negative eigenvalues: it is a density of maxima. 
	means[nu_samp,:]=ret[1:,0]
	sdevs[nu_samp,:]=ret[1:,1]


"""
Write these arrays out to an npz file for plotting.
"""

np.savez("densities_"+str(nfields)+"fields_"+str(npoints_arg)+".npz",npoints=npoints,nfields=nfields,s0=s0,s1=s1,gamma=gamma,nu=nu,signed_analytic=signed_analytic,signed_mean=signed_mean,signed_sdev=signed_sdev,means=means,sdevs=sdevs)


