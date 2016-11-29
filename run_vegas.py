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
nfields=sys.argv[2]

print "NPTS: ",npoints_arg
print "NFIELDS: ",npoints_arg

"""
Define the true expectation value of lam1*lam2*lam3 without the heaviside
function in lam3, but with ordering in lambdas. This is to be used as
a verification that the vegas algorithm is performing well. Beacuse the integrand
is identical to the real case, but only the limits of integration change,
I expect good performance of vegas relative to this analytic form to be
a strong predictor of accuracy in the true estimates.
"""

def truint(myN,mynu,mys0,mys1):
	ret = ((mys1**3)/(mys0**3))*(mynu**(myN-4.0))*np.exp((-1.0/2.0)*mynu**2) * ((-6.0 +11.0*myN -6.0*myN**2 + myN**3 - 3.0*((myN-1.0)**2)*(mynu**2) + 3.0*myN*(mynu**4) - (mynu**6))/((2**((myN+1.0)/2.0))*3.0*(3.0**(1.0/2.0))*(np.pi**(3.0/2.0)) * special.gamma(myN/2.0)))
	return ret


"""
Define parameters, as well as the set of values nu should loop over,
and the resolution.
"""
s0=2.0
s1=2.0
gamma=0.5
s2 = (s1**2)/(gamma*s0)
N=float(nfields)
print "NFIELDS: ",N


nsamples=10
nu=np.linspace(0.05,7.0,nsamples)
#nu=np.linspace(1,5.0,nsamples)

npoints=int(float(npoints_arg))
print "NPOINTS: ",npoints


"""
Initialize output vectors.
"""

means=np.zeros(nsamples)
sdevs=np.zeros(nsamples)

true_vals=np.zeros(nsamples)

pct_sdevs=np.zeros(nsamples)
errs=np.zeros(nsamples)
abs_errs=np.zeros(nsamples)
pct_errs=np.zeros(nsamples)
ratios=np.zeros(nsamples)

"""
Run the vegas algorithm over varying nu.
"""

for samp in tqdm(range(nsamples)):
	ret = integrate.integrate(N,nu[samp],gamma,s1,s2,npoints)
	
	signed = 0.0
	for i in range(1,ret.shape[0]):
		signed += (-1.0)**(i-1)*ret[i,0]
		
	means[samp]=ret[0,0]
	#means[samp]=signed
	sdevs[samp]=ret[0,1]


	true_vals[samp]=truint(N,nu[samp],s0,s1)
	pct_sdevs[samp] = abs(100.0*(sdevs[samp])/means[samp])
	errs[samp] = (means[samp]-true_vals[samp])
	abs_errs[samp] = abs(means[samp]-true_vals[samp])
	pct_errs[samp] = abs(100.0*(means[samp]-true_vals[samp])/true_vals[samp])
	ratios[samp] = means[samp]/true_vals[samp]


"""
Format an output table for spot-checking!
"""

t = PrettyTable(['Nu','Means','True','Sdev','PctSdev','Err','AbsErr','PctErr','Ratio'])

for i in range(len(nu)):
	rowstrs=[]
	rowstrs.append('{:06.4f}'.format(nu[i]))
	rowstrs.append('{:06.4f}'.format(means[i]))
	rowstrs.append('{:06.4f}'.format(true_vals[i]))
	rowstrs.append('{:06.4f}'.format(sdevs[i]))
	rowstrs.append('{:06.4f}'.format(pct_sdevs[i]))
	rowstrs.append('{:06.4f}'.format(errs[i]))
	rowstrs.append('{:06.4f}'.format(abs_errs[i]))
	rowstrs.append('{:06.4f}'.format(pct_errs[i]))
	rowstrs.append('{:06.4f}'.format(ratios[i]))
	t.add_row(rowstrs)

print t



"""
Write these arrays out to an npz file for plotting.
"""

np.savez("signed_"+str(nfields)+"field_"+str(npoints_arg)+".npz",npoints=npoints,N=N,s0=s0,s1=s1,gamma=gamma,nu=nu,means=means,sdevs=sdevs,true_vals=true_vals,pct_sdevs=pct_sdevs,errs=errs,abs_errs=abs_errs,pct_errs=pct_errs)


