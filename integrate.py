#! /usr/bin/python

import pyximport; pyximport.install()
import vegas
import integrand as INT
import numpy as np
#import time

import math
import scipy
from scipy import special

pi=np.pi

"""
Defining prefactors
"""

def multigamma(x, N):
	return scipy.exp(special.multigammaln(x,N))

def Wn_calc(N):
	return (pi**(3.0*(N-1.0)/2.0))/multigamma((N-1.0)/2.0,3)

def Pn_calc(N,gamma):
	ret = (1/((2*pi)**(2.0*N+3))) * ((5.0**(5.0/2.0))*(3**(3.0*N/2.0+3.0))/(2*(1.0-gamma**2)**(1.0/2.0)))
	return ret

def An_calc(N):
	ret = (2*pi**(N/2.0))/special.gamma(N/2.0)
	return ret

def Sn_calc(N,gamma,nu,s1,s2):
	ret = ((s2**3)/(s1**3))*(nu**(N-1.0)) * ((nu/gamma)**((3.0/2.0)*(N-1.0))) * scipy.exp(-1.0*(nu**2)/2.0)
	return ret

def integrate(N,nu,gamma,s1,s2,npoints):	

	"""
	Calculating prefactors
	"""
	An = An_calc(N)
	Pn = Pn_calc(N,gamma)


	Wn = Wn_calc(N) 
	Sn = Sn_calc(N,gamma,nu,s1,s2)


	domain=[[-pi/2,pi/2],[-pi/2,pi/2],[-pi/2,pi/2],[0,pi/2],[0,pi/2],[0,pi/2],[-pi/2,pi/2],[-pi/2,pi/2],[-pi/2,pi/2]]

	f=INT.f_cython(dim=9,N=N,nu=nu,gamma=gamma)
		
	integ = vegas.Integrator(domain, nhcube_batch=1000)
	
	integ(f, nitn=10, neval=npoints)
	vecresult = integ(f, nitn=10, neval=npoints)

	retlist = np.zeros((5,2))
	prefac = (Wn*8.0/6.0)*(2*pi**2*An*Pn*Sn)
	for i in range(0,vecresult.shape[0]):
		retlist[i,0] = prefac*vecresult[i].mean
		retlist[i,1] = prefac*vecresult[i].sdev

	return retlist
