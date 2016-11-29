#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
from scipy import special



def truint(myN,mynu,mys0,mys1):
    ret = ((mys1**3)/(mys0**3))*(mynu**(myN-4.0))*np.exp((-1.0/2.0)*mynu**2) * ((-6.0 +11.0*myN -6.0*myN**2 + myN**3 - 3.0*((myN-1.0)**2)*(mynu**2) + 3.0*myN*(mynu**4) - (mynu**6))/((2**((myN+1.0)/2.0))*3.0*(3.0**(1.0/2.0))*(np.pi**(3.0/2.0)) * special.gamma(myN/2.0)))
    return ret



from pylab import rcParams
rcParams['figure.figsize'] = 10, 8
plt.rc('font', size=16)

#files = ["scalar_min_5field_1e6.npz","min_5field_1e6.npz","vectorized_min_0train_5field_1e6.npz"]
files = ["hess_min_5field_1e5.npz", "vectorized_realmin_mintrain_5field_1e5.npz", "vectorized_realmin_signtrain_5field_1e5.npz"] 

signed=False


f=[]
for filename in files:
	f.append(np.load(filename))

N=f[0]['N']
s0=f[0]['s0']
s1=f[0]['s1']
gamma=f[0]['gamma']

npoints=f[0]['npoints']

nu=f[0]['nu']
#pct_sdevs=f['pct_sdevs']

#xedge=[0.09,2.1]
xedge=[0.0,5.05]

def nucut(nu,xlo,xhi):
	return xlo<=nu and nu <= xhi

nucut_vec = np.vectorize(nucut)

cut_indices=np.where(nucut_vec(nu,xedge[0],xedge[1]))
print cut_indices

nu=f[0]['nu'][cut_indices]


sdevs=[]
means=[]
errs=[]
pct_errs=[]
abs_errs=[]

for el in f:
	sdevs.append(el['sdevs'][cut_indices])
	means.append(el['means'][cut_indices])
	errs.append(el['errs'][cut_indices])
	pct_errs.append(el['pct_errs'][cut_indices])
	abs_errs.append(el['abs_errs'][cut_indices])


"""
print pct_errs
print np.max(pct_errs)
"""

if signed:
	truex_nsamp=1000
	truex=np.linspace(xedge[0],xedge[1],1000)
	truef_vec=np.vectorize(truint)
	truey=truef_vec(N,truex,s0,s1)


axes=[]


axes.append(plt.subplot2grid((20,1),(6,0),rowspan=14,colspan=1))
axes[0].errorbar(nu,means[0],yerr=sdevs[0],marker="o",color='red',linestyle='None',lw=1.0,label="Scalar")
axes[0].errorbar(nu,means[1],yerr=sdevs[1],marker="o",color='green',linestyle='None',lw=1.0,label="Vectorized. Min Training")
axes[0].errorbar(nu,means[2],yerr=sdevs[2],marker="o",color='blue',linestyle='None',lw=1.0,label="Vectorized. Signed Training")
if signed:
	axes[0].plot(truex,truey,ls='-',lw=1.5,color='black',label="Analytic")
#axes[0].plot(xax_vector,pert_vec0,color="red",label="Standard Vacuum Solution",lw=2,ls='-')

#ax.set_yscale('log')
#ax.set_ylim([1e-1,1])
#ax.set_ylim([0,1.0])
axes[0].set_xlabel(r'$\nu$'+" (Field Standard Deviations)",fontsize=18)
axes[0].set_ylabel(r'$\langle \mathcal{N}_{min}(\nu)\rangle$',fontsize=18)
axes[0].set_xlim(xedge)
#axes[0].set_ylim([-10.0,3.0])
#axes[0].yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False))
#axes[0].yaxis.get_major_formatter().set_powerlimits((2, 3))

if sign
axes.append(plt.subplot2grid((20,1),(0,0),rowspan=5,colspan=10,sharex=axes[0]))

#axes[1].plot(xax_vector,np.absolute(np.subtract(pert_vec1,num_vec)),color="green",lw=2,ls='--')
axes[1].errorbar(nu,errs[0],yerr=sdevs[0],linestyle='None',marker="^",color="red",label="True Error")
axes[1].errorbar(nu,errs[1],yerr=sdevs[1],linestyle='None',marker="^",color="green")
axes[1].errorbar(nu,errs[2],yerr=sdevs[2],linestyle='None',marker="^",color="blue")
#axes[1].errorbar(nu,errs_uo,yerr=sdevs_uo,linestyle='None',marker="^",color="orange",label="True Error Un-Ordered")
axes[1].plot(truex,np.zeros(len(truex)),color="black")

#axes[1].set_title("Expectation Comparison: Analytic vs. Vegas")
axes[1].set_xlim(xedge)
axes[1].set_ylabel("Error")
#axes[1].set_ylim([0,100])
axes[1].yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False))
axes[1].yaxis.get_major_formatter().set_powerlimits((2, 3))
plt.setp(axes[1].get_xticklabels(), visible=False)






# Now add the legend with some customizations.
#legend = ax.legend(loc='center right', shadow=False)
#legend = ax.legend(bbox_to_anchor=(0.9, 0.6), bbox_transform=plt.gcf().transFigure, shadow=False)

h1, l1 = axes[0].get_legend_handles_labels()
h2, l2 = axes[1].get_legend_handles_labels()

#legend = axes[0].legend(h1+h2,l1+l2,bbox_to_anchor=(0.925, 0.4),shadow=False,fancybox=True)
legend = axes[0].legend(h1+h2,l1+l2,loc='lower right',shadow=False,fancybox=True)
#legend = axes[0].legend(h1+h2,l1+l2,bbox_to_anchor=(0.01, 0.7),bbox_transform=plt.gcf().transFigure, shadow=False,fancybox=True)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('1.0')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize(16)

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width

for ax in axes:
	ax.tick_params(axis='x', which='major', labelsize=16)
	ax.tick_params(axis='y', which='major', labelsize=16)


# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='white', edgecolor='black')




textlist=[
r'$\sigma_0=%.2f$'%(s0)+'\n',
r'$\sigma_1=%.2f$'%(s1)+'\n',
r'$\gamma=%.2f$'%(gamma)+'\n',
r'$N=%1d$'%(N)]
textstr=''.join(textlist)

print textstr

#axes[0].text(0.415, 0.2375, textstr, transform=axes[0].transAxes, fontsize=18,verticalalignment='top', bbox=props)
axes[0].text(0.9, 0.95, textstr, transform=axes[0].transAxes, fontsize=18,verticalalignment='top', bbox=props)
"""
Calculate Chi2
"""

def chi2_el(myn,mynu,mymean,mysig):
	return ((truint(myn,mynu,s0,s1)-mymean)/mysig)**2

chi2_f=np.vectorize(chi2_el)

chi2=(1.0/float(len(nu)))*np.sum(chi2_f(N,nu,means,sdevs))
print "CHI2 Ordered: ",chi2

print len(nu)

#props = dict(boxstyle='round', facecolor='blue', edgecolor='black',alpha=0.4)
"""
props = dict(boxstyle='round', facecolor='white', edgecolor='black')
axes[0].text(0.69, 0.96, r'$\chi^2/d.o.f.$'+'={:03.2f}'.format(chi2), transform=axes[0].transAxes, fontsize=18,verticalalignment='top', bbox=props)
"""
#visible_labels = [lab for lab in axes[0].get_yticklabels() if lab.get_visible() is True and lab.get_text() != '']
visible_labels = axes[1].get_yticklabels()
plt.setp(visible_labels[1::2], visible=False)

plt.show()
#print sp.AtmosphericNeutrinoOscillationProbability(1,1,100*param.GeV,param.PI,param)
