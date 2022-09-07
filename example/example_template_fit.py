'''
Created on Sep 6, 2022

@author: Matthias Plum
@email:  matthias.plum@sdsmt.edu
'''
import numpy as np
from scipy.stats import truncnorm

import pylab as mpl

from template_analysis.Template_Fit import Template_Analysis
####Setup example 
x = np.linspace(0,np.log(56),200)
bins = np.linspace(0,np.log(56),50)
fit_range=(0,np.log(56))

binned = True
run_minos = False
seed = 42
rng = np.random.RandomState(seed)
### Number of events per elementary group
nH  = 500
nHe = 100
nO  = 200
nFe = 500

### Template shape generation
bound_a   = 0.
bound_b   = np.log(56)

mean_H    = 0
scale_H   = 0.75
mean_He    = np.log(4)
scale_He   = .75
mean_O    = np.log(14)
scale_O   = .75
mean_Fe    = np.log(56)
scale_Fe   = 0.75

H_a, H_b = (bound_a - mean_H) / scale_H, (bound_b - mean_H) / scale_H
He_a, He_b = (bound_a - mean_He) / scale_He, (bound_b - mean_He) / scale_He
O_a, O_b = (bound_a - mean_O) / scale_O, (bound_b - mean_O) / scale_O
Fe_a, Fe_b = (bound_a - mean_Fe) / scale_Fe, (bound_b - mean_Fe) / scale_Fe

rv_H = truncnorm(H_a, H_b,loc=mean_H, scale=scale_H)
rv_He = truncnorm(He_a, He_b,loc=mean_He, scale=scale_He)
rv_O = truncnorm(O_a, O_b,loc=mean_O, scale=scale_O)
rv_Fe = truncnorm(Fe_a, Fe_b,loc=mean_Fe, scale=scale_Fe)

### Sampling data set
H_data=rv_H.rvs(nH,random_state=seed)
He_data=rv_He.rvs(nHe,random_state=seed)
O_data=rv_O.rvs(nO,random_state=seed)
Fe_data=rv_Fe.rvs(nFe,random_state=seed)

data = [H_data,He_data,O_data,Fe_data]

data_flat = np.concatenate(data)

###Plotting PDFs
fig1, ax = mpl.subplots(1, 1)
ax.plot(x, rv_H.pdf(x),'r-', lw=2, alpha=0.6, label='H pdf')
ax.plot(x, rv_He.pdf(x),'y-', lw=2, alpha=0.6, label='He pdf')
ax.plot(x, rv_O.pdf(x),'g-', lw=2, alpha=0.6, label='O pdf')
ax.plot(x, rv_Fe.pdf(x),'b-', lw=2, alpha=0.6, label='Fe pdf')
ax.legend(loc=0)

### Ploting histogram and create subplot for fit results
fig2, (ax1,ax2) = mpl.subplots(2, 1)

ax1.hist(data, bins=bins, density=False, stacked=True,histtype='step', color=['r','orange','g','b'])
ax2.hist(data_flat, bins=bins,density=True)

### Create list of the template PDFs or functions
template_pdfs = [rv_H.pdf,rv_He.pdf,rv_O.pdf,rv_Fe.pdf]

### Run template fitting method binned or unbinned
template = Template_Analysis(minos=run_minos)
template.join_pdfs(template_pdfs)
if binned:
  template.template_binned_likelihood(data_flat, bins, fit_range)
  template.likelihood.show(template.minuit, ax=ax2, parts=True);
else:
  template.template_unbinned_likelihood(data_flat)
  template.likelihood.show(template.minuit, ax=ax2, bins=len(bins), parts=True);

mpl.show()

