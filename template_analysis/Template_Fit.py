'''
Created on Sep 6, 2022

@author: Matthias Plum
@email:  matthias.plum@sdsmt.edu
'''
import probfit
import iminuit
import pkg_resources

class Template_Analysis:

  def __init__(self, minos=False, binned=False, strategy=0):
    if pkg_resources.parse_version(iminuit.__version__) < pkg_resources.parse_version("2.16.0"):
      raise RuntimeError("Found iminiut version "+iminuit.__version__+" is to old. Please upgrade.")
    
    self.minos      = minos
    self.binned     = binned
    self.strategy   = strategy
    self.num_pdfs   = 1
    self.pdf        = None
    self.likelihood = None
    self.minuit     = None

  def join_pdfs(self,template_pdfs):
    ### Set number of PDFs in template analysis
    self.num_pdfs = len(template_pdfs)
    ### Convert PDFs to extended PDFs
    pdfs_list = []
    for i, pdf in enumerate(template_pdfs):
      pdfs_list.append(probfit.Extended(pdf, extname='N'+str(i+1)))

    ### Define an extended PDF consisting of four components
    self.pdf = probfit.AddPdf(*pdfs_list)

  def template_likelihood(self, data, set_bins, set_fitrange, weights=None):
    ### Start fit-parameter (Assume worst case)
    pars = {}
    for i in range(self.num_pdfs):
      pars['N'+str(i+1)]=len(data)/self.num_pdfs
      
    if self.binned:
      self.likelihood = probfit.BinnedLH(self.pdf, data, bins=len(set_bins), extended=True, bound=set_fitrange, weights=weights)
    else:
      self.likelihood = probfit.UnbinnedLH(self.pdf, data, weights=weights, extended=True)

    self.minuit = iminuit.Minuit(self.likelihood, **pars)
    self.minuit.limits["N1", "N2", "N3",'N4'] = (0, len(data))
    self.minuit.strategy = self.strategy
    self.minuit.errordef = 0.5
    self.minuit.print_level = 0
    
    self.minuit.migrad()
    self.minuit.hesse()
    if self.minos:
      self.minuit.minos()
      
    print('fmin:')
    print(self.minuit.fmin)
    print('error matrix:')
    print(self.minuit.covariance)
    ### or the correlation matrix
    print('correlation matrix:')
    print(self.minuit.covariance.correlation())
    if self.minos:
      print(self.minuit.merrors)
