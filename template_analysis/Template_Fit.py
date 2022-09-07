'''
Created on Sep 6, 2022

@author: Matthias Plum
'''
import probfit
import iminuit

class CR_Template_Analysis:


    def __init__(self,minos=False):
      self.minos = minos
      self.num_pdfs = 1
      self.pdf = None
      self.likelihood = None
      self.minuit = None

    def join_pdfs(self,template_pdfs):
      ### Set number of PDFs in template analysis
      self.num_pdfs = len(template_pdfs)
      ### Convert PDFs to extended PDFs
      pdfs_list = []
      for i, pdf in enumerate(template_pdfs):
        pdfs_list.append(probfit.Extended(pdf, extname='N'+str(i+1)))

      ### Define an extended PDF consisting of four components
      self.pdf = probfit.AddPdf(*pdfs_list)
    
    def template_binned_likelihood(self, data,set_bins,set_fitrange):
      ### Start fit-parameter
      ### Assume worst case
      pars = {}
      for i in range(self.num_pdfs):
        pars['N'+str(i+1)]=len(data)/self.num_pdfs

      self.likelihood = probfit.BinnedLH(self.pdf, data, bins=len(set_bins), extended=True, bound=set_fitrange)

      self.minuit = iminuit.Minuit(self.likelihood, pedantic=False, print_level=0,**pars)
      self.minuit.migrad()
      self.minuit.hesse()
      if self.minos:
        self.minuit.minos()
        
      print('fmin:')
      print(self.minuit.fmin)
      print('error matrix:')
      print(self.minuit.matrix())
      # or the correlation matrix
      print('correlation matrix:')
      print(self.minuit.matrix(correlation=True))

    def template_unbinned_likelihood(self, data):
      ### Start fit-parameter
      ### Assume worst case
      pars = {}
      for i in range(self.num_pdfs):
        pars['N'+str(i+1)]=len(data)/self.num_pdfs

      self.likelihood = probfit.UnbinnedLH(self.pdf, data, extended=True)

      self.minuit = iminuit.Minuit(self.likelihood, pedantic=False, print_level=0,**pars)
      self.minuit.migrad()
      self.minuit.hesse()
      if self.minos: 
        self.minuit.minos()
        
      print('fmin:')
      print(self.minuit.fmin)
      print('error matrix:')
      print(self.minuit.matrix())
      # or the correlation matrix
      print('correlation matrix:')
      print(self.minuit.matrix(correlation=True))