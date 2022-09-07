'''
Created on Sep 6, 2022

@author: Matthias Plum
'''
import probfit
import iminuit

class CR_Template_Analysis:


    def __init__(self,minos=False):
      self.minos = minos
      self.pdf = None
      self.likelihood = None
      self.minuit = None

    def join_pdfs(self,template_pdfs):
      ### Convert PDFs to extended PDFs
      ext_pdf1 = probfit.Extended(template_pdfs[0], extname='N1')
      ext_pdf2 = probfit.Extended(template_pdfs[1], extname='N2')
      ext_pdf3 = probfit.Extended(template_pdfs[2], extname='N3')
      ext_pdf4 = probfit.Extended(template_pdfs[3], extname='N4')
      
      ### Define an extended PDF consisting of four components
      self.pdf = probfit.AddPdf(ext_pdf1,ext_pdf2,ext_pdf3,ext_pdf4)
    
    def template_binned_likelihood(self, data,set_bins,set_fitrange):
      ### Start fit-parameter
      ### Assume worst case
      pars = dict(N1=len(data)/4,
                  N2=len(data)/4,
                  N3=len(data)/4,
                  N4=len(data)/4,
                  )

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
      pars = dict(N1=len(data)/4,
                  N2=len(data)/4,
                  N3=len(data)/4,
                  N4=len(data)/4,
                  )

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