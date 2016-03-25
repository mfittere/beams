from constants import *
import numpy as np

def gammarel(EGeV,m0=mp):
  """returns the relativistic gamma
   input: energy E [GeV], m0 [MeV]"""
  return EGeV*1.e3/m0

def betarel(EGeV,m0=mp):
  g=gammarel(EGeV,m0=m0)
  return np.sqrt(1-1/g**2)


def emitrms(epsn,EGeV,m0=mp):
  """returns rms emittance in [mum].
  input: epsn [mum], E [GeV], m0 [MeV]
  """
  gamma= gammarel(EGeV,m0)
  beta = betarel(EGeV,m0)
  return epsn/(gamma*beta)

def eta(alpha_c,EGeV,m0=mp):
  """calculates the phase slip factor
  and gamma transition where alpha_c
  is the momentum compaction factor"""
  g=gammarel(EGeV,m0=m0)
  eta=alpha_c-1/g**2
  gt=np.sqrt(1/np.abs(alpha_c))
  return eta,gt

def frev(EGeV,C=26658.8832,m0=mp):
  """calculate the harmonic number"""
  b=betarel(EGeV=EGeV,m0=m0)
  return (b*clight)/C

def nus(h,eta,EGeV,VMV,phirf=0,m0=mp):
  """calculate synchrotron tune
  h: harmonic number
  eta: phase slip factor
  EGeV: beam energy [GeV]
  VMV: rf voltage in [MV]
  phirf: rf phase [rad]
  m0: particle mass [MeV]
  """
  return np.sqrt(h*np.abs(eta)*VMV*1.e6/(2*np.pi*betarel(EGeV,m0=m0)**2*EGeV*1.e9))
