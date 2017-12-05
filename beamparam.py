from constants import *
import numpy as np

def gammarel(EGeV,m0=mp):
  """returns the relativistic gamma
  Parameters:
  -----------
  EGeV: kinetic energy [GeV]
  m0: restmass [MeV]
  """
  return (EGeV*1.e3+m0)/m0

def betarel(EGeV,m0=mp):
  g=gammarel(EGeV,m0=m0)
  return np.sqrt(1-1/g**2)

def brho(EGeV,m0=mp):
  """returns the magnetic rigidity [T/m]
  Parameters:
  -----------
  EGeV: kinetic energy [GeV]
  m0: restmass [MeV]
  """
  beta=betarel(EGeV,m0)
  return 10/2.998*beta*EGeV
  
def frev(EGeV,C=26658.8832,m0=mp):
  """calculate the revolution frequency"""
  b=betarel(EGeV=EGeV,m0=m0)
  return (b*clight)/C

def emitnorm(eps,EGeV,m0=mp):
  """returns normalized emittance in [mum].
  input: eps [mum], E [GeV], m0 [MeV]
  """
  gamma= gammarel(EGeV,m0)
  beta = betarel(EGeV,m0)
  return eps*(gamma*beta)

def emitrms(epsn,EGeV,m0=mp):
  """returns rms emittance in [mum].
  input: epsn [mum], E [GeV], m0 [MeV]
  """
  gamma= gammarel(EGeV,m0)
  beta = betarel(EGeV,m0)
  return epsn/(gamma*beta)

def emitnorm(eps,EGeV,m0=mp):
  """returns normalized emittance in [mum].
  input: eps [mum], E [GeV], m0 [MeV]
  """
  gamma= gammarel(EGeV,m0)
  beta = betarel(EGeV,m0)
  return eps*(gamma*beta)

def eta(alpha_c,EGeV,m0=mp):
  """calculates the phase slip factor
  and gamma transition where alpha_c
  is the momentum compaction factor"""
  g=gammarel(EGeV,m0=m0)
  eta=alpha_c-1/g**2
  gt=np.sqrt(1/np.abs(alpha_c))
  return eta,gt

def nus(h,eta,EGeV,VMV,phirf=0,m0=mp):
  """calculate synchrotron tune
  h: harmonic number
  eta: phase slip factor
  EGeV: beam energy [GeV]
  VMV: rf voltage [MV]
  phirf: rf phase [rad]
  m0: particle mass [MeV]
  """
  return np.sqrt(h*np.abs(eta)*VMV*1.e6/(2*np.pi*betarel(EGeV,m0=m0)**2*EGeV*1.e9))

def bucket_area(h,eta,EGeV,VMV,frf,m0=mp):
  """returns the bucket area [eVs]
  h: harmonic number
  eta: phase slip factor
  EGeV: beam energy [GeV]
  VMV: rf voltage [MV]
  frf: rf frequency [Hz]
  m0: particle mass [MeV]
  """
  b=betarel(EGeV=EGeV,m0=m0)
  return (8*b)/(frf*np.pi)*np.sqrt((EGeV*1.e9*VMV*1.e6)/(2*np.pi*h*eta))

def bucket_height(h,eta,EGeV,VMV,m0=mp):
  """returns the bucket half heigth dE/E
  full bucket height reaches from [-dE/E,+dE/E]
  h: harmonic number
  eta: phase slip factor
  EGeV: beam energy [GeV]
  VMV: rf voltage [MV]
  m0: particle mass [MeV]
  """
  b=betarel(EGeV=EGeV,m0=m0)
  return b*np.sqrt((2*VMV*1.e6)/(np.pi*h*eta*EGeV*1.e9))

def bucket_length(frf):
  """returns the bucket length l=clight/frf [m]
  frf: rf frequency [Hz]
  """
  return clight/frf

def beta_z(frev,eta,nus):
  """returns the longitudinal beta function [m]
  where beta_z=sigz/sigp with bunch length sigz
  and momentum spread sigp.
  Note that dp/p=1/beta**2 dE/E where beta is
  the relativistic beta.

  Parameters:
  -----------
  frev: revolution frequency [Hz]
  eta: phase slip factor
  nus: synchrotron tune
  """
  return eta*clight/(2*np.pi*nus*frev)

def stored_beam_energy(EGeV,nb=2808,np=1.15e11):
  """returns the beam enery [MJ]
  Parameters:
  -----------
  EGeV: beam energy
  nb: number of bunches
  np: number of particles per bunch
  """
  return EGeV*1.e9*nb*np*echarge*1.e-6

def beam_current(nb,np,frev=11245):
  """returns the beam current [A], assumes
  that particles have charge +/-e.
  Parameters:
  ----------
  nb: number of bunches
  np: number of particles per bunch
  frev: revolution frequency [Hz]
  """
  return echarge*nb*np*frev
