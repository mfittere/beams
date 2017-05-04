"""python module to calculate e-lens parameters"""
from constants import *
from beamparam import *
import numpy as np

def hel_thetamax(r2,I,l,EkeV_e,EGeV_p,direction='opposite'):
  """
  calculate maximum kick from uniform distribution
  in r=sqrt(x**2+y**2) for e-lens:
  
  Parameters:
  -----------
  r2     : outer radius [m]
  I      : current of e-beam [A]
  l      : length e-lens [m]
  EkeV_e : e-beam energy [keV]
  EGeV_p : p-beam energy [keV]
  direction: direction of e-beam and p-beam,
    options are 'same' and 'opposite'.
    For 'same' (v_e*v_p > 0) the magnetic and electric force
    point in opposite directions, for 'opposite' (v_e*v_p < 0)
    the magnetic and electric force add up
  
  Returns:
  --------
  thetamax: maximum kick [rad]
  """
  beta_e=betarel(EkeV_e*1.e-6, m0=me)
  beta_p=betarel(EGeV_p, m0=mp)
  brho_p=brho(EGeV_p, m0=mp)
  print beta_e,beta_p,brho_p
  if (direction=='opposite'):
    return 2*l*I*(1+beta_e*beta_p)/(4*np.pi*eps0_vacuum*brho_p*beta_e
                                    *beta_p*clight**2)*1/r2
  if (direction=='same'):
    return 2*l*I*(1-beta_e*beta_p)/(4*np.pi*eps0_vacuum*brho_p*beta_e
                                    *beta_p*clight**2)*1/r2

def hel_kick(r,r1,r2,thetamax):
  """
  return the kick from uniform distribution
  in r=sqrt(x**2+y**2) for e-lens:

  Parameters:
  -----------
  r        : proton particle position r=sqrt(x**2+y**2)
  r1       : inner radius [m]
  r2       : outer radius [m]
  thetamax : maximum kick [m] (see hel_thetamax)
  """
  # project -r on r
  r = np.abs(r)
  # r > r2 -> kick =1
  kick = np.where(r>=r2,1,0)
  # overwrite intermediate region
  kick_aux = np.where((r>r1) & (r<r2),(r**2-r1**2)/(r2**2-r1**2),0)
  return (kick+kick_aux)*(r2/r)*thetamax

