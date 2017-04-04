"""python module to calculate e-lens parameters"""
from constants import *
from beamparam import *
import numpy as np

def hel_thetamax(r2,I,l,EkeV_e,EGeV_p,direction='opposite'):
  """
  calculate maximum kick from ideal e-lens:
  
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
