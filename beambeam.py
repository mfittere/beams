from constants import *
from beamparam import gammarel
import numpy as np
                       
#def gammarel(E,m0=mp):
#  """returns the relativistic gamma E [GeV]"""
#  return E*1.e3/m0

def bbparam(Np,beta0x,beta0y,epsnx,epsny,r0=rp):
  """returns the incoherent beam-beam parameter (as default for protons)
  assuming the same beam parameters for both colliding beams.
  Np=particles per bunch, rp=classical proton radius,
  beta0x/y=hor./vert. beta function at IP [m],
  epsnx/y=hor/vert. normalized emittance [mum]"""
  xix=Np*r0*beta0x/(2*np.pi*np.sqrt(beta0x*epsnx)*(np.sqrt(beta0x*epsnx)+np.sqrt(beta0y*epsny))*1.e-6) 
  xiy=Np*r0*beta0y/(2*np.pi*np.sqrt(beta0y*epsny)*(np.sqrt(beta0x*epsnx)+np.sqrt(beta0y*epsny))*1.e-6)
  return (xix,xiy)

def reductionlumioffset(dx,sigx,dy,sigy,sigz,phi):
  """reduction factor R due to offset collision
       L=R*L_0
  dx,sigx,dy,sigy in mum, phi in murad,sigz in cm
  dx, dy = beam separation = x(b1)-x(b2)
  sigx,sigy = beam size at IP, sigz=bunchlength
  phi = half crossing angle"""
  return np.exp(-((dx*1.e-6/2)**2/((sigx*1.e-6)**2*(np.cos(phi*1.e-6))**2+(sigz*1.e-2)**2*(np.sin(phi*1.e-6))**2))-(dy*1.e-6/(2*sigy*1.e-6))**2)

