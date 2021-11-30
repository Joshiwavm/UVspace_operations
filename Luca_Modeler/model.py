try:
  from numba import float32, jit
  from numba import float64, vectorize
  nonumba = False
except:
  printInfo('Numba not found: using Numpy vectorization and no JIT compilation.')
  nonumba = True

import scipy.integrate
import numpy as np

# ======================================================================
# Support functions
# ======================================================================
# Line-of-sight eccentricity
# ----------------------------------------------------------------------
def elos(e): return e/np.sqrt(2.00-e*e)

# ======================================================================
# Integrated profiles
# ======================================================================
# When using the gnfw, remove the first point from grid (0.0)
# ----------------------------------------------------------------------
def generateProfile(grid,mguess,mtype,c500=1.00,mass=1.00,limdist=np.inf,epsrel=1.00E-06,freeLS=None):
  if   mtype=='betaPressure': 
    return betaProfile(grid,mguess[6],mguess[2],mguess[3],mguess[4],mguess[8],limdist,epsrel,freeLS)
  elif mtype=='gnfwPressure': 
    return gnfwProfile(grid,mguess[6],mguess[2],mguess[3],mguess[4],mguess[8],mguess[9],mguess[10],limdist,epsrel,freeLS)
  elif mtype=='A10Pressure':
    return a10Profile(grid,mguess[6],mguess[2],mguess[3],mguess[4],mguess[8],mguess[9],mguess[10],mguess[12],c500,mass,limdist,epsrel,freeLS)

# 3D A10 model profile
# ----------------------------------------------------------------------
if nonumba:
  def a10ProfileIntegrand(x,xi,alpha,beta,gamma,ap,c500,mass): 
    return (x**(-gamma))*((1.00+(x**alpha))**((gamma-beta)/alpha))*(x/((x*x-xi*xi)**0.50))*(mass**((ap+0.10)/(1.00+(2.00*x/c500)**3.00)))
else:
  @jit(forceobj=True,cache=True)
  def a10ProfileIntegrand(x,xi,alpha,beta,gamma,ap,c500,mass): 
    return (x**(-gamma))*((1.00+(x**alpha))**((gamma-beta)/alpha))*(x/((x*x-xi*xi)**0.50))*(mass**((ap+0.10)/(1.00+(2.00*x/c500)**3.00)))

# A10 model integrale
# ----------------------------------------------------------------------
if nonumba:
  def _a10ProfileIntegral(x,alpha,beta,gamma,ap,c500,mass,limdist=np.inf,epsrel=1.00E-06): 
    return 2.00*scipy.integrate.quad(a10ProfileIntegrand,x,limdist,args=(x,alpha,beta,gamma,ap,c500,mass),epsrel=epsrel)[0]
  a10ProfileIntegral = np.vectorize(_a10ProfileIntegral)
else:
  @vectorize([float32(float32,float32,float32,float32,float32,float32,float32,float32,float32),
              float64(float64,float64,float64,float64,float64,float64,float64,float64,float64)],forceobj=True)
  def a10ProfileIntegral(x,alpha,beta,gamma,ap,c500,mass,limdist=np.inf,epsrel=1.00E-06): 
    return 2.00*scipy.integrate.quad(a10ProfileIntegrand,x,limdist,args=(x,alpha,beta,gamma,ap,c500,mass),epsrel=epsrel)[0]

# Integrated elliptical A10 model profile
# ----------------------------------------------------------------------
def a10Profile(grid,offset,amp,major,e,alpha,beta,gamma,ap,c500,mass,limdist=np.inf,epsrel=1.00E-06,freeLS=None): 
  integral = np.zeros_like(grid,dtype=np.float64)
  integral[grid<=limdist] = a10ProfileIntegral(grid[grid<=limdist],alpha,beta,gamma,ap,c500,mass,limdist,epsrel)
  ellipse = np.sqrt(1.00-elos(e)**2) if freeLS is None else freeLS
  return offset+amp*major*ellipse*integral

# 3D gNFW model profile
# ----------------------------------------------------------------------
if nonumba:
  def gnfwProfileIntegrand(x,xi,alpha,beta,gamma): 
    return (x**(-gamma))*((1.00+(x**alpha))**((gamma-beta)/alpha))*x/((x*x-xi*xi)**0.50)
else:
  @jit(forceobj=True,cache=True)
  def gnfwProfileIntegrand(x,xi,alpha,beta,gamma): 
    return (x**(-gamma))*((1.00+(x**alpha))**((gamma-beta)/alpha))*x/((x*x-xi*xi)**0.50)

# gNFW model integral
# ----------------------------------------------------------------------
if nonumba:
  def _gnfwProfileIntegral(x,alpha,beta,gamma,limdist=np.inf,epsrel=1.00E-06): 
    return 2.00*scipy.integrate.quad(gnfwProfileIntegrand,x,limdist,args=(x,alpha,beta,gamma),epsrel=epsrel)[0]
  gnfwProfileIntegral = np.vectorize(_gnfwProfileIntegral)
else:
  @vectorize([float32(float32,float32,float32,float32,float32,float32),
              float64(float64,float64,float64,float64,float64,float64)],forceobj=True)
  def gnfwProfileIntegral(x,alpha,beta,gamma,limdist=np.inf,epsrel=1.00E-06): 
    return 2.00*scipy.integrate.quad(gnfwProfileIntegrand,x,limdist,args=(x,alpha,beta,gamma),epsrel=epsrel)[0]

# Integrated elliptical gNFW model profile
# ----------------------------------------------------------------------
def gnfwProfile(grid,offset,amp,major,e,alpha,beta,gamma,limdist=np.inf,epsrel=1.00E-06,freeLS=None): 
  integral = np.zeros_like(grid,dtype=np.float64)
  integral[grid<=limdist] = gnfwProfileIntegral(grid[grid<=limdist],alpha,beta,gamma,limdist,epsrel)
  ellipse = np.sqrt(1.00-elos(e)**2) if freeLS is None else freeLS
  return offset+amp*major*ellipse*integral

# 3D beta model profile
# ----------------------------------------------------------------------
if nonumba:
  def betaProfileIntegrand(x,xi,alpha,beta,gamma): 
    return ((1.00+(x**2.00))**(-1.50*beta))*x/((x*x-xi*xi)**0.50)
else:
  @jit(forceobj=True,cache=True)
  def betaProfileIntegrand(x,xi,alpha,beta,gamma): 
    return ((1.00+(x**2.00))**(-1.50*beta))*x/((x*x-xi*xi)**0.50)

# Beta model integral
# ----------------------------------------------------------------------
if nonumba:
  def _betaProfileIntegral(x,alpha,beta,gamma,limdist=np.inf,epsrel=1.00E-06): 
    return 2.00*scipy.integrate.quad(betaProfileIntegrand,x,limdist,args=(x,beta),epsrel=epsrel)[0]
  betaProfileIntegral = np.vectorize(_betaProfileIntegral)
else:
  @vectorize([float32(float32,float32,float32,float32,float32,float32),
              float64(float64,float64,float64,float64,float64,float64)],forceobj=True)
  def betaProfileIntegral(x,alpha,beta,gamma,limdist=np.inf,epsrel=1.00E-06): 
    return 2.00*scipy.integrate.quad(betaProfileIntegrand,x,limdist,args=(x,beta),epsrel=epsrel)[0]

# Integrated elliptical beta model profile 
# ----------------------------------------------------------------------
def betaProfile(grid,offset,amp,major,e,beta,limdist=np.inf,epsrel=1.00E-06,freeLS=None):
  integral = np.zeros_like(grid,dtype=np.float64)
  integral[grid<=limdist] = betaProfileIntegral(grid[grid<=limdist],beta,limdist,epsrel)
  ellipse = np.sqrt(1.00-elos(e)**2) if freeLS is None else freeLS
  return offset+amp*major*ellipse*integral
