import numpy as np

from astropy import units as u
import astropy.constants as const
from astropy.cosmology import FlatLambdaCDM

class getinfo:
      def __init__(self,
                   reffreq=1e11,
                   parline=[-5,5,100],
                   limdist=np.inf,
                   limepsr=1.00E-06,
                   fb=0.175,
                   mu=0.590,
                   mue=1.140,
                   cosmo=None):

        self.reffreq = reffreq
        self.limdist = limdist
        self.limepsr = limepsr

        self.linmesh = np.array([])

        if len(parline):
            if len(parline)!=3: 
                printError('Provide correct numbers of line parameters [log scale]')
            
        self.linmesh = np.append(0.0,np.logspace(parline[0],parline[1],parline[2],dtype=np.float64))

        self.ysznorm = const.sigma_T/const.m_e/const.c**2
        self.ysznorm = self.ysznorm.to(u.cm**3/u.keV/u.Mpc)

        self.cosmo = FlatLambdaCDM(H0=70.00,Om0=0.30) if cosmo is None else cosmo

        self.fb = fb
        self.mu = mu
        self.mue = mue
        
        
    
def yszCorrect(freq,cdelt,order):
    return mt.comptonToJyPix(freq,cdelt[0],cdelt[1])*mt.comptonRelativ(freq,order)

# Provide SZ relativistic correction terms
def computeFlatCompton(freq,cdelt,order):
    if (freq[0]!=freq[1]):
        return scintegr.quad(yszCorrect,freq[0],freq[1],args=([cdelt[0],cdelt[1]],order))[0]/(freq[1]-freq[0])
    else: 
        return yszCorrect(freq[0],cdelt,order)
