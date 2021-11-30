import numpy as np
import scipy.integrate as scintegr
import scipy.interpolate as scinterp
import galario

import os
import sys; sys.path.append('/afs/mpa/data/lucadim/projects/multitool/')
import multitool as mt

def imtouv(uvinp,uvout,pbfile,imfile,Tout=0.00,osz=4,iscompton=False):

  uvout = uvout.replace('.npz','')
  print(uvinp)
  if (not os.path.exists(uvinp)): mt.printError('\nERROR - UVfile not found\n\n')
  if (not os.path.exists(pbfile)): mt.printError('\nERROR - Primay beam .fits file not found\n\n')

  imdata, imhead = mt.reduceStandardFits(imfile)
  
  if not iscompton:
    imfreq = [imhead['FREQ']-0.50*imhead['BAND'],imhead['FREQ']+0.50*imhead['BAND']]
    imdelt = [np.abs(imhead['CDELT1']),np.abs(imhead['CDELT2'])]
    imconv = mt.comptonCorrect(np.array([computeFlatCompton(imfreq,imdelt,order)for order in range(1+osz)]),Tout)
    imdata = imdata/imconv

  pbdata, pbhead = mt.reduceCasaFits(pbfile)
  pbdata[np.where(np.isnan(pbdata))] = 0.0

  uvload = np.load(uvinp,fix_imports=True,encoding='bytes')
  uvdata = np.array([np.copy(uvload[uvload.files[0]][0].flatten()),
                     np.copy(uvload[uvload.files[0]][1].flatten()),
                     np.copy(uvload[uvload.files[0]][2].flatten())])
  uvload.close(); del uvload

  uvimage = np.ascontiguousarray(np.flip(np.multiply(imdata,pbdata),axis=0))
  uvmodel = galario.double.sampleImage(uvimage,pbhead['CDELT2']*np.pi/180.,uvdata[0],uvdata[1])

  uvconv  = mt.comptonCorrect(np.array([mt.comptonToJyPix(uvdata[2],pbhead['CDELT1'],pbhead['CDELT2'])*mt.comptonRelativ(uvdata[2],order) for order in range(1+osz)]),Tout)
  uvmodel = np.multiply(uvconv,uvmodel)
 
  np.savez_compressed(uvout,np.array([uvmodel.real,uvmodel.imag]))



#####################################################################################
# UV model generator - From parameters
#####################################################################################
def partouv(uvinp,uvout,pbfile,uvmod,uvpar,reffreq=1.00E+11,applypb=True):

  uvout = uvout.replace('.npz','')

  if (not os.path.exists(uvinp)): mt.printError('\nERROR - UVfile not found\n\n')
  if (not os.path.exists(pbfile)): mt.printError('\nERROR - Primay beam .fits file not found\n\n')

  pbdata, pbhead = mt.reduceCasaFits(pbfile)
  pbdata[np.where(np.isnan(pbdata))] = 0.0

  uvload = np.load(uvinp,fix_imports=True,encoding='bytes')
  uvdata = np.array([np.copy(uvload[uvload.files[0]][0].flatten()),
                     np.copy(uvload[uvload.files[0]][1].flatten()),
                     np.copy(uvload[uvload.files[0]][2].flatten())])
  uvload.close(); del uvload

  spwpoints = np.zeros(np.shape(uvdata[0]),dtype=np.complex128)
  for m, mod in enumerate(uvpar):
    if (uvmod[m][0]=='pointSource'):
      spwRA   = pbhead['CRPIX1']+(uvpar[m][0]-pbhead['CRVAL1'])*np.cos(np.deg2rad(pbhead['CRVAL2']))/pbhead['CDELT1'] 
      spwDec  = pbhead['CRPIX2']+(uvpar[m][1]-pbhead['CRVAL2'])/pbhead['CDELT2']

      spwairy = scinterp.RectBivariateSpline(np.linspace(0,pbhead['NAXIS2'],pbhead['NAXIS2']),
                                             np.linspace(0,pbhead['NAXIS1'],pbhead['NAXIS1']),
                                             pbdata,kx=1,ky=1,s=0)
      spwbeam = spwairy.ev(spwDec,spwRA)
                
      if uvmod[m][1]=='powerLaw':
        spwpoint = spwbeam*uvpar[m][2]*((uvdata[2]/reffreq)**(uvpar[m][3]))+0j
      elif uvmod[m][1]=='powerLawMod':
        spwpoint = spwbeam*uvpar[m][2]*((uvdata[2]/reffreq)**(uvpar[m][3]+uvpar[m][4]*np.log(uvdata[2]/reffreq)))+0j
      
      spwpoints += galario.double.apply_phase_vis((np.deg2rad(uvpar[m][0]-pbhead['CRVAL1']))*np.cos(np.deg2rad(pbhead['CRVAL2'])),
                                                   np.deg2rad(uvpar[m][1]-pbhead['CRVAL2']),
                                                   uvdata[0], uvdata[1], spwpoint)

  np.savez_compressed(uvout,np.array([spwpoints.real,spwpoints.imag]))




#####################################################################################

def yszCorrect(freq,cdelt,order):
  return mt.comptonToJyPix(freq,cdelt[0],cdelt[1])*mt.comptonRelativ(freq,order)

# Provide SZ relativistic correction terms
def computeFlatCompton(freq,cdelt,order):
  if (freq[0]!=freq[1]):
    return scintegr.quad(yszCorrect,freq[0],freq[1],args=([cdelt[0],cdelt[1]],order))[0]/(freq[1]-freq[0])
  else: return yszCorrect(freq[0],cdelt,order)


#> # Provide SZ relativistic correction terms
#> def computeFlatCompton(freq,cdelt):
#>   if (freq[0]!=freq[1]):
#>     return scintegr.quad(mt.comptonToJyPix,freq[0],freq[1],args=(cdelt[0],cdelt[1]))[0]/(freq[1]-freq[0])
#>   else: return mt.comptonToJyPix(freq[0],cdelt[0],cdelt[1])
