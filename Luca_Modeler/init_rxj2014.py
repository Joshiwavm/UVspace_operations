import numpy as np

from astropy.io import fits

import glob

import tempfile; tempfile.tempdir = '/afs/mpa/temp/lucadim/tmp'
import imp

import os
import sys
sys.path.append('/afs/mpa/data/lucadim/projects/MODELER/v0.0.4/')
sys.path.append('/afs/mpa/data/lucadim/projects/multitool/')

steps = [3,4]

cwddir = os.getcwd()
os.chdir('/afs/mpa/temp/lucadim/alma/RXJ2014/')

model = '02'
array = 'c07m'
if   array=='c12m': imsize, cell, rawspw =  512, '0.30arcsec', '0,1,2'
elif array=='c07m': imsize, cell, rawspw =  256, '1.50arcsec', '1,2,3'

spws = [0,1,2]
fields = [0]

def prunebeam(imname):
  exportfits('{0}.pb'.format(imname),'{0}.pbeam.fits'.format(imname),overwrite=True)
  exportfits('{0}.psf'.format(imname),'{0}.psf.fits'.format(imname),overwrite=True)
  exportfits('{0}.image'.format(imname),'{0}.image.fits'.format(imname),overwrite=True)

  # os.system('rm -rf {0}.pb'.format(imname))
  # os.system('rm -rf {0}.psf'.format(imname))
  # os.system('rm -rf {0}.image'.format(imname))
  os.system('rm -rf {0}.model'.format(imname))
  os.system('rm -rf {0}.sumwt'.format(imname))
  os.system('rm -rf {0}.weight'.format(imname))
  os.system('rm -rf {0}.residual'.format(imname))

  hdu = fits.open('{0}.pbeam.fits'.format(imname))[0]
  hdu.data[np.isnan(hdu.data)] = 0.0

  data = np.copy(hdu.data[0,0])
  diff0a = np.ones(np.shape(data))
  diff0b = np.ones(np.shape(data))
  diff1a = np.ones(np.shape(data))
  diff1b = np.ones(np.shape(data))

  diff0a[1:,   :  ] = np.diff(data,axis=0); diff0a[int(diff0a.shape[0]/2)+1:,:] *= -1.0 
  diff0b[ :-1, :  ] = np.diff(data,axis=0); diff0b[int(diff0b.shape[0]/2):,:] *= -1.0 
  diff1a[ :  ,1:  ] = np.diff(data,axis=1); diff1a[:,int(diff1a.shape[1]/2)+1:] *= -1.0 
  diff1b[ :  , :-1] = np.diff(data,axis=1); diff1b[:,int(diff1b.shape[1]/2):] *= -1.0 

  mask0 = np.logical_or(diff0a<0,diff0b<0)
  mask1 = np.logical_or(diff1a<0,diff1b<0)
  mask = np.logical_or(mask0,mask1)
  mask = np.logical_or(mask,data==0)
  data[mask] = np.nan

  hdu.data[0,0] = np.copy(data)
  hdu.writeto('{0}.pbeam.fits'.format(imname),overwrite=True)



####################################################
# Saving visibility points
####################################################
if (1 in steps):
  import msmanager
  imp.reload(msmanager)

  rawname = '{0}/calibrated/rxj2014_calibrated_{0}.bin.ms'.format(array)

  if not os.path.exists(rawname):
    split(vis='{0}/calibrated/rxj2014_calibrated_{0}.raw.ms'.format(array),outputvis=rawname,keepflags=False,datacolumn='data',spw=rawspw)#,timebin='30s')

  visname = '{0}/calibrated/rxj2014_calibrated_{0}.sub.{1}.ms'.format(array,model)
  if not os.path.exists(visname): split(vis=rawname,outputvis=visname,keepflags=False,datacolumn='data')

  visname = '{0}/calibrated/rxj2014_calibrated_{0}.mod.{1}.ms'.format(array,model)
  if not os.path.exists(visname): split(vis=rawname,outputvis=visname,keepflags=False,datacolumn='data')

  if not os.path.exists('{0}/models/numpy/uv/'.format(array,model)):
    os.system('mkdir -p {0}/models/numpy/uv/'.format(array,model))

  npzname = '{0}/models/numpy/uv/rxj2014_calibrated_{0}.mod.{1}.field-fid.spw-sid.uv.npz'.format(array,model)
  msmanager.uvsaver(visname=visname,npzname=npzname,spws=spws,fields=fields)


####################################################
# Generate model visibilities
####################################################
if (2 in steps):
  import modeler
  imp.reload(modeler)

  import dill
  import dynesty.utils
  import scipy.special

  uvsamp = '/afs/mpa/data/lucadim/projects/eszee/v0.0.3/rxj2014/rxj2014_{0}_a10_up_2000_stat_auto/rxj2014_{0}_a10_up_2000_stat_auto_pickle'.format(model)
  with open(uvsamp,'rb') as uvopen: uvload = dill.load(uvopen)
  uvwght = uvload['logwt']-scipy.special.logsumexp(uvload['logwt']-uvload['logz'][-1])
  uvmod = []; uvwght = np.exp(uvwght-uvload['logz'][-1])

  for p in range(uvload.samples.shape[-1]):
    uvmod.append(dynesty.utils.quantile(uvload.samples[:,p],[0.16,0.50,0.84],weights=uvwght))
  uvmod = np.transpose(np.asarray(uvmod))
  
  uvpar = np.array([[uvmod[1][4],uvmod[1][5],uvmod[1][6],uvmod[1][7]]])

  for s, spw in enumerate(spws):
    for f, field in enumerate(fields):

      uvinp = '{0}/models/numpy/uv/rxj2014_calibrated_{0}.mod.{3}.field-{1}.spw-{2}.uv.npz'.format(array,field,spw,model)
      uvout = '{0}/models/numpy/uv/rxj2014_calibrated_{0}.mod.{3}.field-{1}.spw-{2}.ri.npz'.format(array,field,spw,model)
      pbfile = '{0}/splits/rxj2014_calibrated_{0}.bin.im.field-{1}.spw-{2}.pbeam.fits'.format(array,field,spw)

      modeler.partouv(uvinp=uvinp,uvout=uvout,pbfile=pbfile,uvmod=[['pointSource','powerLaw'] for m, mod in enumerate(uvpar)],uvpar=uvpar,reffreq=9.20E+10,applypb=True)


####################################################
# Load model into a MS file and compute residuals
####################################################
if (3 in steps):
  import msmanager
  imp.reload(msmanager)

  visname = '{0}/calibrated/rxj2014_calibrated_{0}.sub.{1}.ms'.format(array,model)
  npzname = '{0}/models/numpy/uv/rxj2014_calibrated_{0}.mod.{1}.field-fid.spw-sid.ri.npz'.format(array,model)
        
  msmanager.uvloader(visname=visname,npzname=npzname,todo='subtract',spws=spws,fields=fields)

  visname = '{0}/calibrated/rxj2014_calibrated_{0}.mod.{1}.ms'.format(array,model)
  npzname = '{0}/models/numpy/uv/rxj2014_calibrated_{0}.mod.{1}.field-fid.spw-sid.ri.npz'.format(array,model)
  msmanager.uvloader(visname=visname,npzname=npzname,todo='replace',spws=spws,fields=fields)



####################################################
# Additional operations
####################################################
if (4 in steps):
  for s, spw in enumerate(spws):
    for f, field in enumerate(fields):
      if not os.path.exists('{0}/models/rxj2014_calibrated_{0}.mod.012x.{3}.field-{1}.spw-{2}.im.m100.image'.format(array,field,spw,model)):
        tclean(vis        = '{0}/calibrated/rxj2014_calibrated_{0}.mod.012x.{1}.ms'.format(array,model),
               imagename  = '{0}/models/rxj2014_calibrated_{0}.mod.012x.{3}.field-{1}.spw-{2}.im.m100'.format(array,field,spw,model),
               spw        = '{0}'.format(spw),
               datacolumn =     'data',
               imsize     =     imsize, 
               cell       =       cell,
               gridder    = 'standard',
               specmode   =      'mfs',
               pblimit    =       0.00,
               weighting  =   'briggs',
               robust     =      -1.00,
               niter      =          0)

      prunebeam('{0}/models/rxj2014_calibrated_{0}.mod.012x.{3}.field-{1}.spw-{2}.im.m100'.format(array,field,spw,model))

      if not os.path.exists('{0}/models/rxj2014_calibrated_{0}.raw.012x.field-{1}.spw-{2}.im.m100.image'.format(array,field,spw,model)):
        tclean(vis        = '{0}/calibrated/rxj2014_calibrated_{0}.raw.012x.ms'.format(array,model),
               imagename  = '{0}/models/rxj2014_calibrated_{0}.raw.012x.field-{1}.spw-{2}.im.m100'.format(array,field,spw,model),
               spw        = '{0}'.format(spw),
               datacolumn =     'data',
               imsize     =     imsize, 
               cell       =       cell,
               gridder    = 'standard',
               specmode   =      'mfs',
               pblimit    =       0.00,
               weighting  =   'briggs',
               robust     =      -1.00,
               niter      =          0)

      prunebeam('{0}/models/rxj2014_calibrated_{0}.raw.012x.field-{1}.spw-{2}.im.m100'.format(array,field,spw,model))

      imagename = ['{0}/models/rxj2014_calibrated_{0}.raw.012x.field-{1}.spw-{2}.im.m100.image'.format(array,field,spw,model),
                   '{0}/models/rxj2014_calibrated_{0}.mod.012x.{3}.field-{1}.spw-{2}.im.m100.image'.format(array,field,spw,model)]
      outfile = '{0}/models/rxj2014_calibrated_{0}.res.012x.{3}.field-{1}.spw-{2}.im.m100.image'.format(array,field,spw,model)

      immath(imagename=imagename,outfile=outfile,expr='IM0-IM1')
      exportfits(outfile,'{0}.fits'.format(outfile),overwrite=True)

####################################################
# Additional operations
####################################################
if (5 in steps):
  if not os.path.exists('{0}/models/rxj2014_calibrated_{0}.mod.{3}.im.p200.image'.format(array,0,0,model)):
    tclean(vis        = '{0}/calibrated/rxj2014_calibrated_{0}.mod.{1}.ms'.format(array,model),
           imagename  = '{0}/models/rxj2014_calibrated_{0}.mod.{3}.im.p200'.format(array,0,0,model),
           spw        =    '0,1,2',
           datacolumn =     'data',
           imsize     =     imsize, 
           cell       =       cell,
           gridder    = 'standard',
           specmode   =      'mfs',
           pblimit    =       0.00,
           weighting  =   'briggs',
           robust     =       2.00,
           niter      =          0)

  exportfits('{0}/models/rxj2014_calibrated_{0}.mod.{3}.im.p200.image'.format(array,0,0,model),
             '{0}/models/rxj2014_calibrated_{0}.mod.{3}.im.p200.image.fits'.format(array,0,0,model),
             overwrite=True)

  if not os.path.exists('{0}/models/rxj2014_calibrated_{0}.sub.{3}.im.p200.image'.format(array,0,0,model)):
    tclean(vis        = '{0}/calibrated/rxj2014_calibrated_{0}.sub.{1}.ms'.format(array,model),
           imagename  = '{0}/models/rxj2014_calibrated_{0}.sub.{3}.im.p200'.format(array,0,0,model),
           spw        =    '0,1,2',
           datacolumn =     'data',
           imsize     =     imsize, 
           cell       =       cell,
           gridder    = 'standard',
           specmode   =      'mfs',
           pblimit    =       0.00,
           weighting  =   'briggs',
           robust     =       2.00,
           niter      =          0)

  exportfits('{0}/models/rxj2014_calibrated_{0}.sub.{3}.im.p200.image'.format(array,0,0,model),
             '{0}/models/rxj2014_calibrated_{0}.sub.{3}.im.p200.image.fits'.format(array,0,0,model),
             overwrite=True)

  if not os.path.exists('{0}/models/rxj2014_calibrated_{0}.bin.im.p200.image'.format(array,0,0,model)):
    tclean(vis        = '{0}/calibrated/rxj2014_calibrated_{0}.bin.ms'.format(array,model),
           imagename  = '{0}/models/rxj2014_calibrated_{0}.bin.im.p200'.format(array,0,0,model),
           spw        =    '0,1,2',
           datacolumn =     'data',
           imsize     =     imsize, 
           cell       =       cell,
           gridder    = 'standard',
           specmode   =      'mfs',
           pblimit    =       0.00,
           weighting  =   'briggs',
           robust     =       2.00,
           niter      =          0)

  exportfits('{0}/models/rxj2014_calibrated_{0}.bin.im.p200.image'.format(array,0,0,model),
             '{0}/models/rxj2014_calibrated_{0}.bin.im.p200.image.fits'.format(array,0,0,model),
             overwrite=True)

  imagename = ['{0}/models/rxj2014_calibrated_{0}.bin.im.p200.image'.format(array,0,0,model),
               '{0}/models/rxj2014_calibrated_{0}.mod.{3}.im.p200.image'.format(array,0,0,model)]
  outfile = '{0}/models/rxj2014_calibrated_{0}.res.{3}.im.p200.image'.format(array,0,0,model)

  immath(imagename=imagename,outfile=outfile,expr='IM0-IM1')
  exportfits(outfile,'{0}.fits'.format(outfile),overwrite=True)