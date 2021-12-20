from src.Utillities import *

from astropy.constants import c
from astropy.io import fits
import sys
import os

import casatools
from casatasks import *

class MsReader:
    def __init__(self, filename, obs_type, spws = ['0','1','2', '3'], fields = ['0'], band = 'band3'):
        self.ms_file  = filename
        self.obs_type = obs_type
        self.band     = band
        
        self._set_obs(spws, fields)
        
        self.binvis = './output/{1}/output_{0}_{1}.ms.field-fid.spw-sid'.format(band,self.obs_type)
        self.timebin = '0s'
        
        
    def _set_obs(self,  spws, fields):
        self.fields = fields
        self.spws = spws
        
        if self.obs_type=='com07m':
            self.imsize = 256
            self.imcell = '1.50arcsec'

        elif self.obs_type=='com12m':
            self.imsize = 1024
            self.imcell = '0.15arcsec'
            
        else:
            raise print("Wrong array element")
        
    def prep_pb(self):
        
        ms=casatools.ms()
        self.splitms_file = []
        
        for f, field in enumerate(self.fields):
            print('Processing field {0}'.format(field))
            for s, spw in enumerate(self.spws):
                print('- Spectral window {0}'.format(spw))

                outvis = self.binvis.replace('-fid','-{0}'.format(field))
                outvis = outvis.replace('-sid','-{0}'.format(spw))
                    
                tclean(vis         =                  self.ms_file,
                       imagename   =   outvis.replace('.ms','.im'),
                       datacolumn  =                        'data',
                       field       =                         field,
                       spw         =                           spw,
                       niter       =                             0,
                       pblimit     =                           0.0,
                       imsize      =                   self.imsize,
                       cell        =                   self.imcell,
                       gridder     =                    'standard',
                       weighting   =                     'natural',
                       specmode    =                         'mfs',
                       parallel    =                         False)

                exportfits('{0}.pb'.format(outvis.replace('.ms','.im')),
                           '{0}.pbeam.fits'.format(outvis.replace('.ms','.im')),
                           overwrite=True)
            
                os.system('rm -rf {0}.pb'.format(outvis.replace('.ms','.im')))
                os.system('rm -rf {0}.psf'.format(outvis.replace('.ms','.im')))

                os.system('rm -rf {0}.image'.format(outvis.replace('.ms','.im')))
                os.system('rm -rf {0}.model'.format(outvis.replace('.ms','.im')))
                os.system('rm -rf {0}.sumwt'.format(outvis.replace('.ms','.im')))
                os.system('rm -rf {0}.weight'.format(outvis.replace('.ms','.im')))
                os.system('rm -rf {0}.residual'.format(outvis.replace('.ms','.im')))

            
                hdu = fits.open('{0}.pbeam.fits'.format(outvis.replace('.ms','.im')))[0]
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
                mask  = np.logical_or(mask0,mask1)
                mask  = np.logical_or(mask,data==0)
                data[mask] = np.nan

                hdu.data[0,0] = np.copy(data)
                hdu.writeto('{0}.pbeam.fits'.format(outvis.replace('.ms','.im')),overwrite=True)

        os.system('rm -vf *.last')
        os.system('rm -vf *.log')
    
    def uvloader(self, spw, field):

        ms=casatools.ms()
        ms.open(self.ms_file)
        ms.selectinit(reset=True)
        ms.selectinit(datadescid=int(spw))
        ms.select({'field_id': int(field)})

        u = np.copy(ms.getdata(['u'])['u'])
        v = np.copy(ms.getdata(['v'])['v'])
        freqs  = ms.range('chan_freq')['chan_freq'][:,0]

        ms.close()

        uwave  = (u*freqs[0]/const.c.value).flatten()
        vwave  = (v*freqs[0]/const.c.value).flatten()
        uvfreq = (np.zeros(np.shape(u))+freqs[0]).flatten()
                
                
        uvdata = np.array([uwave, vwave, uvfreq])
        return uvdata        
    
    def uvdata_loader(self):

        UVreal = np.empty(0)
        UVimag = np.empty(0)
        uvdist = np.empty(0)
        uvwghts = np.empty(0)
        
        for spw in self.spws:
            for field in self.fields:
                
                ms=casatools.ms()
                ms.open(self.ms_file)
                ms.selectinit(reset=True)
                ms.selectinit(datadescid=int(spw))
                ms.select({'field_id': int(field)})
                
                rec = ms.getdata(['u','v','data','weight'])
                uvreal = ((rec['data'][0][:].real+rec['data'][1][:].real)/2.0)
                uvimag = ((rec['data'][0][:].imag+rec['data'][1][:].imag)/2.0)
                uvwght = 4.0/(1.0/rec['weight'][0]+1.0/rec['weight'][1])
                
                u = rec['u']
                v = rec['v']
                freqs  = ms.range('chan_freq')['chan_freq'][:,0]                
                
                ms.close()
                
                uwave  = (u*freqs[0]/const.c.value).flatten()
                vwave  = (v*freqs[0]/const.c.value).flatten()
                
                uvdist  = np.append(uvdist, (uwave**2 + vwave**2)**0.5*1e-3)
                UVreal  = np.append(UVreal, np.mean(uvreal, axis = 0)*1e3)
                UVimag  = np.append(UVimag, np.mean(uvimag, axis = 0)*1e3)
                uvwghts = np.append(uvwghts, uvwght*1e-6)
        
        return np.array([uvdist, UVreal, UVimag, uvwghts])
