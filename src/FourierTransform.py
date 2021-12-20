from src.Utillities import *

import scipy.integrate as scintegr
import scipy.interpolate as scinterp
from galario.double import apply_phase_vis, sampleImage

class FT():
    def __init__(self, uvdata, pb_head, pb_im):
        self.uvdata   = uvdata
        self.pb_head  = pb_head
        self.pb_im    = pb_im
        
    ################################################################################
        
    def _spw_PointSource(self, model, reffreq):
        spwRA   = self.pb_head['CRPIX1'] # I think CRPIX1 is in pixel
        spwRA  += (model['RA']-self.pb_head['CRVAL1'])*np.cos(np.deg2rad(self.pb_head['CRVAL2']))/self.pb_head['CDELT1'] 
                        
        spwDec  = self.pb_head['CRPIX2']+(model['Dec']-self.pb_head['CRVAL2'])/self.pb_head['CDELT2']
        
        spwairy = scinterp.RectBivariateSpline(np.linspace(0,self.pb_head['NAXIS2'],self.pb_head['NAXIS2']),
                                               np.linspace(0,self.pb_head['NAXIS1'],self.pb_head['NAXIS1']),
                                               self.pb_im,kx=1,ky=1,s=0)
        
        spwbeam = spwairy.ev(spwDec,spwRA)
        return spwbeam
    
    def _spw_GaussSource(self, model, reffreq, spwpoints):
        pass
    
    def _spw_GaussSurface(self, model, reffreq, spwpoints):
        pass
    
    def _spw_powerLaw(self, model, spectrum, reffreq, spwbeam):
        spwpoint  = spwbeam*model['Amplitude']
        
        index     = spectrum['SpecIndex']
        spwpoint *= ((self.uvdata[2]/reffreq)**(index))+0j
        return spwpoint
    
    def _spw_powerLawMod(self, model, spectrum, reffreq, spwbeam):
        spwpoint  = spwbeam*model['Amplitude']
        
        index     = spectrum['SpecIndex']+spectrum['SpecCurv']*np.log(self.uvdata[2]/reffreq)
        spwpoint *= ((self.uvdata[2]/reffreq)**(index))+0j
        return spwpoint
    
    def _spw_powerDust(self, model, spectrum, reffreq, spwbeam):
        pass
    
    ################################################################################
    
        
    def partouv(self, spwpoints, model, spectrum, model_type, spectrum_type):
        
        reffreq   = self.pb_head['RESTFRQ']

        if model_type == 'pointSource':
            spwbeam = self._spw_PointSource(model, reffreq) 
        elif model_type == 'gaussSource':
            spwbeam = self._spw_GaussSource(model, reffreq)
        elif model_type == 'gaussSurface':
            spwbeam = self._spw_GaussSurface(model, reffreq)

        if spectrum_type == 'powerLaw':
            spwpoint = self._spw_powerLaw(model, spectrum, reffreq, spwbeam)
        elif spectrum_type == 'powerLawMod':
            spwpoint = self._spw_powerLawMod(model, spectrum, reffreq, spwbeam)    
        elif spectrum_type == 'powerDust':
            spwpoint = self._spw_powerLawMod(model, spectrum, reffreq, spwbeam)    
        elif spectrum_type == 'tSZ':
            pass

        spwpoints += apply_phase_vis((np.deg2rad(model['RA']-self.pb_head['CRVAL1']))*np.cos(np.deg2rad(reffreq)),
                                      np.deg2rad(model['Dec']-self.pb_head['CRVAL2']),
                                      self.uvdata[0], self.uvdata[1], spwpoint)
        
        return spwpoints
    
    
    def imtouv(self, im):
        pass
