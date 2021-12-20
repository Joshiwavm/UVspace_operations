from src.Wrapper import *
from src.MsManager import *
from src.UnitTransform import *
from src.FourierTransform import *

from astropy.io import fits
from astropy.wcs import WCS

"""
To do:            
    - Do the fourier transform of the image/cube
    
    - Get spws and fields in MsManager direcelty from the listobs (don't know if that is possible)
    - put residues and model back into new ms files
"""

class Modeler:
    def __init__(self, 
                 filename_samples, 
                 filename_ms, 
                 obs_type, 
                 outputdir = './output/',
                 save = False):
        
        self.file = filename_samples
        self.msfile = filename_ms
        self.outputdir = outputdir
        self.save = save
        
        reader = FileReader()
        self.popt = reader.execute(self.file)
        self.info = getinfo()
        
        self.obs_type = obs_type

        self.msreader = MsReader(self.msfile, self.obs_type)        
        self.msreader.prep_pb()
    
    def _load_pb(self, spw, field):
        outvis = self.msreader.binvis.replace('-fid','-{0}'.format(field))
        outvis = outvis.replace('-sid','-{0}'.format(spw))
                
        pb, he_pb = fits.getdata('{0}.pbeam.fits'.format(outvis.replace('.ms','.im')), header = True)
        pb_im = pb[0,0] 
        pb_im[np.isnan(pb_im)] = 0.0
        
        self.pb_im = pb_im
        self.pb_hb = he_pb

    ################################################################################
        
    def _make_grid(self):
        wcs_pb = WCS(self.pb_hb, naxis = 2)
        
        x, y  = np.meshgrid(np.arange(0, len(self.pb_im), 1.), np.arange(0, len(self.pb_im), 1.))
        x += 0.5
        y -= 0.5
        
        sky = wcs_pb.pixel_to_world(x, y)
        RA  = sky.ra*u.deg
        Dec = sky.dec*u.deg

        r  =  (Dec-np.mean(Dec))**2 
        r += ((RA - np.mean(RA))*np.cos(np.deg2rad(np.mean(Dec.value))))**2
        r  = r**0.5 
        return r
    
    def _make_modelimage(self):
        grid  = self._make_grid()
        image = np.zeros_like(grid.value)

        for t  in self.popt.keys(): 
            if t == 'Link':
                continue
                
            if COMPONENTS[t]['make_image'] == True:                
                transform = TransformInput(self.popt[t], t)
                input_par = transform.generate() 
                
                coord = grid * self.info.cosmo.angular_diameter_distance(self.popt[t]['z'])
                coord = coord.value/input_par['major']
                
                info = getinfo()
                
                if t == 'gnfwPressure':
                    rs = info.linmesh[1:]
                else:
                    rs = info.linmesh
                    
                integrated_P = COMPONENTS[t]['function'](rs, **input_par)
                image += np.interp(coord, rs, integrated_P) * self.info.ysznorm.value

        if self.save: np.save(self.outputdir + 'noiseless_model', image)
        return image
    
    ################################################################################

    def _partouv(self):
        spwpoints = np.zeros(np.shape(self.uvdata[0]),dtype=np.complex128)
        
        for i in np.unique(self.popt['Link']):
            model_type, spectrum_type = np.array(list(self.popt)[:len(self.popt['Link'])])[(np.array(self.popt['Link']) == i)]
            
            if COMPONENTS[model_type]['make_image'] == True:
                continue
            
            spwpoints = self.FT.partouv(spwpoints, self.popt[model_type], self.popt[spectrum_type], model_type, spectrum_type)
        return spwpoints
    
    
    def _get_analytical_model(self):
        vis_model = np.empty(0)
        datas    = np.empty((3,0))
        
        for spw in self.msreader.spws:
            for field in self.msreader.fields:
                
                self._load_pb(spw, field)
                self.uvdata = self.msreader.uvloader(spw, field) 
                self.FT     = FT(self.uvdata, self.pb_hb, self.pb_im)
                
                vis_model = np.append(vis_model, self._partouv())
                datas     = np.hstack((datas, self.uvdata))
        return vis_model, datas
    
    def run(self):        
        vis_model, datas = self._get_analytical_model()
            
        check = 0
        for key in self.popt.keys():
            if key == 'Link':
                continue
            check += COMPONENTS[key]['make_image'] 

        if check > 0: 
            image = self._make_modelimage()
            
        return vis_model, self.popt['Scaling'], datas
