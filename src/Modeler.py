from src.Wrapper import *
from src.MsManager import *
from src.UnitTransform import *
from src.FourierTransform import *

from astropy.io import fits
from astropy.wcs import WCS

"""
To do:
    - I have one pb for all three spws
    - I have to provide a pb for the data set, so I need to make it compatible for ALMA and ACA
        - I can make the pb from tclean the ms file to make a dirty image
        - Also the pb changes in size over frequency which makes the density profile more inacurate
            
    - Do the fourier transform of the image/cube
    - Do the alpha scaling in the end. Make sure it automatically provides the proper scaling to ALMA and ACA
    - Put it on Github
    
    -put residues and model  back into new ms files
"""


class Modeler:
    def __init__(self, 
                 filename_samples, 
                 filename_ms, 
                 filename_pb, 
                 outputdir = '/scigarfs/home/jvmarrewijk/eszee/plots/', #hardcoded quickfix
                 save = False):
        
        self.file = filename_samples
        self.msfile = filename_ms
        self.pbfile = filename_pb
        self.outputdir = outputdir
        self.save = save

        reader = FileReader()
        self.popt = reader.execute(self.file)

        self._load_pb()
        self.info = getinfo()
        
        msreader = MsReader(self.msfile)
        self.uvdata = msreader.uvloader() 
                
        self.FT = FT(self.uvdata, self.pb_hb, self.pb_im)
    
    def _load_pb(self):
        pb, he_pb = fits.getdata(self.pbfile, header = True)
        pb_im = pb[0,0] #hard coded quickfix
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
        r  = r**0.5 #degrees
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
    
    
    def _partouv(self):
        spwpoints = np.zeros(np.shape(self.uvdata[0]),dtype=np.complex128)
        
        for i in np.unique(self.popt['Link']):
            model_type, spectrum_type = np.array(list(self.popt)[:len(self.popt['Link'])])[(np.array(self.popt['Link']) == i)]
            
            if COMPONENTS[model_type]['make_image'] == True:
                continue
            
            spwpoints = self.FT.partouv(spwpoints, self.popt[model_type], self.popt[spectrum_type], model_type, spectrum_type)
        return spwpoints
            
    def run(self):        
        #Fourier Transform analytical Functions first
        vis_model = self._partouv()

        #check if an image needs to be transformed        
        check = 0
        for key in self.popt.keys():
            if key == 'Link':
                continue
            check += COMPONENTS[key]['make_image'] 

        if check > 0: 
            image = self._make_modelimage()
            
        #apply alpha scaling
            
        return vis_model, self.uvdata