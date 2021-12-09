from settings import COMPONENTS
from wrapper import *
from models import *
from FourierTransform import *
from msmanager import *
from UnitTransform import *

from astropy.io import fits

class Modeler:
    def __init__(self, filename_samples, filename_ms, filename_pb, outputdir = './plots/',  save = False):
        self.file = filename_samples
        self.msfile = filename_ms
        self.pbfile = filename_pb
        self.outputdir = outputdir
        self.save = save

        reader = FileReader()
        self.popt = reader.execute(self.file)

        self.pbeam = self._load_pb()

    def _load_ms(self):
        #use msmanager.py file here
        pass
    
    def _load_pb(self):
        """ For now I use the pb file out of output directory made with eszee. 
        Better would be to make a dirty image from the ms file and grab the created .pb image. 
        """
        pb, he_pb = fits.getdata(self.pbfile, header = True)
        pb_im = pb[0,0]
        pb_im[np.isnan(True)] = 0.0
        
        return pb_im

    def _make_grid(self):
        xx, yy = np
        pass
    
    def _make_modelimage(self):

        #use pb for a grid

        grid  = 
        image = np.zeros_like(grid)

        for t  in self.popt.keys():
            if COMPONENTS[t]['make_image'] == True:
                image += COMPONENTS[t]['function'](grid, **self.popt[t])

        if self.save: np.save(self.outputdir + 'noiseless_model', image)
        return image

    def run(self):        
        
        #check if we need to fourier transform
        check = 0
        for key in self.popt.keys():
            check += COMPONENTS[key]['make_image'] 

        if check > 0: 
            image = self._make_modelimage()
            fourier = True
        else: 
            fourier = False
        
        if fourier:
            pass
            # Fourier = FT()
            # out_uv = Fourier.execute(image)
        else:
            pass
            # directly to the visibility modeling

def main(filename_samples, filename_ms, filename_pb): 

    obj = Modeler(filename_samples, filename_ms, filename_pb, save = True)
    out = obj.run()

    print(80*"#")
    print(out)
    print(80 * "#")

if __name__ == '__main__':
    filename_samples  = '/scigarfs/home/jvmarrewijk/eszee/outputdir/RXC_J2014.8_ALMA_pointonly_pickle'
    filename_ms       = 'some/name/to/ms/file'
    filename_pb       = 'some/path'

    main(filename_samples, filename_ms, filename_pb)
