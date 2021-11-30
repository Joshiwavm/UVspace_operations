from settings import COMPONENTS
from wrapper import *
from models import *
from fourier import *
from msmanager import MsReader

class Modeler:
    def __init__(self, types, filename, ms_filename, outputdir = './plots/'):
        self.types = types
        self.file = filename
        self.msfile = ms_filename
        self.outputdir = outputdir

        self._load_models()
        self._load_popt()
        self._load_ms()

    def _load_models(self):
        self.models = []
        for t in self.types:
            if 'function' in COMPONENTS[t]:
                self.models.append(COMPONENTS[t]['function'])
            else:
                self.models.append(None)

    def _load_popt(self):
        self.popt = []
        for t in self.types:
            reader = FileReader([t])
            self.popt.append(reader.execute(self.file))

    def _load_ms(self):
        #use msmanager.py file here
        pass

    def _make_modelimage(self, grid, save):
        image = np.zeros_like(grid)
        for idx, t  in enumerate(self.types):
            if 'function' in COMPONENTS[t]:
                image += self.models[idx](grid, **self.popt[idx])

        if save: np.save(self.outputdir + 'noiseless_model', image)
        return image

    def run(self, grid, save = False):        
        
        image = self._make_modelimage(grid, save)
        

        # Fourier = FourierTransform()
        # out_uv = Fourier.execute(image)

        # self._save(out_uv, obs)
        # return out_uv


    def _save(self):
        pass

def main(filename, ms_name, types): 

    grid = np.arange(0,1024,1)
    obj = Modeler(types, filename, ms_name)
    out = obj.run(grid)

    print(80*"#")
    print(out)
    print(80 * "#")

if __name__ == '__main__':
    filename  = '/scigarfs/home/jvmarrewijk/eszee/outputdir/RXC_J2014.8_ALMA_pointonly_pickle'
    ms_name   = 'some/name/to/ms/file'
    types     = ["pointSource", "powerLaw"] #should be inproper order as modelled

    main(filename, ms_name, types)