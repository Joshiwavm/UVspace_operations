import numpy as np
import astropy.constants as const
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from wrapper import FileReader

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

class TransformInput():
    
    def __init__(self, 
                 popt, 
                 typ, 
                 c500=1.00, 
                 mass=1.00,
                 limdist=np.inf,
                 epsrel=1.00E-06,
                 freeLS=None):
        
        self.popt = popt
        self.type = typ
        self.c500 = c500
        self.mass = mass
        self.limdist = limdist
        self.epsrel = epsrel
        self.freeLS = freeLS
        
    def _A10Pressure(self, info):
        param_input = {
            "offset":self.popt['Offset'],
            "amp":self.popt['log10'],
            "major":self.popt['c500'],
            "e":self.popt['e'],
            "alpha":self.popt['Alpha'],
            "beta":self.popt['Beta'],
            "gamma":self.popt['Gamma'],
            "ap":self.popt['Alpha_p'],
            "c500":self.c500,
            "mass":self.mass,
            "limdist":self.limdist,
            "epsrel":self.epsrel,
            "freeLS":self.freeLS
        }
        
        Hz = info.cosmo.H(self.popt['z'])

        c500 = self.popt['c500']
        m500 = self.popt['log10']
        bias = self.popt['bias']

        r500 = ((3.00/4.00/np.pi/500.00/info.cosmo.critical_density(self.popt['z']))*(1.00-bias)*(10**m500)*u.solMass)**(1.00/3.00)

        param_input['major'] = r500.to(u.Mpc).value/c500
        param_input['amp']   = self.popt['P0']*(3.00/8.00/np.pi)*(info.fb*info.mu/info.mue)
        param_input['amp']  *= (((((2.5e2*Hz*Hz)**2.00)*((1.00-bias)*(10**(m500-15.00))*u.solMass)/(const.G**(0.50)))**(2.00/3.00)).to(u.keV/u.cm**3)).value
        param_input['amp' ] *= 1e10

        if self.popt['Alpha_p']==-0.10:
            self.type = 'gnfwPressure'
            return self.generate() #NOT SURE IF THIS IS WORKING!!!!
        else: 
            mass = ((1.00-bias)*(10**(m500-14.00))*(info.cosmo.H0.value/70.00)/3.00)

        limdist = info.limdist*c500 if np.isfinite(info.limdist) else np.inf
        freeLS = self.popt['depth'] if self.type =='A10PressureLS' else None

        param_input['c500'] = c500
        param_input['mass'] = mass
        param_input['limdist'] = limdist
        param_input['freeLS'] = freeLS
        return param_input
    
    def _gnfwPressure(self, info):
        param_input = {
            "offset":self.popt['Offset'],
            "amp":self.popt['Amplitude'],
            "major":np.deg2rad(self.popt['Major']),
            "e":self.popt['e'],
            "alpha":self.popt['Alpha'],
            "beta":self.popt['Beta'],
            "gamma":self.popt['Gamma'],
            "limdist":self.limdist,
            "epsrel":self.epsrel,
            "freeLS":self.freeLS
        }
               
        param_input['major'] *= info.cosmo.angular_diameter_distance(self.popt['P0'])        
        return param_input
    
    def _betaPressure(self, info): #done
        param_input = {
            "offset": self.popt['Offset'] ,
            "amp": self.popt['Amplitude'] ,
            "major": np.deg2rad(self.popt['Major']),
            "e": self.popt['e'] ,
            "beta": self.popt['Beta'],
            "limdist": self.limdist,
            "epsrel": self.epsrel,
            "freeLS": self.freeLS
        }
        
        param_input['major'] *= info.cosmo.angular_diameter_distance(self.popt['z']).value
        return param_input
    
    def generate(self):        
        if self.type=='betaPressure': 
            return self._betaPressure(getinfo())
        
        elif self.type=='gnfwPressure': 
            return self._gnfwPressure(getinfo())
         
        elif (self.type=='A10Pressure') or (self.type == 'A10PressureLS'):
            return self._A10Pressure(getinfo())

def main():
    names    = ['RXC_J2014.8_ALMA_pointonly', 'RXC_J2014.8_ALMA_bothpoint_and_SZ', 'RXC_J2014.8_ACA_and_ALMA_null_test_pointsource']
    name     = names[1]
    filename = '/scigarfs/home/jvmarrewijk/eszee/outputdir/'+name+'_pickle'
    
    reader      = FileReader()
    mapped_data = reader.execute(filename)
    
    typ = 'A10Pressure'
    print(len(mapped_data))
    print(mapped_data[typ])
    
    transform  = TransformInput(mapped_data[typ], typ)
    input_data = transform.generate()
    print()
    print(input_data)
    
if __name__ == '__main__':
    main()