from src.Utillities import *

import casatools
from casatools import ms

class MsReader:
    def __init__(self, filename):
        self.ms_file = filename

    def uvloader(self, spws = np.arange(3), fields = [0]): #hardcoded quickfix
        
        uwave  = np.empty(0)
        vwave  = np.empty(0)
        uvfreq = np.empty(0)

        for spw in spws:
            for field in fields:
                
                ms=casatools.ms()
                ms.open(self.ms_file)
                ms.selectinit(reset=True)
                ms.selectinit(datadescid=spw)
                ms.select({'field_id': field})

                u = np.copy(ms.getdata(['u'])['u'])
                v = np.copy(ms.getdata(['v'])['v'])
                freqs  = ms.range('chan_freq')['chan_freq'][:,0]

                ms.close()

                uwave  = np.append(uwave,   (u*freqs[0]/const.c.value).flatten())
                vwave  = np.append(vwave,   (v*freqs[0]/const.c.value).flatten())
                uvfreq = np.append(uvfreq, (np.zeros(np.shape(u))+freqs[0]).flatten())
                
                
        uvdata = np.array([uwave, vwave, uvfreq])
        return uvdata
    
    def uvdata_loader(self, spws = np.arange(3), fields = [0]):

        UVreal = np.empty(0)
        UVimag = np.empty(0)
        uvdist = np.empty(0)
        uvwghts = np.empty(0)
        
        for spw in spws:
            for field in fields:
                
                ms=casatools.ms()
                ms.open(self.ms_file)
                ms.selectinit(reset=True)
                ms.selectinit(datadescid=spw)
                ms.select({'field_id': field})
                
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

def main():
    filename = '/scigarfs/home/jvmarrewijk/eszee/RXC_J2014.8_Science_Data/RXC_J2014.8_ALMA_timebin30s_SZ_concat.ms'

    reader = MsReader(filename)
    uvdata = reader.uvloader()
    
    print(uvdata.shape)

if __name__ == '__main__':
    main()
    