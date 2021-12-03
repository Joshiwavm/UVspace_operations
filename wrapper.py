from settings import COMPONENTS
import numpy as np
import scipy

class FileReader:

    def _mapper(self, params):

        self.types.append('Scaling')        
        mapped_types = {}
        
        count = 0
        for idx, t in enumerate(self.types):
            mapped_params = {}
            for key, value in zip(self.labels[idx], params[count: count + len(self.labels[idx]) + 1]):
                mapped_params[key] = value
 
            count += len(self.labels[idx])
            mapped_types[t] = mapped_params
        return mapped_types
               
        
    def _read_combinedata(self, results, popt, fixd):
    
        params = []
        count = 0

        for idx, lbl in enumerate(self.labels):
            for jdx, l in enumerate(lbl):
                if fixd[idx][jdx]: 
                    params.append(popt[count])
                    count +=1
                else:
                    if COMPONENTS[self.types[idx]]['spectrum'] == False:
                        comp = 'model'
                    else:
                        comp = 'spectrum'
                        
                    for i in range(len(results['pars'])):
                        if results['pars'][i][comp]['type'] == self.types[idx]:
                            index = i
                       
                    params.append(results['pars'][index][comp]['guess'][jdx])                
                
        return params
        
    def _read_fixedvalues(self, results, popt):
        
        fixed = []
        for var in results['vary'][:-1]:
            for v in var['values']:
                fixed.append(list(var['values'][v]['vary']))
        
        scaling = []
        for scle in results['vary'][-1]['values']['vary']:
            scaling.append(scle)
        
        fixed.append(list(scaling))
        return fixed
    
    def _read_popt(self, results):
        samples = np.copy(results['samples']['samples'])
        weights = results['samples']['logwt']-scipy.special.logsumexp(results['samples']['logwt']-results['samples']['logz'][-1])
        weights = np.exp(weights-results['samples']['logz'][-1])
        popt    = np.average(samples, weights = weights, axis = 0)
        return popt
    
    
    def _read_labels(self,results):
        self.labels = []
        for t in self.types:
            self.labels.append(COMPONENTS[t]["variables"])
    
        scaling = []
        for i in range(len(results['vary'][-1]['values']['vary'])):
            scaling.append('scale'+str(i))
            
        self.labels.append(scaling)
        
    def _read_types(self, results):
        self.types = []
        for res in results['pars']:
            for r in res:
                self.types.append(res[r]['type'])

    def _read(self, filename):
        pickled_file = np.load(filename, allow_pickle=True)   
        
        self._read_types(pickled_file) 
        self._read_labels(pickled_file) 

        if pickled_file['type'] == 'dynesty':
            popt = self._read_popt(pickled_file)
            fixd = self._read_fixedvalues(pickled_file, popt)
            params = self._read_combinedata(pickled_file, popt, fixd)
            return params
        else:
            printError('Provide correct sampler type')

    def execute(self, filename):
        params = self._read(filename)
        mapped_data = self._mapper(params)
        return mapped_data

def main():
    names    = ['RXC_J2014.8_ALMA_pointonly', 'RXC_J2014.8_ALMA_bothpoint_and_SZ', 'RXC_J2014.8_ACA_and_ALMA_null_test_pointsource']
    name     = names[1]
    filename = '/scigarfs/home/jvmarrewijk/eszee/outputdir/'+name+'_pickle'
    
    reader      = FileReader()
    mapped_data = reader.execute(filename)
    print(mapped_data)

if __name__ == '__main__':
    main()
