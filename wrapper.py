from settings import COMPONENTS
import numpy as np
import scipy

class FileReader:
    def __init__(self, types):
        self.labels = []
        for t in types:
            self.labels.append(COMPONENTS[t]["variables"])
        
    def _mapper(self, data):
        # Need to account for fixed parameters in the sampling.
        mapped_data = {}
        for label in self.labels:
            for key, value in zip(label, data):
                mapped_data[key] = value
        return mapped_data

    def _read(self, filename):
        # Need to account for sampling method. For now based on nested sampling
        results = np.load(filename, allow_pickle=True)

        samples = np.copy(results['samples'])
        weights = results['logwt']-scipy.special.logsumexp(results['logwt']-results['logz'][-1])
        weights = np.exp(weights-results['logz'][-1])
        popt    = np.average(samples, weights = weights, axis = 0)

        return popt

    def execute(self, filename):
        data = self._read(filename)
        mapped_data = self._mapper(data)
        return mapped_data

def main():
    filename = '/scigarfs/home/jvmarrewijk/eszee/outputdir/RXC_J2014.8_ALMA_pointonly_pickle'
    types  = ["pointSource", "powerLaw"]
    reader = FileReader(types)
    mapped_data = reader.execute(filename)
    print(mapped_data)

if __name__ == '__main__':
    main()