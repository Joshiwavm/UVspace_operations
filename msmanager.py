#import casacore packages and ms file managements for the fourier transform

class MsReader:
    def __init__(self) -> None:
        pass

    def _read(self, filename):
        pass
    
    def execute(self, filename):
        pass
    
def main():
    filename = '/path/to/filename'

    reader = MsReader()
    mapped_data = reader.execute(filename)
    print(mapped_data)

if __name__ == 'main':
    main()