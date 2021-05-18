from time import time
from itertools import product



class Datetime:

    '''This class extents python time module'''
    
    
    
    def __init__(self):
        import time as time
        self.year  = time.gmtime().tm_year
        self.month = time.gmtime().tm_mon
        self.day = time.gmtime().tm_mday
        self.hour = time.gmtime().tm_hour 
        self.minu = time.gmtime().tm_min
        self.sec = time.gmtime().tm_sec
        

    def __str__(self):
        
        dateTime = f"{self.day}{self.year}{self.month}-{self.hour}{self.minu}{self.sec}"
        return dateTime
        

def check_and_install(package):
    import importlib
    try:
        importlib.util.find_spec(package)
    except ImportError:
        import pip
        pip.main(['install', package])


def timeit(function):
    '''This wrapper is used to evaluate execution times'''
    def wrapper (*args, **kwargs):
        start = time()
        result = function(*args, **kwargs)
        print ("time elapsed: {:.2f} s".format(time()-start))
        return result

    return wrapper

def GenerateKmers(kmer_size):
    '''This method generates all possible kmers of given size k'''
    p = product(['A','T','G','C'],repeat=kmer_size)
    a = map("".join,p)
    return list(a)

def WriteToPhylipDM(labels, matrix, filename):
    '''This method writes a phylip formatted distance matrix'''
    NumOfOtus = len(labels)
    with open(filename, 'w') as f:
        f.write("\t"+ str(NumOfOtus)+"\n")
        for i, j in zip(labels, matrix):
            FormatToString = '{}    {}\n'.format(i,"\t".join(j.astype('str')))
            f.write(FormatToString)

def WriteDicToFile(sorted_dict):
    
    '''This method writes dictionary keys to a file'''
    with open(f'kmers{Datetime()}','w') as file:
        file.write('\n'.join(sorted_dict))
