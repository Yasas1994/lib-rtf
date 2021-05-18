
required_packages = ['biopython', 'numpy', 'dendropy', 'ete3', 'pandas']

#for package in required_packages:
#	check_and_install(package)

import numpy as np
from Bio import SeqIO 
import re 
#from skbio import io
#from skbio import DistanceMatrix
#from skbio.tree import nj
from time import time
#import skbio
import subprocess
import sys
import argparse
#import gensim 
#from gensim.models import Word2Vec
#from gensim.test.utils import common_texts, get_tmpfile
from itertools import product
import os
import psutil
import pandas as pd
import dendropy
from ete3 import Tree

process = psutil.Process(os.getpid())



class ArgumentsParser(argparse.ArgumentParser):

    def __init__(self,*args, **kwargs):
        super(ArgumentsParser, self).__init__()


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

class Distance:

    '''Add more distance functions -> 
       Add capability to take average distance of 2 or more distance measures
       Add capability to take average distance for a range of k values
    '''

    
    def __init__(self,a,b):
        
        self.a = np.array(a,'float16') + 1.0e-17
        self.b = np.array(b,'float16') + 1.0e-17
        #print(self.a.shape)
        
        if len(self.a.shape)==1:
            self.a = self.a.reshape(1,-1)
            self.b = self.b.reshape(1,-1)
        else:
            pass
        
    def cosine_similarity_(self):

        '''calculates cosine similarity'''
        from sklearn.metrics.pairwise import cosine_similarity
        return cosine_similarity(self.a,self.b)
    
    
    def angular_distance(self):

        cs = self.cosine_similarity_()
        '''angular_distance = cos-1(cosine_similarity)/pi'''
        return (2*np.arccos(np.round(cs, 6)))/np.pi
    
    def euclidean_distance(self):

        from sklearn.metrics.pairwise import euclidean_distances
        return euclidean_distances(self.a,self.b)
        
class Sequences:
     
    
    def __init__(self, filename, format_='fasta'):
        
        self.filename  = filename
        self.format_ = format_
        

    def __iter__(self):
        
        fasta_sequences = SeqIO.parse(open(self.filename),self.format_)
        for seq in fasta_sequences:

            sequence = str(seq.__dict__['_seq'])
            id_ = seq.__dict__['id']
            matchobj1 = re.match(r'(\S+)', id_)
            #matchobj = re.match(r'(\S+\|\S+)',id)
            matchobj2 = re.search(r'[^ATGC]', sequence)

            if matchobj2:

                print('WARNING ! ' + id_ + ' has invalid character/s ', matchobj2.group())
                class_ = matchobj1.group(1)
                yield sequence.upper(), class_

            if matchobj1 and not matchobj2 :
    
                class_ = matchobj1.group(1)
                yield sequence.upper(), class_

class ReturnTimeCounts:
    '''
    Allow the users to prebuild the kmer index and to use user defined k-mer indices \
    output the k-mer index 
    
    '''
    '''
    simplified algorithm

    this algorithm iterates twice through the sequence dataset.
    In the first iteration, it creates a hash table for all possible return times for all sequences present in the sequence dataset
    in the second iteration, a hash table per each sequence is populated by counting return-times by sliding a window along the entirety of the sequnce

    '''
    @staticmethod
    def gapped_kmer(kmer, positions):
        t = str()
        for i in range(0,len(positions)):
            if positions[i] == '1':
                t = t + kmer[i]
        return t
    
    
    @staticmethod
    @timeit
    def cal_rt_vectors(sequences, k, positions = None):
        '''this function returns return-time counts for nucleotide sequences
        Parameters
        ----------

        sequences:(str[list]) 
        contains nucleotide sequences 

        k:(int)
        sets the kmer size for return times count calculation

        positions:str
        defines on and off positions of gapped kmers
        '''

        dic_global = {}
        out_put_vectors = [] #stores output return times vectors
        label_output = [] #stores sequence labels

        '''first iteration
        builds a hash table of all kmers present in the dataset
        '''

        for seq,label in sequences: # iterate through every sequence in the dataset
            dic_indices = {}
            for i in range(0,len(seq)-(k + 1)): # from sequence index 0 to seq.length - (kmer_size + 1)

                kmer = seq[i:i+k] # extracts kmer substring from the seq
                if positions != None:
                    kmer = ReturnTimeCounts.gapped_kmer(kmer = kmer, positions = positions) # returns the gapped kmer
                #print(kmer)

                if kmer in dic_indices: # if the dic has already seen the kmer 

                    dic_global[kmer+str(i-(dic_indices[kmer])-1)] = 1 #creates hash keys ;i-(dic[kmer])-1 is the return time for example AC10 - AC occuring after an offset of 10 nucleotides
                dic_indices[kmer] = i # this dictionary keeps track of the index at which a particular kmer is last seen 

        #sort dic_global and create new hash tables to be used in counting return times
        sorted_global = sorted(dic_global)

        '''second iteration
        populates the hash table built in the first iteration with return times counts
        '''
        for seq,label in sequences:
            mem = {} # this dictionary keeps track of the index at which a particular kmer is last seen 
            rt_counts = dict.fromkeys(sorted_global,0) # this hash table is populated with return times in the next iteration
            label_output.extend([label]) # appends sequence labels to an array

            for i in range(0,len(seq)-(k + 1)):

                kmer = seq[i:i+k]
                if positions != None:
                    kmer = ReturnTimeCounts.gapped_kmer(kmer = kmer, positions = positions)

                #print(kmer)
                if kmer in mem:
                    rt_counts[kmer+"{}".format(i-(mem[kmer])-1)] += 1 #i-(dic[kmer])-1 is the return time
                mem[kmer] = i

            out_put_vectors.append(np.array([*rt_counts.values()])) # appends the counts to output vector
            del mem
            del rt_counts
        return np.array(out_put_vectors),label_output,sorted_global



class ReturnTimeDistribution:

    @staticmethod
    def cal_rtd_vectors(sequences, k, positions = None):
        '''this function returns return-time counts for nucleotide sequences
        Parameters
        ----------
        
        sequences:(str[list]) 
        contains nucleotide sequences 
        
        k:(int)
        sets the kmer size for return times count calculation
        
        positions:str
        defines on and off positions of gapped kmers
        '''
        
        dic_global = {}
        out_put_vectors_mean = []
        out_put_vectors_sd = []
        label_output = []
        
        '''first iteration
        builds a hash table of all kmers present in the dataset
        '''
            
        for seq,label in sequences: # iterate through every sequence in the dataset
            dic_indices = {}
            for i in range(0,len(seq)-(k + 1)): # from sequence index 0 to seq.length - (kmer_size + 1)

                kmer = seq[i:i+k] # extracts kmer substring from the seq
                if positions != None:
                    kmer = ReturnTimeCounts.gapped_kmer(kmer = kmer, positions = positions) # returns the gapped kmer
                #print(kmer)

                if kmer in dic_indices: # if the dic has already seen the kmer 

                    dic_global[kmer] = 1 #creates hash keys ;i-(dic[kmer])-1 is the return time for example AC10 - AC occuring after an offset of 10 nucleotides
                dic_indices[kmer] = i # this dictionary keeps track of the index at which a particular kmer is last seen 

        #sort dic_global and create new hash tables to be used in counting return times
        sorted_global = sorted(dic_global)
        
        '''second iteration
        populates the hash table built in the first iteration with return times counts
        '''
        for seq,label in sequences:
            mem = {} # this dictionary keeps track of the index at which a particular kmer is last seen 
            rt_counts = dict(zip(sorted_global,[[] for i in range(len(sorted_global))])) # this hash table is populated with return times in the next iteration
            label_output.extend([label]) # appends sequence labels to an array
            
            for i in range(0,len(seq)-(k + 1)):

                kmer = seq[i:i+k]
                if positions != None:
                    kmer = ReturnTimeCounts.gapped_kmer(kmer = kmer, positions = positions)

                #print(kmer)
                if kmer in mem:
                    rt_counts[kmer].append(i-(mem[kmer])-1) #i-(dic[kmer])-1 is the return time
                mem[kmer] = i
                
            temp_values_std = [] # stores standard deviations for a single sequence
            temp_values_mean = [] # stores means for a single sequence
            
            for i in rt_counts: # iterates through keys of return time counts dictionary
                if rt_counts[i]: # checks whether the list is not empty

                    temp_values_std.append(np.std(np.asarray(rt_counts[i],'int32'))) # updates the temporary array with calculated SDs
                    temp_values_mean.append(np.average(np.asarray(rt_counts[i],'int32'))) #updates the temporary array with calculated means           
                                        
                elif not rt_counts[i]: # checks whether the list is empty

                    temp_values_mean.append(0) # if the list is empty append zero         
                    temp_values_std.append(0) 
                    
                    
            temp_values_std = np.asarray(temp_values_std) # converts temp array to a numpy array
            temp_values_mean = np.asarray(temp_values_mean) # converts temp array to a numpy array
            
            out_put_vectors_mean.append(temp_values_mean) # appends the counts to output vector
            out_put_vectors_sd.append(temp_values_std)

        # 0 - out_put_vectors_mean, 1 -  out_put_vectors_sd, 2 - label_output, 3 - sorted_global
        return np.asarray(out_put_vectors_mean),np.asarray(out_put_vectors_sd),label_output,sorted_global

class GenerateTree:

    @staticmethod
    @timeit
    def ReturnTimeCountsTree(RTCVector, labels, filename):

        DistanceFile = 'DMatrix_RTF{}'.format(filename.split(".")[0])

        print("Length of RTC vector: ", len(RTCVector[0]))
        dist = Distance(RTCVector,RTCVector)
        matrix = np.nan_to_num(dist.angular_distance())
        del dist
        
        WriteToPhylipDM(labels, matrix, filename = DistanceFile)
        a = subprocess.run(['rapidnj_x64.exe', '-d','-x',filename,'-i','pd',DistanceFile], capture_output=True)


        if a.returncode == 0:
            print ("done")
        else:
            print(a)

    
    @staticmethod
    def ReturnTimeDistributionTree(RTDVector_mean, RTDVector_sd, labels, filename):
        
        DistanceFile = 'DMatrix_RTD{}'.format(filename.split(".")[0])

        dist_mean = Distance(RTDVector_mean,RTDVector_mean)
        dist_sd = Distance(RTDVector_sd,RTDVector_sd)

        matrix_mean = np.nan_to_num(dist_mean.euclidean_distance())
        matrix_sd = np.nan_to_num(dist_sd.euclidean_distance())
        matrix = np.nan_to_num((matrix_mean + matrix_sd)**0.5)
        
        dm = pd.DataFrame(matrix, columns=labels, index=labels)
        dm.to_csv(DistanceFile, sep='\t')

        #newick_str = nj(dm, result_constructor=str)
        pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(src=open(DistanceFile), delimiter="\t")
        nj_tree = pdm.nj_tree()
        newick_str = str(nj_tree) + ";"
        #writes the tree to a file
        file = open('{}.nwk'.format(filename),'w')
        file.write(newick_str)
        file.close()

class dna2vec:


    @staticmethod
    def cal_dna2vec_vectors(sequences, model = "dna2vec_3.model", k = 3):
        
        model = Word2Vec.load(model)
        outputVector = []
        outputLabels = []
        
        for sequence, label in sequences:

            dna2vecVector = np.zeros(300)
                
            for i in range(len(sequence)):
                
                substring = sequence[i:i+3]
                
                if substring in model.wv:
                    
                    vector = model.wv[substring]
                    
                    dna2vecVector += vector
                        
            dna2vecVector = dna2vecVector
            outputVector.append(dna2vecVector)
            outputLabels.append(label)
            #clearprint(dna2vecVector)
        return np.array(outputVector),outputLabels

if __name__ == "__main__":


    parser = ArgumentsParser(description = 'Return Times based phylogeny')
    parser.add_argument('-f','--filename', type = str, required = True, metavar = ' ', help = 'path to file cotaining input sequences')
    parser.add_argument('-k', '--ksize', type = str, required = True, metavar = ' ', help = 'set k-mer size')
    parser.add_argument('-p', '--positions', type = str, metavar = ' ', help = 'binary pattern which defines the on and off positions of k-mers')
    parser.add_argument('-o', '--out', type = str, required = True, metavar = ' ', help = 'path to output')
    parser.add_argument('-t', '--type', type = str, required = True, metavar = ' ', help = 'specify the sequence type. i.e Nuceotides - N, Proteins - P')
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    args = parser.parse_args()
    #print(ReturnTimeCounts.gapped_kmer('ATGC','1111'))
    #X_mean,X_sd,y_,sorted_global = ReturnTimeDistribution.cal_rtd_vectors(Sequences('ecoli.fasta','fasta'),k=21,positions=None)
    #GenerateTree.ReturnTimeDistributionTree(X_mean, X_sd, y_, 'ecoli_rtd_k21')
    
    X, y_, sorted_global = ReturnTimeCounts.cal_rt_vectors(Sequences(args.filename,'fasta'), k = int(args.ksize), positions = args.positions)
    #X, y_ = dna2vec.cal_dna2vec_vectors(Sequences(args.filename,'fasta'))
    #GenerateTree.ReturnTimeCountsTree(X, y_, args.out)
    #m0 = process.memory_info().rss/(1024*1024) # in bytes 
    #X, y_, sorted_global = ReturnTimeCounts.cal_rt_vectors(Sequences('Hep.txt','fasta'), k = int(1))
    #m1 = process.memory_info().rss/(1024*1024)  # in bytes 
    GenerateTree.ReturnTimeCountsTree(X, y_, args.out)
    #m2 = process.memory_info().rss/(1024*1024)
    #print(process, m2 - m0)

