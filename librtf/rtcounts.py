from librtf.rtutils import timeit
import numpy as np

class ReturnTimeCounts:
    '''
    Allow the users to prebuild the kmer index and to use user defined k-mer indices 
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
        '''
        Creates a gap reduced k-mer for a given binary pattern of NO and OFF positions

        kmer:(str)
        a nucleotide or protein k-mer e.x: ATTGCATG

        potitions:(str)
        a binary pattern string defining ON and OFF positions of the string e.x 111011101

        example:

        ATTGCATG + 11011011 ---> ATGCTG (note that nucleotides in the positions 3 and 6 have been removed)

        '''
        t = str() # stores the reduced gapped k-mer 
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
