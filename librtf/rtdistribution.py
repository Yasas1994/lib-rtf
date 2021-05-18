from librtf.rtutils import timeit
from librtf.rtcounts import ReturnTimeCounts
import numpy as np

class ReturnTimeDistribution:

    @staticmethod
    @timeit
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
        return np.asarray(out_put_vectors_mean),np.asarray(out_put_vectors_sd)