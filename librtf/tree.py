import pandas as pd
import dendropy
import numpy as np
from librtf.rtutils import timeit, WriteToPhylipDM
import subprocess
from librtf.distance import Distance

class GenerateTree:

    @staticmethod
    @timeit
    def ReturnTimeCountsTree(RTCVector, labels, filename, nj = None):

        DistanceFile = 'DMatrix_{}'.format(filename.split(".")[0])

        print("Length of RTC vector: ", len(RTCVector[0]))
        dist = Distance(RTCVector,RTCVector)
        matrix = np.nan_to_num(dist.angular_distance())
        del dist

        if nj == False: 
            WriteToPhylipDM(labels, matrix, filename = DistanceFile)
            a = subprocess.run(['rapidnj_x64.exe', '-d','-x',filename,'-i','pd',DistanceFile], capture_output=True)


            if a.returncode == 0:
                print ("................................")
            else:
                print(a)
        else:
            dm = pd.DataFrame(matrix, columns=labels, index=labels)
            dm.to_csv(DistanceFile, sep='\t')

            pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(src=open(DistanceFile), delimiter="\t")
            nj_tree = pdm.nj_tree()
            newick_str = str(nj_tree) + ";"
            #writes the tree to a file
            file = open('{}.nwk'.format(filename),'w')
            file.write(newick_str)
            file.close()
            print ("...................................")


    
    @staticmethod
    @timeit
    def ReturnTimeDistributionTree(RTDVector_mean, RTDVector_sd, labels, filename, nj=None):
        
        DistanceFile = 'DMatrix_{}'.format(filename.split(".")[0])
        print("Length of RTD vector: ", len(RTDVector_mean[0]))
        dist_mean = Distance(RTDVector_mean,RTDVector_mean)
        dist_sd = Distance(RTDVector_sd,RTDVector_sd)

        matrix_mean = np.nan_to_num(dist_mean.euclidean_distance())
        matrix_sd = np.nan_to_num(dist_sd.euclidean_distance())
        matrix = np.nan_to_num((matrix_mean + matrix_sd)**0.5)
        
        if nj == False:

            WriteToPhylipDM(labels, matrix, filename = DistanceFile)
            a = subprocess.run(['rapidnj_x64.exe', '-d','-x',filename,'-i','pd',DistanceFile], capture_output=True)

            if a.returncode == 0:
                print ("....................")
            else:
                print(a)

        else:
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
            print ("....................")
