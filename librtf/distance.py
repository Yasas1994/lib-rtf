import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics.pairwise import pairwise_distances


class Distance:

    '''Add more distance functions -> 
       Add capability to take average distance of 2 or more distance measures
       Add capability to take average distance for a range of k values
        [‘cityblock’, ‘cosine’, ‘euclidean’, ‘l1’, ‘l2’, ‘manhattan’]
    '''

    
    def __init__(self,a,b):
        
        self.a = np.array(a,'float16') + 1.0e-17 # adds a small amount to avoid zerodivision error
        self.b = np.array(b,'float16') + 1.0e-17
        #print(self.a.shape)
        
        if len(self.a.shape)==1:
            self.a = self.a.reshape(1,-1)
            self.b = self.b.reshape(1,-1)
        else:
            pass
        
    def cosine_similarity_(self):

        '''calculates cosine similarity'''
        #from sklearn.metrics.pairwise import cosine_similarity
        return cosine_similarity(self.a,self.b)
    
    
    def angular_distance(self):

        cs = self.cosine_similarity_()
        '''angular_distance = cos-1(cosine_similarity)/pi'''
        return (2*np.arccos(np.round(cs, 6)))/np.pi
    
    def euclidean_distance(self):

        #from sklearn.metrics.pairwise import euclidean_distances
        return euclidean_distances(self.a,self.b)


    def cityblock_distance(self):

        return pairwise_distances(self.a,  metric='cityblock', n_jobs=1)

    def l1_distance(self):

        return pairwise_distances(self.a,  metric='l1', n_jobs=1)

    def l2_distance(self):

        return pairwise_distances(self.a,  metric='l2', n_jobs=1)
    
    def manhattan_distances(self):

        return pairwise_distances(self.a,  metric='manhattan', n_jobs=1)
