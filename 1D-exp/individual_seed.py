from switch import switch
import constants as c
import random
import math
import numpy as np
import sys
import pickle

class INDIVIDUAL:
    def __init__(self, i):
        self.N_light = 10
        self.N = 20
        self.genome = random.sample(range(0, self.N), self.N_light) #np.random.randint(0, high=c.GRID_SIZE-1, size=(c.FIBERS*2, 4), dtype='int')
        self.fitness = 0
        self.ID = i
    
    def Compute_Fitness(self, show=False):
        # wait for the simulation to end and get the fitness
        self.fitness = switch.evaluate(self.genome)
        if show:
            switch.showPacking(self.genome)
            print("fitness is:")
            print(self.fitness)
        return self.fitness
    
    def Mutate(self):
        geneToMutate = random.randint(0, self.N_light-1)
        candidate = random.randint(0, self.N-1)
        # check for uniqueness
        while candidate in self.genome:
            candidate = random.randint(0, self.N-1)
        self.genome[geneToMutate] = candidate
        
    
    def Print(self):
        print('[', self.ID, self.fitness, ']', end=' ')
        
    
    def Save(self):
        f = open('savedFitnessSeed.dat', 'ab')
        pickle.dump(self.fitness , f)
        f.close()
    
    def SaveBest(self):
        f = open('savedBestsSeed.dat', 'ab')
        pickle.dump(self.genome , f)
        f.close()