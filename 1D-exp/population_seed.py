from individual_seed import INDIVIDUAL
import copy
import random
import pickle

class POPULATION:
    def __init__(self, popSize):
        self.p = {}
        self.popSize = popSize

    def Initialize(self):
        for i in range(0, self.popSize):
            self.p[i] = INDIVIDUAL(i)
            
    def Print(self):
        for i in self.p:
            if (i in self.p):
                self.p[i].Print()
        print()
    
    def Evaluate(self):
        for i in self.p:
            self.p[i].Compute_Fitness()
            
    def Mutate(self):
        for i in self.p:
            self.p[i].Mutate()
            
    def ReplaceWith(self, other):
        for i in self.p:
            if ( self.p[i].fitness < other.p[i].fitness ):
                self.p[i] = other.p[i]
    
    def Fill_From(self, other):
        self.Copy_Best_From(other)
        #self.Print()
        self.Collect_Children_From(other)
        #self.Print()
    
    def Copy_Best_From(self, other):
        for i in other.p:
            if i == 0:
                highest_fitness = other.p[i].fitness
                highest_index = i
            elif (other.p[i].fitness >= highest_fitness):
                highest_fitness = other.p[i].fitness
                highest_index = i
        
        self.p[0] = copy.deepcopy(other.p[highest_index])
        return highest_index
    
    def Collect_Children_From(self, other):
        for i in range(1,self.popSize):
            winner = other.Winner_Of_Tournament_Selection()
            self.p[i] = copy.deepcopy(winner)
            self.p[i].Mutate()
            
    def Winner_Of_Tournament_Selection(other):
        p1 = random.randint(0, other.popSize-1)
        p2 = random.randint(0, other.popSize-1)
        while p1 == p2:
            p2 = random.randint(0, other.popSize-1)
        if other.p[p1].fitness > other.p[p2].fitness:
            return other.p[p1]
        else:
            return other.p[p2]
    
    def SaveBestIndiv(self):
        for i in self.p:
            if i == 0:
                highest_fitness = self.p[i].fitness
                highest_index = i
            elif (self.p[i].fitness >= highest_fitness):
                highest_fitness = self.p[i].fitness
                highest_index = i
        self.p[highest_index].SaveBest()
        self.p[highest_index].Save()
    
    def ComputeAverage(self):
        f_sum = 0
        for i in self.p:
            f_sum += self.p[i].fitness 
        avg = f_sum / self.popSize
        
        f = open('savedAvgsSeed.dat', 'ab')
        pickle.dump(avg , f)
        f.close()
        return avg