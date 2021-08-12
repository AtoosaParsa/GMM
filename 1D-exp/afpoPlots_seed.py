import pickle
import matplotlib.pyplot as plt
from switch import switch
import constants as c
import numpy

runs = c.RUNS
gens = c.numGenerations
fitnesses = numpy.zeros([runs, gens])
temp = []
individuals = []
with open('savedRobotsAfpoSeed.dat', "rb") as f:
    for r in range(1, runs+1):
        for g in range(1, gens+1):
            try:
                if g == gens:
                    individuals.append(pickle.load(f))
                    temp.append(individuals[0].fitness)
                else:
                    temp.append(pickle.load(f).fitness)
            except EOFError:
                break
        fitnesses[r-1] = temp
        temp = []
f.close()

mean_f = numpy.mean(fitnesses, axis=0)
std_f = numpy.std(fitnesses, axis=0)

plt.figure(figsize=(6.4,4.8))
plt.plot(list(range(1, gens+1)), mean_f, color='blue')
plt.fill_between(list(range(1, gens+1)), mean_f-std_f, mean_f+std_f, color='cornflowerblue', alpha=0.2)
plt.xlabel("Generations")
plt.ylabel("Best Fitness")
plt.title("Fitness of the Best Individual in the Population - AFPO", fontsize='small')
plt.grid(color='skyblue', linestyle=':', linewidth=0.5)
plt.tight_layout()
#plt.legend(['two robot', 'three robots'], loc='upper left')
plt.show()
plt.savefig("Afpo-bandgap.jpg", dpi=300)


# running the best individuals

bests = numpy.zeros([runs, gens])
temp = []
rubish = []
with open('savedRobotsLastGenAfpoSeed.dat', "rb") as f:
    for r in range(1, runs+1):
        # population of the last generation
        temp = pickle.load(f)
        # best individual of last generation
        best = temp[0]
        switch.showPacking(best.indv.genome)
        print("genome is:")
        print(best.indv.genome)
        print("fitness is:")
        print(switch.evaluate(best.indv.genome, True))
        temp = []
f.close()
