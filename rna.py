#!/usr/bin/env python3
import math
import RNA
from numpy.random import default_rng
from matplotlib import pyplot as plt
import pandas

TEMPLATE = ".((.((....))))."
MRATE = 0.05
BASES = ["A", "U", "G", "C"]
LENGTH = len(TEMPLATE)
POPULATION_SIZE = 50
STARTER_SIZE = int(POPULATION_SIZE / 5)
NUMBER_OF_SIMULATIONS = 10

RNG = default_rng()


def remove(collection, item):
    out = []
    for i in collection:
        if i == item:
            continue
        out.append(i)
    return out


def mutate(seq):
    newSeq = ""
    for base in seq:
        r = RNG.random()
        if r < MRATE:
            newSeq += RNG.choice(remove(BASES, base))
        else:
            newSeq += base
    return newSeq


def getPopulation(size, starterPop = []):
    population = []

    if len(starterPop) != 0:
        index = 0
        while index < size:
            seq = mutate(starterPop[index % len(starterPop)])
            population.append(seq)
            index += 1
    else:
        while len(population) < size:
            seq = ""
            for index in range(len(TEMPLATE)):
                seq += RNG.choice(BASES)
            population.append(seq)
    return population


def getStructure(seq):
    (ss, mfe) = RNA.fold(seq)
    return ss


def getDistances(structure):

    distance = RNA.string_edit_distance(
        RNA.Make_swString(TEMPLATE), RNA.Make_swString(structure))
    return (math.e**-distance)

def recombine(seq1, seq2):
    loc = RNG.integers(0,int(LENGTH*2/3))
    length = RNG.integers(0,int(LENGTH/3))

    seq1_out = seq1[:loc] + seq2[loc:loc+length] + seq1[loc+length:]
    seq2_out = seq2[:loc] + seq1[loc:loc+length] + seq2[loc+length:]

    return (seq1_out, seq2_out)

def simulate(recombination = False):
    population = []
    starterPop = []
    structures = []
    distances = []
    index = 0

    while not(TEMPLATE in structures):
        population = getPopulation(POPULATION_SIZE, starterPop)
        structures = list(map(getStructure, population))
        distances = list(map(getDistances, structures))
        distances = list(map(lambda x: x/sum(distances), distances))
        starterPop = list(RNG.choice(population, STARTER_SIZE, False, distances))

        if recombination and len(starterPop) >= 2:
            (index1, seq1) = RNG.choice(list(enumerate(starterPop)))
            (index2, seq2) = RNG.choice(list(enumerate(starterPop)))

            (seq1, seq2) = recombine(seq1, seq2)
            
            starterPop[int(index1)] = seq1
            starterPop[int(index2)] = seq2

        index += 1
    
    return index


def main():

    with open("noRecombination.csv",'w') as f:
        for i in range(NUMBER_OF_SIMULATIONS):
            generations = simulate(False)
            f.write(f'{generations}\n')

    
    df = pandas.read_csv("noRecombination.csv", names= ["Generations"], )
    df.hist()
    plt.savefig('noRecombination.png')

    with open("withRecombination.csv",'w') as f:
        for i in range(NUMBER_OF_SIMULATIONS):
            generations = simulate(True)
            f.write(f'{generations}\n')

    
    df = pandas.read_csv("withRecombination.csv", names= ["Generations"], )
    df.hist()
    plt.savefig('withRecombination.png')


if __name__ == '__main__':
    main()
