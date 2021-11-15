import time
import os
import argparse
import copy
import sys
import numpy as np

class Chromosome:

    def __init__(
        self,
        FirstLowerBound,
        FirstUpperBound,
        SecondLowerBound,
        SecondUpperBound,
        preferred,
        fitness = 0
        ):
        self.FirstLowerBound = FirstLowerBound
        self.FirstUpperBound = FirstUpperBound
        self.SecondLowerBound = SecondLowerBound
        self.SecondUpperBound = SecondUpperBound
        self.preferred = preferred
        self.fitness = fitness

    #grabs the input file and calculates the fitness later
    def calcFitness(self, input_data):
        self.fitness = 0
        fitnesses = []
        for line in input_data:
            match = int(self.FirstLowerBound <= line[0] <= self.FirstUpperBound
                        and self.SecondLowerBound <= line[1]
                        <= self.SecondUpperBound)  # sets last num to 0 or 1 to see if it fits within the parameters
            profit = int(2 * (int((2 * self.preferred - 1)) * line[2] > 0) - 1)  
            # equals -1 or 1 depending on if a profit is made
            new_fitness = match * profit * round(abs(line[2]), 2)  # if not a match, new_fitness equals 0, if a profit is made, it's positive, if there was a loss, it's negative.
            self.fitness += new_fitness
            fitnesses.append(new_fitness)
        self.fitness = round(sum(fitnesses), 2)

    def __str__(self):
        return str(self.FirstLowerBound) + '\t' + str(self.FirstUpperBound) + '\t' + str(self.SecondLowerBound) + '\t' + str(self.SecondUpperBound) + '\t' + str(self.preferred) + '\tfitness: ' + str(self.fitness) + '\n'

def arg_parse():
    if len(sys.argv) != 8:
        raise ValueError("Enter the name of the file you'd like to read with the number of chromosomes, number of generations, your choice of elitist or tournament, the percentage of what should be formed using selection in this generation and finally the mutation rate.")
        #parses arguments that you need to put in, otherwise you will get the error above
    return (
        sys.argv[1],
        int(sys.argv[2]),
        int(sys.argv[3]),
        sys.argv[4],
        float(sys.argv[5]),
        sys.argv[6],
        float(sys.argv[7]),
        )

def RandomNum():
    return round(np.random.normal(loc=0.0, scale=1.15, size=None), 2)
    #gets a random number to fill in the values of the chromosome

def initialGeneration(chromosome_num):
    #generates the initial generation and seeds in the random number from earlier to split it out of the list of chromosomes that were generaters
    chromosomes = []
    np.random.seed(int(time.time()))
    for i in range(chromosome_num):
        random_val1 = [RandomNum(), RandomNum()]
        random_val1.sort()
        random_val2 = [RandomNum(), RandomNum()]
        random_val2.sort()
        #this will give each chromosome their own assigned value
        new_chromosomes = Chromosome(FirstLowerBound = random_val1[0],
                           FirstUpperBound = random_val1[1],
                           SecondLowerBound = random_val2[0],
                           SecondUpperBound = random_val2[1],
                           preferred = np.random.randint(0, 2))
        chromosomes.append(new_chromosomes)
              #adds the new chromosome to a new list and then return chromosomes
    return chromosomes

#this is for the input file, it will read the file and will return the input file data as an array
def input_data(input_file):
    data = []
    file = open(input_file)
    GenAlgData = file.readlines()
    for line in GenAlgData:
        clean_line = line.split('\n')[0].split('\t')
        for i in range(len(clean_line)):
            clean_line[i] = float(clean_line[i])
        data.append(clean_line)
    return data


def new_offspring(parents, uni_or_kpoint):
  #this function returns the children of the chromosomes, and it will display an error if the input does not include uniform or kpoint
    entry = uni_or_kpoint.lower()
    if entry != 'uniform' and entry != 'kpoint':
      #make the decision between kpoint and uniform, displays an error when input is not lowercase or spelt correctly as seen in the if statment parameters
        raise ValueError("Enter the name of the file you'd like to read with the number of chromosomes, number of generations, your choice of elitist or tournament, the percentage of what should be formed using selection in this generation and finally the mutation rate."
                         )
        sys.exit()
    if uni_or_kpoint.lower() == 'uniform': #uniform function
        new_child = Chromosome(parents[np.random.randint(0, 2)].FirstLowerBound,
                          parents[np.random.randint(0, 2)].FirstUpperBound,
                          parents[np.random.randint(0, 2)].SecondLowerBound,
                          parents[np.random.randint(0, 2)].SecondUpperBound,
                          parents[np.random.randint(0, 2)].preferred)
    else: #this is the kpoint function

        new_child = Chromosome(parents[0].FirstLowerBound,
                          parents[0].FirstUpperBound,
                          parents[1].SecondLowerBound,
                          parents[1].SecondUpperBound, 
                          parents[1].preferred)

    return new_child

def mutation(chromosomes, mutation_rate):
    for i in range(len(chromosomes)):
        Bottom = np.random.uniform(0, 1) <= mutation_rate
        Top = np.random.uniform(0, 1) <= mutation_rate
        Bottom1 = np.random.uniform(0, 1) <= mutation_rate
        Top1 = np.random.uniform(0, 1) <= mutation_rate
        reacclimate = np.random.uniform(0, 1) <= mutation_rate

        #takes the mutation rate and list of chromoosmes, output is a new list which displays the mutation chromosomes, which is entirely based off of the user input in the command      
        if Bottom:
            chromosomes[i].FirstLowerBound = int(Bottom) * RandomNum() + \
                 int(not Bottom) * chromosomes[i].FirstLowerBound
            while chromosomes[i].FirstLowerBound >= chromosomes[i].FirstUpperBound:
                chromosomes[i].FirstLowerBound = int(Bottom) * RandomNum() + \
                     int(not Bottom) * chromosomes[i].FirstLowerBound

        if Top:
            chromosomes[i].FirstUpperBound = int(Top) * RandomNum() + \
                  int(not Bottom) * chromosomes[i].FirstUpperBound
            while chromosomes[i].FirstUpperBound <= chromosomes[i].FirstLowerBound:
                chromosomes[i].FirstUpperBound = int(Top) * RandomNum() + \
                     int(not Top) * chromosomes[i].FirstUpperBound

        if Bottom1:
            chromosomes[i].SecondLowerBound = int(Bottom1) * RandomNum() + \
                 int(not Bottom1) * chromosomes[i].SecondLowerBound
            while chromosomes[i].SecondLowerBound >= chromosomes[i].SecondUpperBound:
                chromosomes[i].SecondLowerBound = int(Bottom1) * RandomNum() + \
                     int(not Bottom1) * chromosomes[i].SecondLowerBound

        if Top1:
            chromosomes[i].FirstUpperBound = int(Top1) * RandomNum() + \
                 int(not Top1) * chromosomes[i].SecondUpperBound
            while chromosomes[i].SecondUpperBound <= chromosomes[i].SecondLowerBound:
                chromosomes[i].SecondUpperBound = int(Top1) * RandomNum() + \
                     int(not Top1) * chromosomes[i].SecondUpperBound

        chromosomes[i].preferred = int(reacclimate) * np.random.randint(0, 2) + \
             int(not  reacclimate) * chromosomes[i].preferred
    return chromosomes


def fitness(chromosomes, input_data):
  #assign the fitness calclations value for the chromosomes 
    for chrom in chromosomes:
        chrom.calcFitness(input_data)
    return chromosomes

def fitnessfunction(n):
  #sorts the chromosomes
    return n.fitness

def selection(chromosomes, elite_or_tourney, x):
  #selects the fittest chrom and check the users input to see if they want to do elitist or tournament and for any possible errors, it has to be in lowercase 
    entry = elite_or_tourney.lower()
    best_fit_chromosomes = []
    if entry != 'elitist' and entry != 'tournament':
        raise ValueError("Enter the name of the file you'd like to read with the number of chromosomes, number of generations, your choice of elitist or tournament, the percentage of what should be formed using selection in this generation and finally the mutation rate."
                         )

        sys.exit()
    if entry == 'elitist':
        chromosomes.sort(reverse=True, key=fitnessfunction)
        best_fit_chromosomes = chromosomes[0:x]
    else:

        while(len(best_fit_chromosomes) < x):
            chromosome1 = chromosomes[np.random.randint(0, len(chromosomes))]
            chromosome2 = chromosomes[np.random.randint(0, len(chromosomes))]
            best_fit_chromosomes.append(max(chromosome1, chromosome2, key=fitnessfunction))
    print("Fittest Chromosomes: ")
    print_chromosomes(best_fit_chromosomes)
    return best_fit_chromosomes


def print_chromosomes(chromosomes):
  #print the list of chromosomes
    fitnesses = []
    for chrom in chromosomes:
        fitnesses.append(chrom.fitness)
        print(chrom)
    print('\n')
    print('max: ' + str(max(fitnesses)))
    print('min: ' + str(min(fitnesses)))
    print('mean:' + str(round(sum(fitnesses) / len(fitnesses), 2)))


def new_gen(half_chromosome, chromosome_num, uni_or_kpoint):
  #new generation of chromosomes based off of the old generation
    chromosomes = copy.deepcopy(half_chromosome)
    while len(chromosomes) < chromosome_num:
        parents = [half_chromosome[np.random.randint(0, len(half_chromosome))],
                   half_chromosome[np.random.randint(0, len(half_chromosome))]]
        new_child = new_offspring(parents, uni_or_kpoint)
        chromosomes.append(new_child)
    return chromosomes


def main():
    input_file, chromosome_num, num_generations, elitist_or_tournament, selection_rate, uni_or_kpoint,  mutation_rate = arg_parse()
    data = input_data(input_file)
    chromosomes = initialGeneration(chromosome_num)
    num_produced = int(float(chromosome_num*selection_rate))
    chromosomes = fitness(chromosomes, data)
    print_chromosomes(chromosomes)

    for i in range(num_generations):
        print('Generation ' + str(i + 1) + ': ')
        chromosomes = selection(chromosomes, elitist_or_tournament, num_produced)
        #selection function used to determine what would be used for new generation
        chromosomes = new_gen(chromosomes, chromosome_num, uni_or_kpoint)
        #makes new chromosomes based off the previous chromosome generation 
        chromosomes = mutation(chromosomes, mutation_rate)
        #mutation function in order to mutate the genes that needed to be mutated, based off of the user input
        chromosomes = fitness(chromosomes, data)
        #fitness function to determine the fitness of each chromosome
        print_chromosomes(chromosomes)

if __name__ == '__main__':
    main()
