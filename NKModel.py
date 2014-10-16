import random

class NK_Model(object):
    '''
    A Kauffman NK Model, with tuneable ruggedness. 



    '''
    
    def __init__(self, n, k=2, random_fitness=True, fitness_vector = None):
        '''
        Initialize a new NK Model
        '''
        self.n = n
        self.k = k
        self.random_fitness = random_fitness
        self.graph = None
        if self.random_fitness:
            # Build a random fitness dictionary
            self.locus_fitness = {}
            for i in range(2**self.k):
                locus = self._num_to_tuple(i, self.k)
                self.locus_fitness[locus] = random.randrange(-5, 5)

        elif len(fitness_vector) == self.k:
            #self.fitness_vector = fitness_vector
            self.locus_fitness = fitness_vector
        else:
            raise Exception("Fitness vector wrong length!")


    # ======================================
    #         Fitness evaluation
    # ======================================
    
    def _local_fitness(self, locus):
        '''
        Evaluate the local fitness of a locus
        '''
        if locus not in self.locus_fitness:
            raise Exception("Incorrect locus!")

        return self.locus_fitness[locus]
    
    def eval_fitness(self, genome):
        '''
        Evaluate the fitness of an entire genome.
        '''
        if len(genome) != self.n:
            raise Exception("Wrong length!")

        fitness = 0
        for i in range(self.n):
            locus = []
            for d in range(i-(self.k-1), i+1): # Loop from i-(k-1) to i, inclusive
                locus.append(genome[d])
            fitness += self._local_fitness(tuple(locus))
        return fitness

    # ======================================
    #         Optimization
    # ======================================
    @staticmethod
    def _num_to_tuple(num, length):
        '''
        Convert a number to a tuple of its binary representation.
        '''
        binary = bin(num)[2:]
        while len(binary) < length:
            binary = "0" + binary
        if len(binary) > length:
            raise Exception("Number too big to belong here!")
        binary = list(binary)
        bits = [int(b) for b in binary]
        return tuple(bits)
        
    def enumerate_space(self):
        '''
        Exhaustively enumerate the fitness of each point in the model space.
        '''
        model_space = {}
        for i in range(2**self.n):
            # Convert number to bit string
            genome = self._num_to_tuple(i, self.n)
            model_space[genome] = self.eval_fitness(genome)
        return model_space