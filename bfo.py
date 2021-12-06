import numpy as np
from random import random


class cell(object):

    def __init__(self):

        self.cell_pos = []
        self.best_pos = []

    def set_best_pos(self, best):
      """ set the current best cell position """
      self.best_pos = best

    def add_cells(self, cells):

        self.cell_pos.append([list(i) for i in cells])

    def get_cells(self):
        """returns all the bacterial cells"""
        return self.cell_pos

    def get_best_pos(self):
        """return best"""

        return list(self.best_pos)


class bfo(cell):

    def __init__(self,function,S, bounds, n, Nre,
                 Nc=2, Ns=12, Ci=0.2, Ped=1.15):
        """
        S: number of bacterium
        n: dimension of the search space
        function: math function
        bounds: [lb ,ub]
        Nre: reproductive steps
        Nc: chemotactic 
        Ns: swim steps
        Ci: the chemotactic step size during each run or tumble.
        Ped: probability of elimination
        """

        super(bfo, self).__init__()
        lb, ub = bounds
        #  Initialize randomly the bacteria foraging optimization population
        self.cells = np.random.uniform(lb, ub, (S, n))
        self.add_cells(self.cells)

        #  Calculate the fitness of each cell
        fitness = np.array([function(x) for x in self.cells])
        #  Set global best cell to best cell
        new_best = self.cells[fitness.argmin()]
        best = new_best

        C_list = [Ci - Ci * 0.9 * i / Nre for i in range(Nre)]
        Ped_list = [Ped - Ped * 0.5 * i / Nre for i in range(Nre)]

        las_fitness = fitness[::1]
        # for the reproductive steps
        for t in range(Nre):

            chem_fitness = [fitness[::1]]
            # for the chemotactic steps range
            for j in range(Nc):
              # for each  cell
                for i in range(S):
                    # For this dimension, move in that direction by the step distance from
                    # where the cell currently is
                    dell = np.random.uniform(-1, 1, n)
                    self.cells[i] += C_list[t] * np.linalg.norm(dell) * dell
                    # calc the fitness of the moved cell
                    c_fitness = function(self.cells[i])
                    # for swimming range
                    for m in range(Ns):
                        if c_fitness < las_fitness[i]:
                            las_fitness[i] = fitness[i]
                            self.cells[i] += C_list[t] * np.linalg.norm(dell) \
                                                * dell
                        else:# mmove cell to random direction
                            dell = np.random.uniform(-1, 1, n)
                            self.cells[i] += C_list[t] * np.linalg.norm(dell) \
                                                * dell
                # calc the fitness of each cell
                fitness = np.array([function(x) for x in self.cells])
                chem_fitness += [fitness]

            chem_fitness = np.array(chem_fitness)

            # calculatefitness function of all chemotactic loops 
            hcells = [(sum(chem_fitness[:, i]), i) for i in range(S)]
            hcells.sort()

            # split the cells
            acells = []
            for i in hcells:
                acells += [list(self.cells[i[1]])]

            n_is_even = True if n & 1 else False

            if n_is_even:
                acells = 2*acells[:S//2]
                self.cells = np.array(acells)
            else:
                acells = 2*acells[:S//2] +\
                                [acells[S//2]]
                self.cells = np.array(acells)

            # if not the last reproductive steps
            if t < Nre - 2:
                # for each  cell
                for i in range(S):
                    r = random()
                    #replace cell with new random generated cells with some P
                    if r >= Ped_list[t]:
                        self.cells[i] = np.random.uniform(lb, ub, n)

            fitness = np.array([function(x) for x in self.cells])
            self.add_cells(self.cells)

            # Calculate the fitness of each cell
            new_best = self.cells[fitness.argmin()]
            # print(new_best)
            if function(new_best) < function(best):
                # Update the best  cell
                best = new_best

        self.set_best_pos(best)


