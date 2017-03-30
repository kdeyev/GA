import random
import copy
import numpy
import os
#import subprocess
#import math

#import numpy as np
#import math
#import matplotlib.pyplot as plot
#import mpl_toolkits.mplot3d.axes3d as axes3d
#import global_font

def nmo_FMM (g, tt, win, fast):
    win_samp = int(win/g.dt/2)

    import model_FD
    g_nmo = model_FD.gather(win_samp*2+1, g.dt, g.dh, g.sPos(), g.rPos())

#    fast = False
    if fast:
        for i in range (g.ntr):
            j = int (tt[i]/g.dt)
            for k in range(-win_samp, win_samp+1):
                in_samp = j+k
                if in_samp < 0 or in_samp >= g.nt:
                    continue
	            
                out_samp = k + win_samp
                g_nmo.v[out_samp][i] = g.v[in_samp][i]
	
        return g_nmo
    else:	
        import scipy.interpolate        
        v_t = numpy.transpose(g.v)
        samples = [i for i in range(g.nt)]
	
        for i in range (g.ntr):
            j = tt[i]/g.dt
            interpolator = scipy.interpolate.interp1d(samples, v_t[i])
            for k in range(-win_samp, win_samp+1):
                in_samp = j+k           
                if in_samp < 0 or in_samp >= g.nt-1:
                    continue
	            
                out_samp = k + win_samp
                g_nmo.v[out_samp][i] = interpolator(in_samp)
	            
        return g_nmo
    
def calc_fit_FMM (g, tt, win, fast, image_path=None):
    g_nmo = nmo_FMM (g, tt, win, fast)
    if image_path != None:
        g_nmo.draw('NMO', figure_name = image_path + '_nmo.png', cmap = 'gray')
        
    sum_ = numpy.zeros((g_nmo.nt))
    sum_sq = numpy.zeros((g_nmo.nt))
    nn = numpy.zeros((g_nmo.nt))
    energy = 0
    for i in range (g_nmo.ntr):
        for k in range(g_nmo.nt):
            amp = g_nmo.v[k][i]
            energy += amp ** 2
            nn[k] += 1
            if amp == 0:
                continue
            
            sum_[k] += amp
            sum_sq[k] += amp**2

         
#    print ('sum_',sum_)
#    print ('sum_sq',sum_sq)
#    print ('nn',nn)
    semb = numpy.zeros((g_nmo.nt))
    for k in range(g_nmo.nt):
        if sum_sq[k] != 0:
            semb[k] = (sum_[k]**2)/(sum_sq[k]*nn[k])
            
#    print ('semb',semb)
    
    aver_semb = numpy.average(semb)
    return aver_semb, energy

def put_spike_FMM (g, tt):
    
    import model_FD
    g_spike = model_FD.gather(g.nt, g.dt, g.dh, g.sPos(), g.rPos())

    for i in range (g.ntr):
        in_samp = int (tt[i]/g.dt)
        if in_samp < 0 or in_samp >= (g.nt -1):
            continue
	           
#        print (i, in_samp)
        g_spike.v[in_samp][i] = 1
        g_spike.v[in_samp+1][i] = 1

    g_spike.norm_ampl = 1
	
    return g_spike


def calcTT_FMM_(vel, source_pos):
    import numpy as np
#    import pylab as plt
    import skfmm
    
    [six, siz] = vel.getIndex(source_pos[0], source_pos[1])
    six = int (six)
    siz = int (siz)
    
    assert (vel.dx == vel.dz)
    dx = vel.dx
      
    phi = np.ones((vel.nx, vel.nz))
    phi[six][siz] = -1
#    print ('phi', phi)
 
    speed = vel.v

#    dist = skfmm.distance(phi, dx)

    tt = skfmm.travel_time(phi, speed, dx)

#    if image_path != None:
#        vel.draw ('', image_path +'_model.png', min_ = 1500, max_=5000)
#        
#        dist_t = np.transpose(dist)
#        phi_t = np.transpose(phi)
#        speed_t = np.transpose(speed) 
#        tt_t = np.transpose(tt)
#            
#        X, Z = np.meshgrid(np.linspace(-1,1,vel.nx), np.linspace(0,1,vel.nz))
##        print ("X", X)
##        print ("Z", Z)
#        
#        plt.subplot(221)
#        plt.title("Zero-contour of phi")
#        plt.contour(X, Z, phi_t, [0], colors='black', linewidths=(3))
#        #plt.gca().set_aspect(1)
#        plt.xticks([]); plt.yticks([])
#    
#            
#        plt.subplot(222)
#        plt.title("Distance")
#        plt.contour(X, Z, phi_t, [0], colors='black', linewidths=(3))
#        plt.contour(X, Z, dist_t , 15)
#        #plt.gca().set_aspect(1)
#        plt.xticks([]); plt.yticks([])
#        
#        plt.subplot(223)
#        plt.title("Vel")
#        plt.contour(X, Z, phi_t, [0], colors='black', linewidths=(3))
#        plt.contour(X, Z, speed_t , 15)
#        #plt.gca().set_aspect(1)
#        plt.xticks([]); plt.yticks([])
#        
#        plt.subplot(224)
##        print (X[0])
##        print (Z[0])
##        print (tt[0])
#        plt.title("Travel time")
#        plt.contour(X, Z, phi_t, [0], colors='black', linewidths=(3))
#        plt.contour(X, Z, tt_t, 15)
#        #plt.gca().set_aspect(1)
#        plt.xticks([]); plt.yticks([])
#    
##        plt.show()
#        plt.savefig(image_path + '_fmm.png',bbox_inches='tight', dpi=100)
    
    m_tt = copy.deepcopy(vel)
    m_tt.v = tt
    return m_tt
    
def calcTT_FMM (g, vel, fast):
    # take source position
    [sx, sz] = g.sPos()
    
    #calc fast marching travel times
    m_tt = calcTT_FMM_(vel, [sx, sz])
    
    if fast:
        rec_pos = g.rPos()
        tt = []
        for r in rec_pos:
            t = m_tt.getValue(r[0], r[1])
            tt.append(t)
	    
        return tt
    else:
        interpolator = m_tt.getInterp()
	        
        rec_pos = g.rPos()
        tt = []
        for r in rec_pos:
	#        print(r)
           t = interpolator(r)
           tt.append(t)
	    
        return tt

def weighted_choice_(w, power = 1, remove_average = 0):
  
  average = 0
  min_val = w[0]
  for i in range(len(w)):
      w[i] = max(0, w[i])
      min_val = min (min_val, w[i])
      average += w[i]

  average = average/len(w)

#  print ('average', average)
#  print ('min_val', min_val)
  
  if remove_average == 1:
      for i in range(len(w)):
          w[i] = w[i] - min_val

  if remove_average == 2:
      for i in range(len(w)):
          w[i] = w[i] - average

  for i in range(len(w)):
      w[i] = max(0, w[i])


  w = [v**power for v in w]
      
  weight_total = sum((v for v in w))
  
#  l = []
#  for item, weight in items:
#    l.append(weight/weight_total*100)
#  print ('distrib', l)
  
  n = random.uniform(0, weight_total)
  
  for i in range(len(w)):
    weight = w[i]
    if n < weight:
      return i
    n = n - weight
    
  print ("weighted_choice problem")
  exit (0)
  return -99
  
    
def weighted_choice(items, power = 1, remove_average = 0):
  """
  Chooses a random element from items, where items is a list of tuples in
  the form (item, weight). weight determines the probability of choosing its
  respective item. Note: this function is borrowed from ActiveState Recipes.
  """
  w = [item[1] for item in items]
  num = weighted_choice_(w, power, remove_average)
  return items[num][0]

def bound_random(v, dv, v1, v2):
    rv = 100000
    if v == None:
        v = int((v1 + v2)/2)
    if dv == None:
        dv = int((v1 + v2)/2)
        
#    print ("bound_random", "v", v, "rv", rv, "v1", v1, "v2", v2, "dv", dv)
    while v + rv not in range(v1, v2) : 
        rv = random.randrange(-dv, dv, 1)
        
    return v+rv
    
def bound_random_new(v1, v2):
    return random.randrange(v1, v2, 1)
  
def round_v (v1, dv1 = 100):
    return round(v1/dv1)*dv1

def round_z (v1, dv1 = 5):
    return round(v1/dv1)*dv1
        
def random_v(v1 = None):
    return round_v(bound_random (v1, 3000, 1500, 6000))    

class GA_constraint ():
    def applyConstraints (self, dna):
        print ('apply_constraint not implemented')
        
       
class GA_helper ():   
    
    def init(self, c, gathers, m, win, fast):
        self.c = c
        self.gathers = gathers
        self.m = m
        self.win = win
        self.fast = fast
        self.draw_gathers = False
        self._constraints = []
        self._cache = {}
        self._gatherCache = {}

        self.max_min = -1 # maximization or minimization
        
        
    def addConstraint(self, constraint):
        self._constraints.append(constraint)
        
    def applyConstraints(self, dna):
        for constraint in self._constraints:
            dna = constraint.applyConstraint(dna)
        return dna
        
    def define_FMM_energy(self):
        self.rt_energy_semb = 0

    def define_FMM_semb(self):
        self.rt_energy_semb = 1 
        
    def define_FMM_energy_semb(self):
        self.rt_energy_semb = 2 
        
    def print_info (self):
        print ('Fast Marching Method')
        if self.rt_energy_semb == 0:
            self.max_min = 0
            print ('Energy')
        if self.rt_energy_semb == 1:
            print ('Semblance')
            self.max_min = 0
        if self.rt_energy_semb == 2:
            print ('Energy * Semblance')
            self.max_min = 0
                
        if self.max_min == 0:
            print ('Maximization')
        if self.max_min == 1:
            print ('Minimization')

    def _random_dna(self):
        print ('random_dna not implemented')            

    def random_dna(self):
        dna = self._random_dna()
        dna = self.applyConstraints(dna)
        return dna
        
    def random_population(self, pop_size):
        pop = []
        for i in range(pop_size):       
            dna = self.random_dna()
            pop.append(dna)
        return pop
        
    def getTT_FMM(self, g, dna_m):
        tt = calcTT_FMM(g, dna_m, self.fast)   
#        print (tt)
        return tt   
        
    def putSpikeToGaters (self, dna):
        dna_m = self.getModel_FD(dna)            
        for shot in range(len(self.gathers)):
            g = self.gathers[shot]
            tt = self.getTT_FMM(g, dna_m);
#            gather_image_path = None
#            if image_path != None and self.draw_gathers:
#                gather_image_path = image_path + 'gather_' + str(shot)
            self.gathers[shot] = put_spike_FMM (g, tt) 
         

    def fitness(self, dna, image_path=None):
        key = str(dna)
        if self._cache.has_key (key):
#            print ("bingo!")
            return self._cache.get (key)
        
        fit = self.fitness_(dna, image_path)
        self._cache [key] = fit
        return fit
        
    def fitness_(self, dna, image_path=None):
        fd_info = None
        fitness = self.__fitness_FMM(dna, image_path)
            
        return fitness, fd_info
        
    def fitnessGather(self, dna, dna_m, shot, image_path=None):
        return self.__fitness_FMM_gather(dna, dna_m, shot, image_path)
        
    def __fitness_FMM_gather(self, dna, dna_m, shot, image_path=None):
        g = self.gathers[shot]

        gather_key = self.getGatherCache (shot, dna)
        if gather_key != None and self._gatherCache.has_key (gather_key):
            return self._gatherCache.get(gather_key)
            
        gather_image_path = None
        if image_path != None and self.draw_gathers:
            gather_image_path = image_path + 'gather_' + str(shot)
        tt = self.getTT_FMM(g, dna_m);
        semb, energy = calc_fit_FMM (g, tt,self.win, self.fast, gather_image_path) 
        gather_fitness = None
        if self.rt_energy_semb == 0:
            gather_fitness = energy
        if self.rt_energy_semb == 1:
            gather_fitness = semb
        if self.rt_energy_semb == 2:
            gather_fitness = semb * energy

        self._gatherCache[gather_key] = gather_fitness
                
        return gather_fitness
 

    def __fitness_FMM(self, dna, image_path=None):
        fitness = 0
        dna_m = self.getModel_FD(dna)            
        for shot in range(len(self.gathers)):
            fitness += self.__fitness_FMM_gather (dna, dna_m, shot, image_path)
                
        return fitness
            
    def getModel(self, dna):
        return self.getModel_FD(dna)            
            
    def getModel_FD(self, dna):
        print ('getModel not implemented')
    
    def _mutate(self, dna, mutation_chance):        
        print ('mutate not implemented')
        
    def mutate(self, dna, mutation_chance):        
        dna = self._mutate(dna, mutation_chance)
        dna = self.applyConstraints(dna)
        return dna
    
    def weight (self, population):
        if self.max_min == 0:
            return self.__weight_max (population)
        if self.max_min == 1:
            return self.__weight_min (population)
        
    def __weight_max (self, population):
        weighted_population = []

        # Add individuals and their respective fitness levels to the weighted
        # population list. This will be used to pull out individuals via certain
        # probabilities during the selection phase. Then, reset the population list
        # so we can repopulate it after selection.
        for individual in population:
            fitness_val, info = self.fitness(individual)

            # Generate the (individual,fitness) pair, taking in account whether or
            # not we will accidently divide by zero.
            pair = (individual, fitness_val, info)

            weighted_population.append(pair)
        return weighted_population
        
    def __weight_min (self, population):
        weighted_population = []

        # Add individuals and their respective fitness levels to the weighted
        # population list. This will be used to pull out individuals via certain
        # probabilities during the selection phase. Then, reset the population list
        # so we can repopulate it after selection.
        for individual in population:
            max_weight = 100000000000000000
            fitness_val, info = self.fitness(individual)
            
            
            # not we will accidently divide by zero.
            if fitness_val == 0:
                pair = (individual, max_weight, info)
            else:
                assert (1.0/fitness_val<max_weight)
                pair = (individual, 1.0/fitness_val, info)
 
            weighted_population.append(pair)
        return weighted_population
    
        
    def getBest(self, weighted_population):
        (best_ind, maximum_weight, best_fitness) = weighted_population[0]
                
        for i in range(len(weighted_population)):
            if weighted_population[i][1] >= maximum_weight:
                (best_ind, maximum_weight, best_fitness) = weighted_population[i]
                
        return (best_ind, maximum_weight, best_fitness)
        
#    def print_weight (self, weighted_population) :      
#        for i in range(len(weighted_population)):
#            print (weighted_population[i])
#            
#        (best_ind, maximum_weight, best_fitness) = self.getBest(weighted_population)
#        
#        weight_total = sum((item[1] for item in weighted_population))
#
#        print ("Best individual", best_ind, maximum_weight, best_fitness)
#        print ("Total fitness", weight_total)

    def draw (self, individual, images_path):
        #print ('individual',individual)
#        dna_m = self.getModel_FD(individual)
        fitness, info = self.fitness (individual, image_path=images_path)
        
        dna_m = self.getModel_FD(individual)
        dna_m.draw ('', images_path +'_model.png', min_ = 1500, max_=5000)
        
        if self.draw_gathers:
            for shot in range(len(self.gathers)):
                g = self.gathers[shot]
                tt = self.getTT_FMM(g, dna_m)
                if tt!= None:
                    gather_image_path = images_path + 'gather_' + str(shot)
                    g.draw (tt = tt, figure_name = gather_image_path + '.png')
                    
        return fitness, info 
  

    def _crossover(self, dna1, dna2):
        print ('_crossover not implemented')
        
    def crossover(self, dna1, dna2):  
        childs = self._crossover(dna1, dna2)
        for i in range(len(childs)):
            childs[i] = self.applyConstraints(childs[i])
        return childs
           
class GA_helperI1 (GA_helper):   

    def __init__(self, c, gathers, m, win, fast):
        self.init(c,gathers,m, win,fast)
        
        self.start_v1 = 1500
        self.end_v1 = 5000
        self.dv1 = 10
        self.start_v2 = 1500
        self.end_v2 = 5000
        self.dv2 = 10

        self.start_z = 50
        self.end_z = 250
        self.dz = 5
        
        self.nv1 = int((self.end_v1 - self.start_v1)/self.dv1 + 1)
        self.nv2 = int((self.end_v2 - self.start_v2)/self.dv2 + 1)
        self.nz = int((self.end_z - self.start_z)/self.dz + 1)
        
        self.cube = numpy.zeros((self.nz, self.nv2, self.nv1))
        
        self.z = numpy.arange(self.start_z, self.end_z + self.dz, self.dz)
        self.v1 = numpy.arange(self.start_v1, self.end_v1 + self.dv1, self.dv1)
        self.v2 = numpy.arange(self.start_v2, self.end_v2 + self.dv2, self.dv2)
        
        
    def random_v1(self,v1 = None):
        return round_v(bound_random (v1, 1000, self.start_v1, self.end_v1), self.dv1)
  
    def random_v2(self,v2 = None):
        return round_v(bound_random (v2, 1000, self.start_v2, self.end_v2), self.dv2)
        
    def random_z1(self,z1 = None):
        return round_z(bound_random (z1, 100, self.start_z, self.end_z), self.dz)
                       # WARNINNG min depth 50m

    def _random_dna(self):
        z1 = self.random_z1(125)
        v1 = self.random_v1(2000)
        v2 = self.random_v2(3500)
        
        dna = [z1, v1, v2]
        return dna

    ## PROXY
    def fitness(self, dna, image_path=None):            
        z = dna[0]
        v1 = dna[1]
        v2 = dna[2]
        z = int((z-self.start_z)/self.dz)
        v1 = int((v1-self.start_v1)/self.dv1)
        v2 = int((v2-self.start_v2)/self.dv2)

        # TODO: need  to uncomment     fitness_val for optimization
#        fitness_val = self.cube[z][v2][v1]
        fitness_val = 0
        
        info = []
        if fitness_val == 0:
            fitness_val, info = self.fitness_(dna, image_path)
            self.cube[z][v2][v1] = fitness_val


#        if image_path != None:
#            figure_name__ = image_path+'cube.png'
#            
#            plotcube(self.cube, 
#                     self.v1,
#                     self.v2,
#                     self.z,
#                      x_label = 'V1',
#                      y_label = 'V2',
#                      z_label = 'Z',
#                      figure_name = figure_name__
#                      )
#            
        return fitness_val, info
        
    @staticmethod
    def fillModel1 (m, z1, v1, v2):  
        for i in range (m.nx):
            for j in range (m.nz):
                z = j*m.dz 
                v = v1
                if (z > z1):
                    v = v2
                    
                m.v[i][j] = v
        
        return m    
    
    def getModel_FD(self, dna):
        import model_FD
        
        m = model_FD.model(self.c.nx, self.c.nz, self.c.dh, self.c.dh)
        dna_m = self.fillModel1(m, dna[0], dna[1], dna[2])
        return dna_m
            
    def _mutate(self, dna, mutation_chance):        
        for c in range(len(dna)):
            if random.random() <= mutation_chance:
                if c == 0:
#                    prev = dna[c]
                    dna[c] = self.random_z1(dna[c])
#                    print ('mutate z1', prev, dna[c])
                if c == 1:
#                    prev = dna[c]
                    dna[c] = self.random_v1(dna[c])
#                    print ('mutate v1', prev, dna[c])
                if c == 2:
#                    prev = dna[c]
                    dna[c] = self.random_v2(dna[c])
#                    print ('mutate v2', prev, dna[c])
        
        return dna
    
class GA_helperI2 (GA_helper):   
    
    def __init__(self, c, gathers, m, win, fast):
        self.init(c,gathers,m, win, fast)

    def empty_model (self):
        import model_FD
        lx = self.c.lx()
        lz = self.c.lz()
        new_dx = 10001
        new_dz = 25
        
        rand_m = model_FD.model(int(lx/new_dx)+1, int(lz/new_dz)+1, new_dx, new_dz)
        return rand_m
        
    def _random_dna(self):
        rand_m = self.empty_model()
        for i in range (rand_m.nx):
            for j in range (rand_m.nz):
                rand_m.v[i][j] = self.random_v(2500)
                    
        return numpy.reshape(rand_m.v,(rand_m.nx*rand_m.nz))
        
    @staticmethod
    def fillModel (m, dna_m):
#        to_interp = []
        for i in range (m.nx):
            for j in range (m.nz):
                x = i*m.dx
                z = j*m.dz
#                point=[x,z]
#                to_interp.append(point)
                m.v[i][j] = dna_m.getValue(x, z)

#        fn = dna_m.getInterp()                
#        m.v =  numpy.reshape(fn(to_interp), (m.nx, m.nz))

                
#        print (m.v)
        return m
    
    def getModel_FD(self, dna):
        dna_m = self.empty_model()
        dna_m.v = numpy.reshape(dna, (dna_m.nx, dna_m.nz))
#        print(dna_m.v)
#        print(dna_m.x_nodes)
#        print(dna_m.y_nodes)
#        print(dna_m.z_nodes)
        dna_m = self.fillModel(self.m, dna_m)
#        print(m.v)
#        print(m.x_nodes)
#        print(m.y_nodes)
#        print(m.z_nodes)
        return dna_m
            
    def _mutate(self, dna, mutation_chance):        
        for c in range(len(dna)):
            if random.random() <= mutation_chance:
                dna[c] = random_v(dna[c])
                
        return dna


class GA_helperI3 (GA_helper):   
    
    def __init__(self, c, gathers, m, win, fast):
        self.init(c,gathers,m, win, fast)
        self.layer_count = 3
        
    def random_z(self, z1 = None):
        lz = int(self.c.lz())
        if z1 == None:
            z1 = int(lz/2)
        return round_z(bound_random (z1, lz, 0, lz))    
        
    def sort_by_z(self, dna):
        for i in range(len(dna)):
            for j in range(i+1, len(dna)):
                if dna[i][0] > dna[j][0]:
                    tmp = dna[i]
                    dna[i] = dna[j]
                    dna[j] = tmp
        return dna
            
    def _random_dna(self):
        dna = [[0, random_v()]]
    
        for l in range(self.layer_count-1):
            dna.append([self.random_z(),random_v()])
        return self.sort_by_z(dna)


    def getModel_FD(self, dna):
        import model_FD
        m = model_FD.model(nx, nz, dx, dz) 
       
        interface_num = 0
        current_v = interfaces[interface_num][1]
        next_z = interfaces[interface_num+1][0]
        for j in range (m.nz):                
            z = j*dz 
            if z > next_z:
                interface_num += 1
                if interface_num < len (interfaces)-1:
                    next_z = interfaces[interface_num+1][0]
                else:
                    next_z = 10000
                    
                current_v = interfaces[interface_num][1]  
            for i in range (m.nx):
                m.v[i][j] = current_v
    
        return m    
            
    def _mutate(self, dna, mutation_chance):        
        for c in range(len(dna)):   
             if random.random() <= mutation_chance:
                 if c != 0: # WARNING special case
                    dna[c][0] = self.random_z(dna[c][0])
             if random.random() <= mutation_chance:
                 dna[c][1] = random_v(dna[c][1])                
               
        dna = self.sort_by_z(dna)
#        print ('mutate dna', dna)
        return dna

    def crossover2(self, dna1, dna2):
        child1 = []
        child2 = []
        for c in range(len(dna1)):
            child_dna1 = []
            child_dna2 = []
            w = random.random()
            if w < 0.5:
                child_dna1.append(dna1[c][0])
                child_dna2.append(dna2[c][0])
            else:
                child_dna1.append(dna2[c][0])
                child_dna2.append(dna1[c][0])
                
            w = random.random()
            if w < 0.5:
                child_dna1.append(dna1[c][1])
                child_dna2.append(dna2[c][1])
            else:
                child_dna1.append(dna2[c][1])
                child_dna2.append(dna1[c][1])
                
            child1.append(child_dna1)
            child2.append(child_dna2)
            
        child1 = self.sort_by_z(child1)
        child2 = self.sort_by_z(child2)
#        print ('cross parents', dna1, dna2)
#        print ('cross childs', child1, child2)
        return child1, child2
        
    def _crossover(self, dna1, dna2):
        dna1, dna2 = self.crossover2(dna1, dna2)
        dna1 = self.sort_by_z(dna1)
        dna2 = self.sort_by_z(dna2)
#        print ('crossover dna1', dna1)
#        print ('crossover dna2', dna2)
        
        return [dna1, dna2]
        
class model_FMM:   
    def __init__(self, nlayer, nx, dx, sx=0, sz=0):
        self.nx = nx
        self.nlayer = nlayer
        self.dx = dx
        self.sx = sx
        self.sz = sz
        self.v = numpy.zeros ((nlayer, nx))
        self.th = numpy.zeros ((nlayer,nx))
        self.x_nodes = [i*self.dx + self.sx for i in range (self.nx)]

    def prepareInterpolators (self):
        import scipy.interpolate
        self.interpolators_th = []
        self.interpolators_v = []
        for i in range (self.nlayer):  
#            print ('x_nodes', self.x_nodes)
#            print ('v', self.v[i])
#            print ('th', self.th[i])
            th = scipy.interpolate.interp1d(self.x_nodes, self.th[i], kind='cubic')
            v = scipy.interpolate.interp1d(self.x_nodes, self.v[i], kind='cubic')
            self.interpolators_th.append(th)
            self.interpolators_v.append(v)
            
    def generateFDModel (self, m):         
        self.prepareInterpolators()
       
        vv = []
        thth = []
        for i in range (self.nlayer):  
            v = numpy.reshape(self.interpolators_v[i](m.x_nodes),(m.nx))
#            print ('v', v)
            vv.append(v)
            
            th = numpy.reshape(self.interpolators_th[i](m.x_nodes),(m.nx))
#            print ('th', th)
            thth.append(th)
        
        for i in range(m.nx):
            layer_num = 0
            
            current_v = vv[layer_num][i]
#            print ('current_v', current_v)
            next_z = self.sz
            th = thth[layer_num][i]
#            print ('th', th)
            next_z = next_z + th
#            print ('next_z', next_z)
            
            for j in range (m.nz):                
                z = j*m.dz
#                print ('z', z)
                if z > next_z:
                    layer_num += 1
                    th = 100000
                    if layer_num < self.nlayer-1:
                        th = thth[layer_num][i]
#                        print ('th', th)
                    next_z = next_z + th
                        
                    current_v = vv[layer_num][i]
#                    print ('current_v', current_v)
                    
                m.v[i][j] = current_v
    
#        for i in range(m.nx): 
#            for j in range (m.nz):  
#                print (i*m.dx,j*m.dz ,m.v[i][j])
                
        return m    
        
class GA_helperI4 (GA_helper):   
    
    def __init__(self, c, gathers, m, win, fast, nlayer, nx, velConstr, thickConstr ):
        self.init(c, gathers, m, win, fast) 
        lx = self.c.lx()
        dx = lx/(nx-1)
        self.fmm_model = model_FMM(nlayer, nx, dx)
        
        self._velConstr = velConstr
        self._thickConstr = thickConstr

#        print (self.fmm_model.dx)
        self._gatherModelIndex = []
        for shot in range(len(self.gathers)):
           g = self.gathers[shot]   
           min_x = min ([r[0] for r in g.rPos ()])   
           max_x = max ([r[0] for r in g.rPos ()])                 
           min_ind = round((min_x - self.fmm_model.sx)/(self.fmm_model.dx))
           max_ind = round((max_x - self.fmm_model.sx)/(self.fmm_model.dx))
#           self._gatherModelIndex.append((min_ind, max_ind))
           
           gather_indices = []
           for j in range(self.fmm_model.nx):
               if j < min_ind or j > max_ind:
                   continue
               gather_indices.append (j)
           self._gatherModelIndex.append (gather_indices)
                   
#           print ("shot", shot, "index", self._gatherModelIndex[shot] )
        
        
    def random_th(self, th1, layer):
        return round_z(bound_random_new (self._thickConstr[layer][0], self._thickConstr[layer][1]))    
        
    def random_v_constr(self, v, layer):
        return round_z(bound_random_new (self._velConstr[layer][0], self._velConstr[layer][1]))  
        
    def _random_dna(self):
        dna = []
        for i in range(self.fmm_model.nlayer):
            layer = []
            for j in range(self.fmm_model.nx):
                layer.append([self.random_th(None, i),self.random_v_constr(None, i)])
            dna.append(layer)
        return dna

    def getModel_FD(self, dna):
#        print ('dna',dna)
        for i in range(self.fmm_model.nlayer):
            for j in range(self.fmm_model.nx):
                th = dna[i][j][0]
                v = dna[i][j][1]

                # TODO copy fmm_model!!!!
                self.fmm_model.v[i][j] = v
                self.fmm_model.th[i][j] = th

#        print ('v',self.fmm_model.v)
#        print ('th',self.fmm_model.th)

        import model_FD
        dna_m = model_FD.model(self.c.nx, self.c.nz, self.c.dh, self.c.dh) 
        dna_m = self.fmm_model.generateFDModel(dna_m)
#        print ('dna_m.v',dna_m.v)
        return dna_m
            
    def _mutate(self, dna, mutation_chance):      
        lz = int(self.c.lz())
        
        for i in range(self.fmm_model.nlayer):
            for j in range(self.fmm_model.nx):
                 if random.random() <= mutation_chance:
                    new_th = self.random_th(dna[i][j][0], i)
                    delta_th = new_th - dna[i][j][0]

#                    print ("_mutate", "th", dna[i][j][0], "new wth",new_th, "delta_th", delta_th )
                    dna[i][j][0] = dna[i][j][0] + delta_th
                    dna[i][j][0] = max (dna[i][j][0], 0)
                    dna[i][j][0] = min (dna[i][j][0], lz)
                    
                    # update layer below
                    if i < (self.fmm_model.nlayer -1):
                        dna[i+1][j][0] = dna[i+1][j][0] - delta_th
                        dna[i+1][j][0] = max (dna[i+1][j][0], 0)
                        dna[i+1][j][0] = min (dna[i+1][j][0], lz)

                 if random.random() <= mutation_chance:
                     dna[i][j][1] = self.random_v_constr(dna[i][j][1], i)                
               
        return dna
        
    def calcFitnessFunc (self, dna):
        dna_m = self.getModel_FD(dna)  
        
        fitness_func = numpy.zeros(self.fmm_model.nx)     
        fitness_gather = []
        for shot in range(len(self.gathers)):
            fitness = self.fitnessGather (dna, dna_m, shot)
            fitness_gather.append (fitness)
            
            gather_indices = self.getGatherIndices (shot)
            for j in gather_indices:
                fitness_func [j] += fitness
            
        return fitness_func, fitness_gather

    def checkChild (self, dna1, dna2, child):
#        fitness_func_child, fintess_gather_child = self.calcFitnessFunc (child)
#        for j in range(self.fmm_model.nx):
#            if fitness_func_child[j] < fitness_func1[j] or fitness_func_child[j] < fitness_func2[j]:
#                print ("PROBLEM function")
#                print ("fitness_func_child", fitness_func_child)
#                print ("fitness_func1", fitness_func1)
#                print ("fitness_func2", fitness_func2)
#                print ("mixed_index", mixed_index)
#                exit (1)
        
        fitness_child = self.fitness(child)[0]
        fitness1 = self.fitness(dna1)[0]
        fitness2 = self.fitness(dna2)[0]
        if fitness_child < fitness1 or fitness_child < fitness2:
#            print ("PROBLEM", fitness_child)
            
#            print ("fitness_child", fitness_child)
#            print ("fitness1", fitness1)
#            print ("fitness2", fitness2)
#            print ("")
#            print ("fitness_func_child", fitness_func_child)
#            print ("fitness_func1", fitness_func1)
#            print ("fitness_func2", fitness_func2)
#            print ("")
#            print ("fitness_func_child", sum(fitness_func_child))
#            print ("fitness_func1", sum(fitness_func1))
#            print ("fitness_func2", sum(fitness_func2))
#            print ("")
#            print ("fintess_gather_child", fintess_gather_child)
#            print ("fitness_gather1", fitness_gather1)
#            print ("fitness_gather2", fitness_gather2)
#            print ("")
#            print ("mixed_index", mixed_index)
#            exit (1)
#            throw (42)
            if fitness1 > fitness2:
                return [dna1]
            else:
                return [dna2]
        else:
#            print ("NORM", fitness_child)
            return [child]      
        
    def crossover3(self, dna1, dna2):
        fitness_func1, fitness_gather1 = self.calcFitnessFunc (dna1)
        fitness_func2, fitness_gather2 = self.calcFitnessFunc (dna2)
#        print ("fitness_func1", fitness_func1)
#        print ("fitness_func2", fitness_func2)
        
        child = copy.deepcopy(dna1)
        mixed_index = []
        for j in range(self.fmm_model.nx):
            parent = None
            if fitness_func1[j] >= fitness_func2[j]:
                parent = 1
            else:
                parent = 2

            mixed_index.append(parent)
            
            for i in range(self.fmm_model.nlayer):                        
                for c in range(len(dna1[i][j])):
                    if parent == 1:
                        child[i][j][c] = dna1[i][j][c]
                    else:
                        child[i][j][c] = dna2[i][j][c]

        return self.checkChild (dna1,dna2,child)

    def crossover4(self, dna1, dna2):
        fitness_func1, fitness_gather1 = self.calcFitnessFunc (dna1)
        fitness_func2, fitness_gather2 = self.calcFitnessFunc (dna2)
#        print ("fitness_func1", fitness_func1)
#        print ("fitness_func2", fitness_func2)
        
        child = copy.deepcopy(dna1)
#        mixed_index = []
        for j in range(self.fmm_model.nx):
            for i in range(self.fmm_model.nlayer):                        
                for c in range(len(dna1[i][j])):            
                    parent = weighted_choice_([fitness_func1[j], fitness_func2[j]])
                    if parent == 0:
                        child[i][j][c] = dna1[i][j][c]
                    if parent == 1:
                        child[i][j][c] = dna2[i][j][c]

        return self.checkChild (dna1,dna2,child)
        
    
    def crossoverPolygam(self, population, power = 1, remove_average = 1):
        
        fitness_func = [ self.calcFitnessFunc (dna) for dna in population]
        
        new_population = []
#        mixed_index = []
        
        while len (new_population) < pop_size:
            child = copy.deepcopy(population[0])
            for j in range(self.fmm_model.nx):
                for i in range(self.fmm_model.nlayer):                        
                    parent = weighted_choice_(fitness_func, power, remove_average)
                    for c in range(len(child[i][j])):            
                        child[i][j][c] = population[parent][i][j][c]

            new_population.append (child)

        return new_population


    def crossover2(self, dna1, dna2):
        child1 = copy.deepcopy(dna1)
        child2 = copy.deepcopy(dna2)
        for i in range(self.fmm_model.nlayer):
            for j in range(self.fmm_model.nx):
                for c in range(len(dna1[i][j])):
                    w = random.random()
                    if w < 0.5:
                        child1[i][j][c] = dna1[i][j][c]
                        child2[i][j][c] = dna2[i][j][c]
                    else:
                        child1[i][j][c] = dna2[i][j][c]
                        child2[i][j][c] = dna1[i][j][c]

        return [child1, child2]
        
    def crossover1(self, dna1, dna2):
        child1 = copy.deepcopy(dna1)
        child2 = copy.deepcopy(dna2)
        for i in range(self.fmm_model.nlayer):
            for j in range(self.fmm_model.nx):
                w = random.random()
                if w < 0.5:
                    child1[i][j] = dna1[i][j]
                    child2[i][j] = dna2[i][j]
                else:
                    child1[i][j] = dna2[i][j]
                    child2[i][j] = dna1[i][j]

        return [child1, child2]
        
    def _crossover(self, dna1, dna2):
        return self.crossover3(dna1, dna2)
                
    def getGatherIndices (self, shot):
        return self._gatherModelIndex [shot]

    def getGatherCache (self, shot, dna):
        gather_indices = self.getGatherIndices (shot)
        submodel = []
        for j in gather_indices:
            a = []
            for i in range(self.fmm_model.nlayer):
                a.append(dna[i][j])
            submodel.append(a)
#        print ("shot", shot, "submodel", submodel )
        return str(shot) + str(submodel)

def GA_test (helper, dna, pop_size, mutation = 0.1):
    correct = helper.fitness(dna)[0]
    correct = helper.applyConstraints(correct)
    print ('Correct answer:', correct)
    
    for i in range(pop_size):
        dna = helper.mutate(dna, mutation)
        
#        print(dna)
        fit = helper.fitness(dna)[0]
#        print (dna, fit)
#        print(i)
        if fit > correct:
            print ('We have a problem:', i, dna, fit)
            exit()
        
 
def MonteCarlo (helper, correct_dna, pop_size, mutation = 0.1):
    correct = helper.fitness(correct_dna)[0]
    print ('Correct answer:', correct)
    
    best_dna = helper.random_dna()
    
    best_fit = helper.fitness(best_dna)[0]
    
    for i in range(pop_size):
        
        if mutation > 0:
            new_dna = helper.mutate(best_dna, mutation)
        else:
            new_dna = helper.random_dna()
            
        new_fit = helper.fitness(new_dna)[0]
#        print (dna, fit)
        if new_fit > best_fit:
            best_dna = new_dna
            best_fit = new_fit
            print ('New best:', i, best_dna, best_fit)
            
#
# Main driver
# Generate a population and simulate GENERATIONS generations.
#

def create_new_population (helper, mutation, weighted_population):
               
    pop_size = len (weighted_population)            
    new_population = []
    # Select two random individuals, based on their fitness probabilites, cross
    # their genes over at a random point, mutate them, and add them back to the
    # population for the next iteration.
    while len (new_population) < pop_size:
        # Selection
        ind1 = weighted_choice(weighted_population)
        ind2 = weighted_choice(weighted_population)
#            print ("after weighted_choice")


        # Crossover
        childs = helper.crossover(ind1, ind2)
#            print (childs)
        for child in childs:
            # Mutate and add back into the population.
            child = helper.mutate(child, mutation)
            new_population.append(child)
#            print ("after crosover")
    return new_population
    

def create_new_population_polygam (helper, mutation, weighted_population):
    new_population = [v[0] for v in weighted_population]
    new_population = helper.crossoverPolygam(new_population)
    
    for i in range(len(new_population)):
        new_population[i] = helper.mutate(new_population[i], mutation)
        
    return new_population;


def GA_run_on_population (helper, images_path, population, 
        generatoin_count = 30, mutation = 0.1):  
    
    weighted_population = helper.weight (population)

    best_generation = 0
    (global_best_ind, global_maximum_weight, global_best_fitness) = helper.getBest(weighted_population)
    
    print ('Init')
#    helper.print_weight (weighted_population)
    print ("global best individual", global_best_ind, global_maximum_weight, global_best_fitness)
    
    # print start point
    helper.draw (global_best_ind, images_path + "start")

    # Simulate all of the generations.
    for generation in range(generatoin_count):
    
        population = create_new_population (helper, mutation, weighted_population)
            
        weighted_population = helper.weight (population)    
        (local_best_ind, local_maximum_weight, best_fitness) = helper.getBest(weighted_population)
        
#        print ("after getBest")
        
        if local_maximum_weight > global_maximum_weight:
            global_best_ind = local_best_ind
            global_maximum_weight = local_maximum_weight
            global_best_fitness = best_fitness
            best_generation = generation
            print ("new global best!")
            helper.draw (local_best_ind, images_path + 'iter_' + str(generation))
#            print ("after draw")

    
        print ('Generation', generation)
#        helper.print_weight (weighted_population)

        print ("global best individual", best_generation, global_maximum_weight, global_best_fitness)
        print ("local best individual", local_maximum_weight, best_fitness)
            
        weight_total = sum((item[1] for item in weighted_population))
        print ("Total weight", weight_total)

    return population
#        helper.fitness(local_best_ind,  image_path=images_path+'gen_'+str(generation))

def GA_run (helper, images_path, correct_dna,
        pop_size = 20, generatoin_count = 30, mutation = 0.1):  
    helper.print_info()

    if not os.path.exists(images_path):
        os.makedirs(images_path)
        
    correct = helper.fitness(correct_dna)
    print ('Correct answer:', correct)
    
    helper.draw (correct_dna, images_path + "correct")
    
    # Generate initial population. This will create a list of POP_SIZE strings,
    # each initialized to a sequence of random characters.
    population = helper.random_population(pop_size)
    population = GA_run_on_population(helper, images_path, population, generatoin_count, mutation)
    return population