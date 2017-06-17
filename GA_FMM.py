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


def getSparkContext ():
#    #######
#    # Spark
#    #######
#    from pyspark import SparkContext, SparkConf
#    
#    conf = SparkConf()
#    conf.setAppName("GA spark")
#    conf.set("spark.ui.enabled", "false" )
#    conf.setMaster("local[15]")
##        if seisspark_config.local_spark:
##    conf.setMaster("local")
#    sc = SparkContext(conf=conf)
##        s_sc.addPyFile("seisspark_config.py")
##        s_sc.addPyFile("seisspark.py")
##        s_sc.addPyFile("segypy.py")
#    return sc
    return None

    
g_sc = getSparkContext ()

def writeArray (filename, dna):
    with open(filename, 'wb') as f:
#            dna = numpy.reshape (dna, (self.fmmModel.nlayer*self.fmmModel.nx*2))
        numpy.save (f, dna)
        
def readArray (filename):
    with open(filename, 'rb') as f:
        dna = numpy.load (f)
#            dna = numpy.reshape (dna, (self.fmmModel.nlayer, self.fmmModel.nx, 2))
        return dna
            
            
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

def put_spike_FMM (g, tt, win=20):
    
    import model_FD
    g_spike = model_FD.gather(g.nt, g.dt, g.dh, g.sPos(), g.rPos())

    for i in range (g.ntr):
        in_samp = int (tt[i]/g.dt)
        if in_samp < 0 or in_samp >= (g.nt):
            continue
	           
        for k in range (max(in_samp-win,0),min(in_samp+win+1,g.nt-1)):
            g_spike.v[k][i] = 1.-abs(float(k-in_samp))/float(win)

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

  return len(w)-1
  print ("weighted_choice problem")
  exit (0)
  return -99
  
    
def weighted_choice(population, power = 1, remove_average = 0):
  """
  Chooses a random element from items, where items is a list of tuples in
  the form (item, weight). weight determines the probability of choosing its
  respective item. Note: this function is borrowed from ActiveState Recipes.
  """
  w = [individ.fitness for individ in population]
  num = weighted_choice_(w, power, remove_average)
  return population[num]

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
        
       
    
class Individ ():
    def __init__ (self, dna):
        self.dna = dna
        self.addToCache = True
        self.key = str(dna)
        self.fitness = None
        self.fitness_PS = 0.0
        self.v = copy.deepcopy(dna)
        for i in range(len(dna)):
            for j in range(len(dna[0])):
                for c in range(len(dna[0][0])):
                    self.v[i][j][c] = 0
        
        self.gather_individs = []
        self.fitness_func = None
        
class GatherIndivid ():
    def __init__ (self, shot, key=None, fitness=None):
        self.shot = shot
        self.addToCache = True
        self.key = key
        self.fitness = fitness
        
class GA_helper ():   
    
    def init(self, modelGeom, gathers, win, fast):
        self.modelGeom = modelGeom
        self.gathers = gathers
        self.win = win
        self.fast = fast
        self.draw_gathers = False
        self._constraints = []
        self._cache = {}
        self._gatherCache = {}
        self._use_cache = False
        
        
    def addConstraint(self, constraint):
        self._constraints.append(constraint)
        
    def createIndivid (self, dna):
        individ = Individ(dna)
        individ.gather_individs = []
        for shot in range(len(self.gathers)):
            key = self.getGatherCacheKey (shot, dna)
            individ.gather_individs.append(GatherIndivid(shot=shot, key=key))
        return individ
        
    def applyConstraints(self, individ):
#        dna = copy.deepcopy(individ.dna)
        for constraint in self._constraints:
            individ.dna = constraint.applyConstraint(individ.dna)
            
        return individ
#        if individ.dna == dna:
##            print ("applyConstraints: individ was not changed")
#            return individ
#        else:
##            print ("applyConstraints: individ was changed")
#            return self.createIndivid(dna)
        
    def define_FMM_energy(self):
        self.rt_energy_semb = 0

    def define_FMM_semb(self):
        self.rt_energy_semb = 1 
        
    def define_FMM_energy_semb(self):
        self.rt_energy_semb = 2 
        
    def print_info (self):
        print ('Fast Marching Method')
        if self.rt_energy_semb == 0:
            print ('Energy')
        if self.rt_energy_semb == 1:
            print ('Semblance')
        if self.rt_energy_semb == 2:
            print ('Energy * Semblance')
                
        print ('Maximization')

    def _random_individ(self):
        print ('random_individ not implemented')            

    def random_individ(self):
        individ = self.createIndivid(self._random_individ())
        individ = self.applyConstraints(individ)
        return individ
        
    def random_population(self, pop_size):
        pop = []
        for i in range(pop_size):       
            individ = self.random_individ()
#            individ = self.fitness(individ)
            pop.append(individ)
        return pop
        
    def getTT_FMM(self, g, dna_m):
        tt = calcTT_FMM(g, dna_m, self.fast)   
#        print (tt)
        return tt   
        
    def putSpikeToGathers (self, dna):
        dna_m = self.getModel_FD(dna)            
        for shot in range(len(self.gathers)):
            g = self.gathers[shot]
            tt = self.getTT_FMM(g, dna_m);
#            gather_image_path = None
#            if image_path != None and self.draw_gathers:
#                gather_image_path = image_path + 'gather_' + str(shot)
            self.gathers[shot] = put_spike_FMM (g, tt) 
            
    def fitness(self, individ, image_path=None):
#        if individ.fitness != None:
##            print ("individ: individ already calculated")
#            individ.addToCache=False
#            return individ
        
        individ.addToCache=True
#        print ("individ: individ not calculated")
        
        dna_m = self.getModel_FD(individ.dna)   
        individ.fitness = 0
        for j in range(len(individ.gather_individs)):
            individ.gather_individs[j] = self.fitnessCalcGather(dna_m, individ.gather_individs[j], image_path)
            individ.fitness += individ.gather_individs[j].fitness
            
        return individ

    def populateCache (self, population):  
        cache_hit = 0
        gather_cache_hit = 0
        for i in range(len(population)):
            individ = population[i]
            if individ.key == None:
                throw (42)

            if individ.addToCache==True and self._use_cache:
                self._cache [individ.key] = individ
                individ.addToCache=False
            else:
                cache_hit = cache_hit + 1
        
            for j in range(len(individ.gather_individs)):
                gather_individ = individ.gather_individs[j]

                if gather_individ.key == None:
                    throw (42)
                
                if gather_individ.addToCache==True and self._use_cache:
                    self._gatherCache [gather_individ.key] = gather_individ
                    gather_individ.addToCache=False
                    individ.gather_individs[j] = gather_individ
                else:
                    gather_cache_hit = gather_cache_hit + 1      
                    
            population[i] = individ 
                
        print ("populateCache: cache_hit", cache_hit, "gather_cache_hit", gather_cache_hit)
        return population

    def fillFromCache (self, population):  
        cache_hit = 0
        gather_cache_hit = 0
        
#        print ('self._cache.keys', self._cache.keys ())
        for i in range(len(population)):
            individ = population[i]
            if individ.key == None:
                throw (42)
                
#            print ("fillFromCache: individ", individ.fitness)
            if individ.fitness==None and self._cache.has_key (individ.key):
                individ = self._cache.get(individ.key)
                cache_hit = cache_hit + 1               
        
                
            for j in range(len(individ.gather_individs)):
                gather_individ = individ.gather_individs[j]
#                print ("fillFromCache: gather_individ", gather_individ.fitness)
                if gather_individ.key == None:
                    throw (42)
                    
                if gather_individ.fitness==None and self._gatherCache.has_key (gather_individ.key):
                    gather_individ = self._gatherCache.get(gather_individ.key)
                    gather_cache_hit = gather_cache_hit + 1               
                    individ.gather_individs[j] = gather_individ

            population[i] = individ
                
        print ("fillFromCache: cache_hit", cache_hit, "gather_cache_hit", gather_cache_hit)
        return population
        
        
#    def fitnessCalc(self, individ, image_path=None):
#        individ.fitness = 0
#        individ.gather_individs = []
#        dna_m = self.getModel_FD(individ)            
#        for shot in range(len(self.gathers)):
#            gather_individ = self.fitnessGather (individ, dna_m, shot, image_path)
#            individ.gather_individs.append (gather_individ)
#            individ.fitness += gather_individ.fitness 
#                
#        return individ
        
#    def fitnessGatherCache(self, dna, gatherCache, dna_m, shot, image_path=None):
#        gather_key = self.getGatherCacheKey (shot, dna)
#        gather_individs = None
#        addToCache = gatherCache.has_key (gather_key)
#        
#        if addToCache:
#            gather_individs = gatherCache.get(gather_key)
#        else:            
#            gather_individs = self.fitnessCalcGather(dna_m, shot, image_path)
#            
#        return GatherIndivid(shot, addToCache, gather_key, gather_individs)
        
        
    def fitnessCalcGather(self, dna_m, gather_individ, image_path=None):
#        if gather_individ.fitness != None:
#            gather_individ.addToCache = False
#            return gather_individ
#            
        g = self.gathers[gather_individ.shot]

        gather_individ.addToCache = True
        gather_image_path = None
        if image_path != None and self.draw_gathers:
            gather_image_path = image_path + 'gather_' + str(gather_individ.shot)
        tt = self.getTT_FMM(g, dna_m);
        semb, energy = calc_fit_FMM (g, tt,self.win, self.fast, gather_image_path) 
        gather_individ.fitness  = None
        if self.rt_energy_semb == 0:
            gather_individ.fitness  = energy
        if self.rt_energy_semb == 1:
            gather_individ.fitness  = semb
        if self.rt_energy_semb == 2:
            gather_individ.fitness  = semb * energy

        return gather_individ
 
    def getModel(self, dna):
        return self.getModel_FD(dna)            
            
    def getModel_FD(self, dna):
        print ('getModel not implemented')
    
    def _mutate(self, individ, mutation_chance):        
        print ('mutate not implemented')
        
    def mutate(self, individ, mutation_chance):        
        individ = self._mutate(individ, mutation_chance)
        individ = self.applyConstraints(individ)
        return individ
    
    def weight (self, population):
        
        print ("cache size ", len(self._cache))
        print ("gather cache size ", len(self._gatherCache))
        
        population = self.fillFromCache (population)

        # Add individuals and their respective fitness levels to the weighted
        # population list. This will be used to pull out individuals via certain
        # probabilities during the selection phase. Then, reset the population list
        # so we can repopulate it after selection.
            
        if g_sc != None:
            par_pop = g_sc.parallelize (population)       
            population = par_pop.map(lambda individ: self.g_broadcastHelper.value.fitness(individ)).collect()
        else:   
            population = [self.fitness(individ) for individ in population]
                              
        population = self.populateCache (population)
        return population 
        
    def getBest(self, population):
        best_individ = population[0]
                
        for i in range(len(population)):
            if population[i].fitness >= best_individ.fitness:
                best_individ = population[i]
                
        return best_individ
        
#    def print_weight (self, population) :      
#        for i in range(len(population)):
#            print (population[i])
#            
#        (best_ind, maximum_weight, best_fitness) = self.getBest(population)
#        
#        weight_total = sum((item[1] for item in population))
#
#        print ("Best individual", best_ind, maximum_weight, best_fitness)
#        print ("Total fitness", weight_total)

    def draw (self, individ, images_path):
        #print ('individual',individual)
#        dna_m = self.getModel_FD(individual)
        individ = self.fitness(individ, image_path=images_path)
        
        dna_m = self.getModel_FD(individ.dna)
        dna_m.draw ('', images_path +'_model.png', min_ = 1500, max_=5000)
        
        if self.draw_gathers:
            for shot in range(len(self.gathers)):
                g = self.gathers[shot]
                tt = self.getTT_FMM(g, dna_m)
                if tt!= None:
                    gather_image_path = images_path + 'gather_' + str(shot)
                    g.draw (tt = tt, figure_name = gather_image_path + '.png')
                    
        return individ 
  

    def _crossover(self, individ1, individ2):
        print ('_crossover not implemented')
        
    def crossover(self, individ1, individ2):  
        childs = self._crossover(individ1, individ2)
        for i in range(len(childs)):
            childs[i] = self.applyConstraints(childs[i])
        return childs
#           
#class GA_helperI1 (GA_helper):   
#
#    def __init__(self, modelGeom, gathers, win, fast):
#        self.init(modelGeom, gathers, win,fast)
#        
#        self.start_v1 = 1500
#        self.end_v1 = 5000
#        self.dv1 = 10
#        self.start_v2 = 1500
#        self.end_v2 = 5000
#        self.dv2 = 10
#
#        self.start_z = 50
#        self.end_z = 250
#        self.dz = 5
#        
#        self.nv1 = int((self.end_v1 - self.start_v1)/self.dv1 + 1)
#        self.nv2 = int((self.end_v2 - self.start_v2)/self.dv2 + 1)
#        self.nz = int((self.end_z - self.start_z)/self.dz + 1)
#        
##        self.cube = numpy.zeros((self.nz, self.nv2, self.nv1))
#        
#        self.z = numpy.arange(self.start_z, self.end_z + self.dz, self.dz)
#        self.v1 = numpy.arange(self.start_v1, self.end_v1 + self.dv1, self.dv1)
#        self.v2 = numpy.arange(self.start_v2, self.end_v2 + self.dv2, self.dv2)
#        
#        
#    def random_v1(self,v1 = None):
#        return round_v(bound_random (v1, 1000, self.start_v1, self.end_v1), self.dv1)
#  
#    def random_v2(self,v2 = None):
#        return round_v(bound_random (v2, 1000, self.start_v2, self.end_v2), self.dv2)
#        
#    def random_z1(self,z1 = None):
#        return round_z(bound_random (z1, 100, self.start_z, self.end_z), self.dz)
#                       # WARNINNG min depth 50m
#
#    def _random_individ(self):
#        z1 = self.random_z1(125)
#        v1 = self.random_v1(2000)
#        v2 = self.random_v2(3500)
#        
#        dna = [z1, v1, v2]
#        return dna
#      
##    def fitness(self, dna, image_path=None):            
##        key = str(dna)
##        if self._cache.has_key (key):
###            print ("bingo!")
##            return self._cache.get (key)
##        
##        fitness = self.fitnessCalc(dna, image_path)
##        self._cache [key] = fitness
##        return fitness
#            
#    @staticmethod
#    def fillModel1 (m, z1, v1, v2):  
#        for i in range (m.nx):
#            for j in range (m.nz):
#                z = j*m.dz 
#                v = v1
#                if (z > z1):
#                    v = v2
#                    
#                m.v[i][j] = v
#        
#        return m    
#    
#    def getModel_FD(self, dna):
#        import model_FD
#        
#        m = model_FD.model(self.modelGeom.nx, self.modelGeom.nz, self.modelGeom.dx, self.modelGeom.dz)
#        dna_m = self.fillModel1(m, dna[0], dna[1], dna[2])
#        return dna_m
#            
#    def _mutate(self, individ, mutation_chance):        
#        dna = copy.deep_copy
#        for c in range(len(dna)):
#            if random.random() <= mutation_chance:
#                if c == 0:
##                    prev = dna[c]
#                    dna[c] = self.random_z1(dna[c])
##                    print ('mutate z1', prev, dna[c])
#                if c == 1:
##                    prev = dna[c]
#                    dna[c] = self.random_v1(dna[c])
##                    print ('mutate v1', prev, dna[c])
#                if c == 2:
##                    prev = dna[c]
#                    dna[c] = self.random_v2(dna[c])
##                    print ('mutate v2', prev, dna[c])
#        
#        return dna
#    
#class GA_helperI2 (GA_helper):   
#    
#    def __init__(self, modelGeom, gathers, win, fast):
#        self.init(modelGeom,gathers, win, fast)
#
#    def empty_model (self):
#        import model_FD
#        lx = self.modelGeom.lx()
#        lz = self.modelGeom.lz()
#        new_dx = 10001
#        new_dz = 25
#        
#        rand_m = model_FD.model(int(lx/new_dx)+1, int(lz/new_dz)+1, new_dx, new_dz)
#        return rand_m
#        
#    def _random_individ(self):
#        rand_m = self.empty_model()
#        for i in range (rand_m.nx):
#            for j in range (rand_m.nz):
#                rand_m.v[i][j] = self.random_v(2500)
#                    
#        return numpy.reshape(rand_m.v,(rand_m.nx*rand_m.nz))
#        
#    @staticmethod
#    def fillModel (m, dna_m):
##        to_interp = []
#        for i in range (m.nx):
#            for j in range (m.nz):
#                x = i*m.dx
#                z = j*m.dz
##                point=[x,z]
##                to_interp.append(point)
#                m.v[i][j] = dna_m.getValue(x, z)
#
##        fn = dna_m.getInterp()                
##        m.v =  numpy.reshape(fn(to_interp), (m.nx, m.nz))
#
#                
##        print (m.v)
#        return m
#    
#    def getModel_FD(self, dna):
#        dna_m = self.empty_model()
#        dna_m.v = numpy.reshape(dna, (dna_m.nx, dna_m.nz))
##        print(dna_m.v)
##        print(dna_m.x_nodes)
##        print(dna_m.y_nodes)
##        print(dna_m.z_nodes)
#        dna_m = self.fillModel(self.m, dna_m)
##        print(m.v)
##        print(m.x_nodes)
##        print(m.y_nodes)
##        print(m.z_nodes)
#        return dna_m
#            
#    def _mutate(self, dna, mutation_chance):        
#        for c in range(len(dna)):
#            if random.random() <= mutation_chance:
#                dna[c] = random_v(dna[c])
#                
#        return dna
#
#
#class GA_helperI3 (GA_helper):   
#    
#    def __init__(self, modelGeom, gathers, win, fast):
#        self.init(modelGeom, gathers, win, fast)
#        self.layer_count = 3
#        
#    def random_z(self, z1 = None):
#        lz = int(self.modelGeom.lz())
#        if z1 == None:
#            z1 = int(lz/2)
#        return round_z(bound_random (z1, lz, 0, lz))    
#        
#    def sort_by_z(self, dna):
#        for i in range(len(dna)):
#            for j in range(i+1, len(dna)):
#                if dna[i][0] > dna[j][0]:
#                    tmp = dna[i]
#                    dna[i] = dna[j]
#                    dna[j] = tmp
#        return dna
#            
#    def _random_individ(self):
#        dna = [[0, random_v()]]
#    
#        for l in range(self.layer_count-1):
#            dna.append([self.random_z(),random_v()])
#        return self.sort_by_z(dna)
#
#
##    def getModel_FD(self, dna):
##        import model_FD
##        m = model_FD.model(self.modelGeom.nx, self.modelGeom.nz, self.modelGeom.dx, self.modelGeom.dz) 
##       
##        interface_num = 0
##        current_v = interfaces[interface_num][1]
##        next_z = interfaces[interface_num+1][0]
##        for j in range (m.nz):                
##            z = j*dz 
##            if z > next_z:
##                interface_num += 1
##                if interface_num < len (interfaces)-1:
##                    next_z = interfaces[interface_num+1][0]
##                else:
##                    next_z = 10000
##                    
##                current_v = interfaces[interface_num][1]  
##            for i in range (m.nx):
##                m.v[i][j] = current_v
##    
##        return m    
#            
#    def _mutate(self, dna, mutation_chance):        
#        for c in range(len(dna)):   
#             if random.random() <= mutation_chance:
#                 if c != 0: # WARNING special case
#                    dna[c][0] = self.random_z(dna[c][0])
#             if random.random() <= mutation_chance:
#                 dna[c][1] = random_v(dna[c][1])                
#               
#        dna = self.sort_by_z(dna)
##        print ('mutate dna', dna)
#        return dna
#
#    def crossover2(self, individ1, individ2):
#        child1 = []
#        child2 = []
#        for c in range(len(individ1.dna)):
#            child_dna1 = []
#            child_dna2 = []
#            w = random.random()
#            if w < 0.5:
#                child_dna1.append(individ1.dna[c][0])
#                child_dna2.append(individ2.dna[c][0])
#            else:
#                child_dna1.append(individ2.dna[c][0])
#                child_dna2.append(individ1.dna[c][0])
#                
#            w = random.random()
#            if w < 0.5:
#                child_dna1.append(individ1.dna[c][1])
#                child_dna2.append(individ2.dna[c][1])
#            else:
#                child_dna1.append(individ2.dna[c][1])
#                child_dna2.append(individ1.dna[c][1])
#                
#            child1.append(child_dna1)
#            child2.append(child_dna2)
#            
#        child1 = self.sort_by_z(child1)
#        child2 = self.sort_by_z(child2)
##        print ('cross parents', individ1, individ2)
##        print ('cross childs', child1, child2)
#        return child1, child2
#        
#    def _crossover(self, individ1, individ2):
#        individ1, individ2 = self.crossover2(individ1, individ2)
#        dna1 = self.sort_by_z(individ1.dna)
#        dna2 = self.sort_by_z(individ2.dna)
##        print ('crossover dna1', dna1)
##        print ('crossover dna2', dna2)
#        
#        return [createIndivid(dna1), createIndivid(dna2)]

class GA_helperI4_Constraint_V (GA_constraint):
    def __init__(self, correct_dna):
        self.correct_dna = correct_dna
        self.nlayer = len(correct_dna)
        self.nx = len(correct_dna[0])
        self.constr_layer = 0
        
    def applyConstraint (self, dna):
        for i in range(self.nlayer):
            for j in range(self.nx):
                if i == self.constr_layer: # first layer
                    # use correct velocity
                    dna[i][j][1] = self.correct_dna[i][j][1]
        
                
        return dna
        
class GA_helperI4_Constraint_Well_Depth (GA_constraint):
    def __init__(self, correct_dna):
        self.correct_dna = correct_dna
        self.nlayer = len(correct_dna)
        self.nx = len(correct_dna[0])
        self.well_pos = int (self.nx/2)
        
    def applyConstraint (self, dna):
        for i in range(self.nlayer):
            for j in range(self.nx):
                if j == self.well_pos: 
                    # use correct depth
                    dna[i][j][0] = self.correct_dna[i][j][0]
        
                
        return dna
        

class GA_helperI4_Constraint_Well (GA_constraint):
    def __init__(self, correct_dna):
        self.correct_dna = correct_dna
        self.nlayer = len(correct_dna)
        self.nx = len(correct_dna[0])
        self.well_pos = int (self.nx/2)
        
    def applyConstraint (self, dna):
        for i in range(self.nlayer):
            for j in range(self.nx):
                if j == self.well_pos: 
                    # use correct depth
                    dna[i][j][0] = self.correct_dna[i][j][0]
                    # use correct vel
                    dna[i][j][1] = self.correct_dna[i][j][1]
        
                
        return dna
        
        
class GA_helperI4_Constraint_Point (GA_constraint):
    def __init__(self, correct_dna, v, d):
        self.correct_dna = correct_dna
        self.nlayer = len(correct_dna)
        self.nx = len(correct_dna[0])
        self.well_pos = int (self.nx/2)
        self.layer = 0
        self.v = v
        self.d = d
        
    def applyConstraint (self, dna):
        for i in range(self.nlayer):
            for j in range(self.nx):
                if j == self.well_pos and i == self.layer:
                    if self.d == True:
                        # use correct depth
                        dna[i][j][0] = self.correct_dna[i][j][0]
                    if self.v == True:
                        # use correct vel
                        dna[i][j][1] = self.correct_dna[i][j][1]
        
                
        return dna
        
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
            
    def generateFDModel (self, modelGeom):         
        self.prepareInterpolators()
        
        import model_FD
        dna_m = model_FD.model(modelGeom.nx, modelGeom.nz, modelGeom.dx, modelGeom.dz) 
       
        vv = []
        thth = []
        for i in range (self.nlayer):  
            v = numpy.reshape(self.interpolators_v[i](dna_m.x_nodes),(dna_m.nx))
#            print ('v', v)
            vv.append(v)
            
            th = numpy.reshape(self.interpolators_th[i](dna_m.x_nodes),(dna_m.nx))
#            print ('th', th)
            thth.append(th)
        
        for i in range(dna_m.nx):
            layer_num = 0
            
            current_v = vv[layer_num][i]
#            print ('current_v', current_v)
            next_z = self.sz
            th = thth[layer_num][i]
#            print ('th', th)
            next_z = next_z + th
#            print ('next_z', next_z)
            
            for j in range (dna_m.nz):                
                z = j*dna_m.dz
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
                    
                dna_m.v[i][j] = current_v
    
#        for i in range(m.nx): 
#            for j in range (m.nz):  
#                print (i*m.dx,j*m.dz ,m.v[i][j])
                
        return dna_m    
        
    def getIndexX (self, x):
        return round( (x - self.sx)/self.dx)
        
class GA_helperI4 (GA_helper):   
    
    def __init__(self, modelGeom, gathers, win, fast, nlayer, nx, velConstr, thickConstr ):
        self.init(modelGeom, gathers, win, fast) 
        lx = self.modelGeom.lx()
        dx = lx/(nx-1)
        self.fmmModel = model_FMM(nlayer, nx, dx)
        
        self._velConstr = velConstr
        self._thickConstr = thickConstr

#        print (self.fmmModel.dx)
        self._gatherModelIndex = []
        for shot in range(len(self.gathers)):
           g = self.gathers[shot]   
           min_x = min ([r[0] for r in g.rPos ()])   
           max_x = max ([r[0] for r in g.rPos ()])                 
           min_ind = self.fmmModel.getIndexX(min_x)
           max_ind = self.fmmModel.getIndexX(max_x)
#           self._gatherModelIndex.append((min_ind, max_ind))
           
           gather_indices = []
           for j in range(self.fmmModel.nx):
               if j < min_ind or j > max_ind:
                   continue
               gather_indices.append (j)
           self._gatherModelIndex.append (gather_indices)
                   
#           print ("shot", shot, "index", self._gatherModelIndex[shot] )
        
        
    def random_th(self, th1, layer):
        return round_z(bound_random_new (self._thickConstr[layer][0], self._thickConstr[layer][1]), self._thickConstr[layer][2])    
        
    def random_v_constr(self, v, layer):
        return round_v(bound_random_new (self._velConstr[layer][0], self._velConstr[layer][1]), self._velConstr[layer][2])  
        
    def _random_individ(self):
        dna = []
        for i in range(self.fmmModel.nlayer):
            layer = []
            for j in range(self.fmmModel.nx):
                layer.append([self.random_th(None, i),self.random_v_constr(None, i)])
            dna.append(layer)
        return dna

    def getModel_FD(self, dna):
        fmm_model = copy.deepcopy(self.fmmModel)
#        print ('dna',dna)
        for i in range(fmm_model.nlayer):
            for j in range(fmm_model.nx):
                th = dna[i][j][0]
                v = dna[i][j][1]

                fmm_model.v[i][j] = v
                fmm_model.th[i][j] = th

#        print ('v',self.fmmModel.v)
#        print ('th',self.fmmModel.th)

        dna_m = fmm_model.generateFDModel(self.modelGeom)
#        print ('dna_m.v',dna_m.v)
        return dna_m
        
    def caclCombinationNum (self, ):
        num = 1
        for i in range(self.fmmModel.nlayer):
            for j in range(self.fmmModel.nx):
                num = num * (self._thickConstr[i][1] - self._thickConstr[i][0])/self._thickConstr[i][2]
                num = num * (self._velConstr[i][1] - self._velConstr[i][0])/self._velConstr[i][2]
        return num
            
    def _mutate(self, individ, mutation_chance):      
        lz = int(self.modelGeom.lz())
        
        dna = copy.deepcopy(individ.dna)
        for i in range(self.fmmModel.nlayer):
            for j in range(self.fmmModel.nx):
                 if random.random() <= mutation_chance:
                    new_th = self.random_th(dna[i][j][0], i)
                    delta_th = new_th - dna[i][j][0]

#                    print ("_mutate", "th", dna[i][j][0], "new wth",new_th, "delta_th", delta_th )
                    dna[i][j][0] = dna[i][j][0] + delta_th
                    dna[i][j][0] = max (dna[i][j][0], 0)
                    dna[i][j][0] = min (dna[i][j][0], lz)
                    
                    # update layer below
                    if i < (self.fmmModel.nlayer -1):
                        dna[i+1][j][0] = dna[i+1][j][0] - delta_th
                        dna[i+1][j][0] = max (dna[i+1][j][0], 0)
                        dna[i+1][j][0] = min (dna[i+1][j][0], lz)

                 if random.random() <= mutation_chance:
                     dna[i][j][1] = self.random_v_constr(dna[i][j][1], i)                
               
        if individ.dna == dna:
#            print ("_mutate: individ was not changed")
            return individ
        else:
#            print ("_mutate: individ was changed")
            return self.createIndivid(dna)
        
    def calcFitnessFunc (self, individ):
        if individ.fitness_func != None:
            return individ
        
        individ = self.fitness (individ)
        
        individ.fitness_func = numpy.zeros(self.fmmModel.nx)     
        for gather_individ in individ.gather_individs:
            
            gather_indices = self.getGatherIndices (gather_individ.shot)
            for j in gather_indices:
                individ.fitness_func [j] += gather_individ.fitness
            
        return individ

        
    def drawCurves (self, population, last):
       
        images_path = 'C:\GA\tests\evgeny\GA_images_evgeny_FMM/'
       
        fff = [self.calcFitnessFunc (individ).fitness_func for individ in population]
        fitness_func = [v[0] for v in fff]
        fitness_gather = [v[1] for v in fff]
        
        import model_FD
        model_FD.draw_convergence (fitness_gather, "Gather", "Fitness", "Fitness of gathers", images_path + 'gather_fit.png', last = last) 
        model_FD.draw_convergence (fitness_func, "Position", "Fitness", "Fitness function of parent", images_path + 'func_fit.png')
        
        exit (0)

    def maxFitness (self, dna):
        dna = numpy.reshape (dna, (self.fmmModel.nlayer, self.fmmModel.nx, 2))
        return 1./self.fitness(self.createIndivid(dna)).fitness
        
    def maximize_nelder_mead (self, dna):
        dna = numpy.reshape (dna, (self.fmmModel.nlayer*self.fmmModel.nx*2))
        
        from scipy.optimize import minimize
        ret_val = minimize(self.maxFitness, dna, method='nelder-mead')
        dna = ret_val.get('x')
        dna = numpy.reshape (dna, (self.fmmModel.nlayer, self.fmmModel.nx, 2))
        return dna
        
    def maximize_basinhopping (self, dna):
        dna = numpy.reshape (dna, (self.fmmModel.nlayer*self.fmmModel.nx*2))
        
        from scipy.optimize import basinhopping
        ret_val = basinhopping(self.maxFitness, dna)
        dna = ret_val.get('x')
        dna = numpy.reshape (dna, (self.fmmModel.nlayer, self.fmmModel.nx, 2))
        return dna
        
        
    def maximize_BFGS (self, dna):
        dna = numpy.reshape (dna, (self.fmmModel.nlayer*self.fmmModel.nx*2))
        
        from scipy.optimize import minimize
        ret_val = minimize(self.maxFitness, dna, method='BFGS')
        dna = ret_val.get('x')
        dna = numpy.reshape (dna, (self.fmmModel.nlayer, self.fmmModel.nx, 2))
        return dna
        
    def get_bounds (self, dna = None) : 
        lb = []
        ub = []
        for i in range(self.fmmModel.nlayer):
            for j in range(self.fmmModel.nx):
                lb.append (self._thickConstr[i][0])
                ub.append (self._thickConstr[i][1])
                lb.append (self._velConstr[i][0])
                ub.append (self._velConstr[i][1])
        
        if dna != None:
            for i in range (len(dna)):
                if dna [i] < lb[i]:
                    dna [i] = lb[i] 
                if dna[i] > ub[i]:
                    dna [i] = ub[i] 
                
        return lb, ub
        
    def maximize_CMA (self, dna):
        dna = numpy.reshape (dna, (self.fmmModel.nlayer*self.fmmModel.nx*2))    
        lb, ub = self.get_bounds (dna)
                
        import cma
#        o = cma.CMAOptions()
#        for k in o:  # returns all possible options
#            print (k, o[k])
#        exit (0)
        options = {
                   'bounds' : [lb, ub], 
                   'maxiter' : 100,
                   'popsize' : 100}

        res = cma.fmin(self.maxFitness, dna, 0.5, options = options)
        dna = res[0]
        dna = numpy.reshape (dna, (self.fmmModel.nlayer, self.fmmModel.nx, 2))
        return dna

    def maximize_swarm (self, images_path):
#        dna = numpy.reshape (dna, (self.fmmModel.nlayer*self.fmmModel.nx*2))
        lb, ub = self.get_bounds ()
                
        from pyswarm import pso
        dna, fopt = pso(self.maxFitness, lb, ub, maxiter=1000, swarmsize = 100)
        dna = numpy.reshape (dna, (self.fmmModel.nlayer, self.fmmModel.nx, 2))
                
        individ = self.fitness(self.createIndivid (dna))
        self.draw_gathers = True
        self.draw (individ, images_path + "swarm_final")
        print ("swarm individual", individ.fitness)        
        
    

    def maximize_differential_evolution (self, images_path):
#        dna = numpy.reshape (dna, (self.fmmModel.nlayer*self.fmmModel.nx*2))
        bounds = []
        for i in range(self.fmmModel.nlayer):
            for j in range(self.fmmModel.nx):
                bounds.append ([self._thickConstr[i][0], self._thickConstr[i][1]])
                bounds.append ([self._velConstr[i][0], self._velConstr[i][1]])
                
        from scipy.optimize import differential_evolution
        ret_val = differential_evolution(self.maxFitness, bounds=bounds, maxiter=300, popsize=15, disp=True)
        dna = ret_val.get('x')
        dna = numpy.reshape (dna, (self.fmmModel.nlayer, self.fmmModel.nx, 2))
        
        
        individ = self.fitness(self.createIndivid (dna))
        self.draw_gathers = True
        self.draw (individ, images_path + "differential_evolution_final")
        print ("maximize_differential_evolutiont individual", individ.fitness)        
        
        return dna
        
        
    def maximize (self, dna):
        return self.maximize_nelder_mead (dna)
 
        
#    def checkChild (self, individ1, individ2, child):
##         fitness_func_child, fintess_gather_child = self.calcFitnessFunc (child)
##        for j in range(self.fmmModel.nx):
##            if fitness_func_child[j] < fitness_func1[j] or fitness_func_child[j] < fitness_func2[j]:
##                print ("PROBLEM function")
##                print ("fitness_func_child", fitness_func_child)
##                print ("fitness_func1", fitness_func1)
##                print ("fitness_func2", fitness_func2)
##                print ("mixed_index", mixed_index)
##                exit (1)
#        
#        fitness_child = self.fitness(child)
#        fitness1 = self.fitness(dna1)
#        fitness2 = self.fitness(dna2)
#        if fitness_child < fitness1 or fitness_child < fitness2:
##            print ("PROBLEM", fitness_child)
#            
##            print ("fitness_child", fitness_child)
##            print ("fitness1", fitness1)
##            print ("fitness2", fitness2)
##            print ("")
##            print ("fitness_func_child", fitness_func_child)
##            print ("fitness_func1", fitness_func1)
##            print ("fitness_func2", fitness_func2)
##            print ("")
##            print ("fitness_func_child", sum(fitness_func_child))
##            print ("fitness_func1", sum(fitness_func1))
##            print ("fitness_func2", sum(fitness_func2))
##            print ("")
##            print ("fintess_gather_child", fintess_gather_child)
##            print ("fitness_gather1", fitness_gather1)
##            print ("fitness_gather2", fitness_gather2)
##            print ("")
##            print ("mixed_index", mixed_index)
##            exit (1)
##            throw (42)
#            if fitness1 > fitness2:
#                return [dna1]
#            else:
#                return [dna2]
#        else:
##            print ("NORM", fitness_child)
#            return [child]      
        
    def crossover3(self, individ1, individ2):
        individ1 = self.calcFitnessFunc (individ1)
        individ2 = self.calcFitnessFunc (individ2)
#        print ("fitness_func1", fitness_func1)
#        print ("fitness_func2", fitness_func2)
        
        child_dna = copy.deepcopy(individ1.dna)
        mixed_index = []
        for j in range(self.fmmModel.nx):
            parent = None
            if individ1.fitness_func1[j] >= individ1.fitness_func2[j]:
                parent = 1
            else:
                parent = 2

            mixed_index.append(parent)
            
            for i in range(self.fmmModel.nlayer):                        
                for c in range(len(individ1.dna[i][j])):
                    if parent == 1:
                        child_dna[i][j][c] = individ1.dna[i][j][c]
                    else:
                        child_dna[i][j][c] = individ2.dna[i][j][c]

        return [self.createIndivid(child_dna)]    
#        return self.checkChild (dna1,dna2,child)

    def crossover4(self, individ1, individ2):
        individ1 = self.calcFitnessFunc (individ1)
        individ2 = self.calcFitnessFunc (individ2)
#        print ("fitness_func1", fitness_func1)
#        print ("fitness_func2", fitness_func2)
        
        child_dna = copy.deepcopy(individ1.dna)
#        mixed_index = []
        for j in range(self.fmmModel.nx):
            parent = weighted_choice_([individ1.fitness_func[j], individ2.fitness_func[j]])
            for i in range(self.fmmModel.nlayer):                        
                for c in range(len(individ1.dna[i][j])):            
                    if parent == 0:
                        child_dna[i][j][c] = individ1.dna[i][j][c]
                    if parent == 1:
                        child_dna[i][j][c] = individ2.dna[i][j][c]

        return [self.createIndivid(child_dna)]   
#        return self.checkChild (dna1,dna2,child)
        
    
    def crossoverPolygam(self, population, power = 1, remove_average = 0):
        
        population = [ self.calcFitnessFunc (individ) for individ in population]
#        print ("fitness_func", fitness_func[0])
                
        pop_size = len(population)
        
        fitness_func_transposed = []
        for j in range(self.fmmModel.nx):
            fitness_func_transposed.append([individ.fitness_func[j] for individ in population])
            
#        print ("fitness_func", fitness_func_transposed[0])
        
        new_population = []
        while len (new_population) < pop_size:
            child_dna = copy.deepcopy(population[0].dna)
            for j in range(self.fmmModel.nx):
                parent = weighted_choice_(fitness_func_transposed[j], power, remove_average)                                            
                for i in range(self.fmmModel.nlayer):  
                    for c in range(len(child_dna[i][j])):            
                        child_dna[i][j][c] = population[parent][i][j][c]

            new_population.append (self.createIndivid(child_dna))

        return new_population


    def crossover2(self, individ1, individ2):
        child_dna1 = copy.deepcopy(individ1.dna)
        child_dna2 = copy.deepcopy(individ2.dna)
        for i in range(self.fmmModel.nlayer):
            for j in range(self.fmmModel.nx):
                for c in range(len(individ1.dna[i][j])):
                    w = random.random()
                    if w < 0.5:
                        child_dna1[i][j][c] = individ1.dna[i][j][c]
                        child_dna2[i][j][c] = individ2.dna[i][j][c]
                    else:
                        child_dna1[i][j][c] = individ2.dna[i][j][c]
                        child_dna2[i][j][c] = individ1.dna[i][j][c]

        return [self.createIndivid(child_dna1), self.createIndivid(child_dna2)]
        
    def crossover1(self, individ1, individ2):
        child_dna1 = copy.deepcopy(individ1.dna)
        child_dna2 = copy.deepcopy(individ2.dna)
        for i in range(self.fmmModel.nlayer):
            for j in range(self.fmmModel.nx):
                w = random.random()
                if w < 0.5:
                    child_dna1[i][j] = individ1.dna[i][j]
                    child_dna2[i][j] = individ2.dna[i][j]
                else:
                    child_dna1[i][j] = individ2.dna[i][j]
                    child_dna2[i][j] = individ1.dna[i][j]

        return [self.createIndivid(child_dna1), self.createIndivid(child_dna2)]
        
    def _crossover(self, individ1, individ2):
        return self.crossover4(individ1, individ2)

    def drawError (self, individ1, individ2, images_path):
        dna_m1 = self.getModel_FD(individ1.dna)   
        dna_m2 = self.getModel_FD(individ2.dna)   
        abs_error_m = copy.deepcopy (dna_m1)
        rel_error_m = copy.deepcopy (dna_m1)
        
        
        rel_error = 0
        abs_error = 0
        for i in range(dna_m1.nx):
            for j in range (dna_m1.nz):                
                abs_error_m.v[i][j] = dna_m2.v[i][j] - dna_m1.v[i][j]
                rel_error_m.v[i][j] = (abs_error_m.v[i][j]) / dna_m1.v[i][j]*100
                rel_error += abs(rel_error_m.v[i][j])
                abs_error += abs(abs_error_m.v[i][j])
                
        if images_path != None:
            abs_error_m.draw ('', images_path +'_abs_error.png', min_ = -500, max_=500, cmap = 'seismic')
            rel_error_m.draw ('', images_path +'_rel_error.png', min_ = -20, max_=20, cmap = 'seismic')
            
        abs_error = abs_error/(dna_m1.nx*dna_m1.nz)
        rel_error = rel_error/(dna_m1.nx*dna_m1.nz)

        return abs_error, rel_error

    def getGatherIndices (self, shot):
        return self._gatherModelIndex [shot]

    def getGatherCacheKey (self, shot, dna):
        gather_indices = self.getGatherIndices (shot)
        submodel = []
        for j in gather_indices:
            a = []
            for i in range(self.fmmModel.nlayer):
                a.append(dna[i][j])
            submodel.append(a)
#        print ("shot", shot, "submodel", submodel )
        return str(shot) + str(submodel)
        
        
    def applyLayerConstr(self, individ):
        for i in range(self.fmmModel.nlayer):
            for j in range(self.fmmModel.nx):
                individ.dna[i][j][0] = max(individ.dna[i][j][0], self._thickConstr[i][0])
                individ.dna[i][j][0] = min(individ.dna[i][j][0], self._thickConstr[i][1])
                
                individ.dna[i][j][1] = max(individ.dna[i][j][1], self._velConstr[i][0])
                individ.dna[i][j][1] = min(individ.dna[i][j][1], self._velConstr[i][1])

        

#def GA_test (helper, dna, pop_size, mutation = 0.1):
#    correct = helper.fitness(dna)
#    correct = helper.applyConstraints(correct)
#    print ('Correct answer:', correct)
#    
#    for i in range(pop_size):
#        dna = helper.mutate(dna, mutation)
#        
##        print(dna)
#        fit = helper.fitness(dna)
##        print (dna, fit)
##        print(i)
#        if fit > correct:
#            print ('We have a problem:', i, dna, fit)
#            exit()
        
 
#def MonteCarlo (helper, correct_dna, pop_size, mutation = 0.1):
#    correct = helper.fitness(correct_dna)
#    print ('Correct answer:', correct)
#    
#    best_dna = helper.random_individ()
#    
#    best_fit = helper.fitness(best_dna)
#    
#    for i in range(pop_size):
#        
#        if mutation > 0:
#            new_dna = helper.mutate(best_dna, mutation)
#        else:
#            new_dna = helper.random_individ()
#            
#        new_fit = helper.fitness(new_dna)
##        print (dna, fit)
#        if new_fit > best_fit:
#            best_dna = new_dna
#            best_fit = new_fit
#            print ('New best:', i, best_dna, best_fit)
#            
#
# Main driver
# Generate a population and simulate GENERATIONS generations.
#
def create_childs (helper, mutation, population):
    # Selection
    individ1 = weighted_choice(population)
    individ2 = weighted_choice(population)

    # Crossover
    childs = helper.crossover(individ1, individ2)
    
    for i in range(len(childs)):
        # Mutate and add back into the population.
        childs[i] = helper.mutate(childs[i], mutation)
        
#        childs[i] = helper.fitness(childs[i])
#        print ("childs[i]", childs[i].dna, childs[i].fitness)
    return childs

def create_new_population (helper, mutation, population):
               
    pop_size = len (population)            
    new_population = []

    # TEST child count
    individ1 = weighted_choice(population)
    individ2 = weighted_choice(population)
    childs = helper.crossover(individ1, individ2)
    childs_count = len (childs)
    
    if g_sc != None:
        par_pop = g_sc.parallelize (range(pop_size/childs_count))
            
        new_population_ = par_pop.map(lambda i: create_childs (helper.g_broadcastHelper.value, mutation, population)).collect()
        new_population = []
        for childs in new_population_:
            for child in childs:
                new_population.append(child)
        return new_population
    else:
        # Select two random individuals, based on their fitness probabilites, cross
        # their genes over at a random point, mutate them, and add them back to the
        # population for the next iteration.
        while len (new_population) < pop_size:
            childs = create_childs (helper, mutation, population)
            for child in childs:
                new_population.append(child)
    #            print ("after crosover")
    
        return new_population
    

def create_new_population_polygam (helper, mutation, population):
    population = helper.crossoverPolygam(population)
    
    for i in range(len(population)):
        population[i] = helper.mutate(population[i], mutation)
        
    return population;

#def smooth (f, N) :
#    return numpy.convolve (f, numpy.ones ((N,))/N, mode='valid')

def get_improvement (func, fitness_smooth_len):
    y = func[-fitness_smooth_len:]
    x = [i for i in range(len(y))]
    
    from scipy import stats
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    return slope/intercept
    

def PS_run_on_population (helper, correct_individ, images_path, population, 
        generatoin_count = 30, fitness_smooth_len = 1000):  
    
    population = helper.weight (population)

    best_generation = 0
    global_best_individ = helper.getBest(population)
    
    print ('Init')
#    helper.print_weight (population)
    print ("global best individual", global_best_individ.fitness)
    
    # print start point
    helper.draw (global_best_individ, images_path + "start")
    helper.drawError (correct_individ, global_best_individ, images_path + 'start')
        

    convergence_best_func = []
#    convergence_best_func.append (0)
#    convergence_best_func.append (global_best_fitness)

    convergence_aver_func = []
#    convergence_aver_func.append (0)
    
    rel_error_func = []
    abs_error_func = []

    omega = 0.5
    cp = 0.5
    cg = 0.5
    
    # Simulate all of the generations.
    for generation in range(generatoin_count):
        for p in population:            
            if generation > 0 and p.fitness_PS == 0: 
                throw (42)
            
            if p.fitness > p.fitness_PS:
                p.fitness_PS = p.fitness
                p.best = p.dna
    
            if p.fitness > global_best_individ.fitness_PS:
                global_best_individ = p
                best_generation = generation                        
                
        new_pop = []


        for p in population:
            rp = random.random()
            rg = random.random()
            
            for i in range(helper.fmmModel.nlayer):
                for j in range(helper.fmmModel.nx):
                    for c in range(2):
                        # check that clone copy v also
                        p.v[i][j][c] = omega*p.v[i][j][c]  \
                                        + cp * rp * (p.best[i][j][c] - p.dna[i][j][c]) \
                                        + cg * rg * (global_best_individ.dna[i][j][c] - p.dna[i][j][c])
  
            for i in range(helper.fmmModel.nlayer):
                for j in range(helper.fmmModel.nx):
                    for c in range(2):
                        p.dna[i][j][c] += p.v[i][j][c]


#            print ("dna", p.dna)
            
#            exit (0)
#            v = p.v + c1 * r1 * (p.best - p.dna) \
#                    + c2 * r2 * (global_best_individ.dna - p.dna)
#            p.dna = p.dna + v

            helper.applyConstraints(p)
            helper.applyLayerConstr(p)
            
#            p = helper.fitness(p)

            new_pop.append(p)
        
        population = helper.weight (new_pop)
        
#        new_best = helper.getBest (population)
#        if new_best.fitness > global_best_individ.fitness:
        print ("new global best!")
        helper.draw (global_best_individ, images_path + 'iter_' + str(generation))
        helper.drawError (correct_individ, global_best_individ, images_path + 'iter_' + str(generation))
    
#            print ("after draw")
        writeArray (images_path + 'global_best', global_best_individ.dna)

    
        print ('Generation', generation)
#        helper.print_weight (population)

        print ("global best individual", best_generation, global_best_individ.fitness)
            
        weight_aver = sum((item.fitness for item in population))/len(population)
        print ("Total weight", weight_aver)

        convergence_best_func.append (global_best_individ.fitness)
        convergence_aver_func.append (weight_aver)
        
        abs_error, rel_error = helper.drawError (correct_individ, global_best_individ, None)
    
        rel_error_func.append(rel_error)
        print ("rel error", rel_error)
        
        abs_error_func.append(abs_error)
        print ("abs error", abs_error)
         
#        convergence_aver_func_smooth = smooth (convergence_aver_func, fitness_smooth_len)
#        convergence_best_func_smooth = smooth (convergence_best_func, fitness_smooth_len)            
        
        
#        if generation % 10 == 0:
        import model_FD
        model_FD.draw_convergence ([convergence_best_func], "Generation", "Fitness", "Best fitness", images_path + 'convergence_best.png')
        model_FD.draw_convergence ([convergence_aver_func], "Generation", "Fitness", "Aver fitness", images_path + 'convergence_aver.png')
        model_FD.draw_convergence ([rel_error_func], "Generation", "Error", "Relative Error of best", images_path + 'rel_error.png')
        model_FD.draw_convergence ([abs_error_func], "Generation", "Error", "Absolute Error of best", images_path + 'abs_error.png')
#           
#            model_FD.draw_convergence ([convergence_aver_func_smooth], "Generation", "Fitness", "Aver fitness", images_path + 'convergence_aver_smooth.png')
#            model_FD.draw_convergence ([convergence_best_func_smooth], "Generation", "Fitness", "Best fitness", images_path + 'convergence_best_smooth.png')

        improuvement = get_improvement (convergence_best_func, fitness_smooth_len)
        print ("improuvement", improuvement)
 
        if generation > fitness_smooth_len and improuvement < 0.001 or generation == generatoin_count-1:                
            final_dna = helper.maximize (global_best_individ.dna)      
            writeArray (images_path + "final", final_dna)           
            final_individ = helper.fitness(helper.createIndivid (final_dna))
            helper.draw_gathers = True
            helper.draw (final_individ, images_path + "final")
            helper.drawError (correct_individ, global_best_individ, images_path + 'final')
            print ("final best individual", final_individ.fitness)        
            return population
#        
    return population
#        helper.fitness(local_best_ind,  image_path=images_path+'gen_'+str(generation))


def GA_run_on_population (helper, correct_individ, images_path, population, 
        generatoin_count = 30, mutation = 0.1, fitness_smooth_len = 1000):  
    
    population = helper.weight (population)

    best_generation = 0
    global_best_individ = helper.getBest(population)
    
    print ('Init')
#    helper.print_weight (population)
    print ("global best individual", global_best_individ.fitness)
    
    # print start point
    helper.draw (global_best_individ, images_path + "start")
    helper.drawError (correct_individ, global_best_individ, images_path + 'start')
        

    convergence_best_func = []
#    convergence_best_func.append (0)
#    convergence_best_func.append (global_best_fitness)

    convergence_aver_func = []
#    convergence_aver_func.append (0)
    
    rel_error_func = []
    abs_error_func = []

#    convergence_aver_func_smooth = []
#    convergence_aver_func_smooth.append (0)
    
    # Simulate all of the generations.
    for generation in range(generatoin_count):
    
        population = create_new_population (helper, mutation, population)
        population = helper.weight (population)    
        
#        print ("population", [individ.fitness for individ in population])

        
        local_best_individ = helper.getBest(population)
        
        if local_best_individ.fitness > global_best_individ.fitness:
            global_best_individ = local_best_individ
            best_generation = generation
            print ("new global best!")
            helper.draw (global_best_individ, images_path + 'iter_' + str(generation))
            helper.drawError (correct_individ, global_best_individ, images_path + 'iter_' + str(generation))
        
#            print ("after draw")
            writeArray (images_path + 'global_best', global_best_individ.dna)

    
        print ('Generation', generation)
#        helper.print_weight (population)

        print ("global best individual", best_generation, global_best_individ.fitness)
        print ("local best individual", local_best_individ.fitness)
            
        weight_aver = sum((item.fitness for item in population))/len(population)
        print ("Total weight", weight_aver)

        convergence_best_func.append (local_best_individ.fitness)
        convergence_aver_func.append (weight_aver)
        
        abs_error, rel_error = helper.drawError (correct_individ, local_best_individ, None)
    
        rel_error_func.append(rel_error)
        print ("rel error", rel_error)
        
        abs_error_func.append(abs_error)
        print ("abs error", abs_error)
         
        population_file = images_path + 'poputalion'
        writeArray (population_file, population)

        
#        if generation % 10 == 0:
        import model_FD
        model_FD.draw_convergence ([convergence_best_func], "Generation", "Fitness", "Best fitness", images_path + 'convergence_best.png')
        model_FD.draw_convergence ([convergence_aver_func], "Generation", "Fitness", "Aver fitness", images_path + 'convergence_aver.png')
        model_FD.draw_convergence ([rel_error_func], "Generation", "Error", "Relative Error of best", images_path + 'rel_error.png')
        model_FD.draw_convergence ([abs_error_func], "Generation", "Error", "Absolute Error of best", images_path + 'abs_error.png')
#           
#            model_FD.draw_convergence ([convergence_aver_func_smooth], "Generation", "Fitness", "Aver fitness", images_path + 'convergence_aver_smooth.png')
#            model_FD.draw_convergence ([convergence_best_func_smooth], "Generation", "Fitness", "Best fitness", images_path + 'convergence_best_smooth.png')

        improuvement = get_improvement (convergence_best_func, fitness_smooth_len)
        print ("improuvement", improuvement)
 
        if generation > fitness_smooth_len and improuvement < 0.00001 or generation == generatoin_count-1:                
            final_dna = helper.maximize (global_best_individ.dna)      
            writeArray (images_path + "final", final_dna)           
            final_individ = helper.fitness(helper.createIndivid (final_dna))
            helper.draw_gathers = True
            helper.draw (final_individ, images_path + "final")
            helper.drawError (correct_individ, global_best_individ, images_path + 'final')
            print ("final best individual", final_individ.fitness)        
            return population
#        
    return population
#        helper.fitness(local_best_ind,  image_path=images_path+'gen_'+str(generation))

def testError (helper, correct_dna, error, images_path = None):
    correct_individ = helper.fitness(helper.createIndivid (correct_dna))
    print ('before error', correct_individ.fitness)
    if images_path != None:
        helper.draw (correct_individ, images_path + "correct")
    
    dna = copy.deepcopy(correct_dna)
    
    for i in range(len(dna)):
        for j in range(len(dna[0])):
            for c in range(len(dna[0][0])):
                dna [i][j][c] += dna [i][j][c]*random.uniform(-error, error)
        
    testLocalOptimization (helper, dna, images_path)

def testLocalOptimization (helper, dna, images_path = None):
    individ = helper.fitness(helper.createIndivid (dna))
    print ('after before local', individ.fitness)

    if images_path != None:
        helper.draw (individ, images_path + "base")
       
    new_dna = helper.maximize_nelder_mead (dna)   
    individ = helper.fitness(helper.createIndivid (new_dna))
    print ('after maximize_nelder_mead', individ.fitness)
    if images_path != None:
        helper.draw (individ, images_path + "nelder_mead")
    
    new_dna = helper.maximize_basinhopping (dna)   
    individ = helper.fitness(helper.createIndivid (new_dna))
    print ('after maximize_basinhopping', individ.fitness)
    if images_path != None:
        helper.draw (individ, images_path + "basinhopping")
        
    new_dna = helper.maximize_BFGS (dna)   
    individ = helper.fitness(helper.createIndivid (new_dna))
    print ('after maximize_BFGS', individ.fitness)
    if images_path != None:
        helper.draw (individ, images_path + "BFGS")
        
    new_dna = helper.maximize_CMA (dna)   
    individ = helper.fitness(helper.createIndivid (new_dna))
    print ('after maximize_CMA', individ.fitness)
    if images_path != None:
        helper.draw (individ, images_path + "CMA")
     
        
def testErrors (helper, correct_dna, images_path):
    for error in xrange(10,15,1):
        testError (helper, correct_dna, error/100., images_path + 'error' + str(error) + '_')
                
def GA_run (helper, images_path, correct_dna,
        pop_size = 20, generatoin_count = 30, mutation = 0.01):  
    helper.print_info()
    
    if g_sc != None:
        helper_spark = copy.deepcopy (helper)
        helper.g_broadcastHelper = g_sc.broadcast(helper_spark)

    if not os.path.exists(images_path):
        os.makedirs(images_path)
        
    correct_individ = helper.fitness(helper.createIndivid(correct_dna))
    print ('Correct answer:', correct_dna.fitness)
    
    helper.draw (correct_individ, images_path + "correct")
    
    # Generate initial population. This will create a list of POP_SIZE strings,
    # each initialized to a sequence of random characters.
    population = helper.random_population(pop_size)
    
    population = GA_run_on_population(helper, correct_individ, images_path, population, generatoin_count, mutation)
    return population