import random
import copy
import numpy
import os
import subprocess
import math

#
def generate1DModel (nx, nz, dx, dz, interfaces):
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
        
    
def calcEnergy_FD (c, vel, image_path = None):
    import model_FD
    
    assert (c.back == 1)
    assert (c.snap == -1)

    test_name = c.path + 'tests/' 
    if not os.path.exists(test_name):
        os.makedirs(test_name)
            
    c.vp_file = test_name + 'vel.bin'
    c.wfl_file = test_name + 'snaps.bin'
    param_name = test_name + 'param_bw.txt'
    
    # write model
    vel.writeToFile (c.vp_file)

    c.writeToFile(param_name)   
    
    # command line args
    args = ['/home/cloudera/TRM/acoustic_FD_TRM/bin',
            param_name]
    
    FNULL = open(os.devnull, 'w')
    subprocess.call(args, stdout=FNULL, stderr=subprocess.STDOUT)
        
    # read snaps
    snaps = model_FD.snapShot(c.wfl_file, vel)
    # only 1 snaps should be red
    assert(len(snaps.spans) == 1)
    
    # take last snap
    final = snaps.spans[0]
    
    # take source position
    [sx, sz] = c.getSourcePosition()
#    print ('sx', sx, 'sz', sz)
    
    
    [six, siz] = final.getIndex(sx, sz)
    six = int (six)
    siz = int (siz)

    total_power = 0
    for i in range (final.nx):
        for j in range (final.nz):
            total_power += final.v[i][j]**2

    #print ('total_power', total_power)
    
    entropy = 0
    for i in range (final.nx):
        for j in range (final.nz):
            p_norm = final.v[i][j]**2/total_power
            entropy += p_norm*math.log(p_norm)
            
    entropy = -entropy
    
    energy = final.v[six][siz]**2
    
#    for i in range (mask.nx):
#        for j in range (mask.nz):
#            [x,z] = final.getCoordByIndex (i, j)
#            dist = math.sqrt((sx- x)**2 + (sz-z)**2)
#            # dist 0 - taper = 1
#            # dist = area - taper = 0.5
#            taper = 1/(dist/area + 1)
#            mask.v[i][j] = taper**4

#    mask.draw ('', 'mask.png')
    
    #take energy in source point
#    energy = 0
#    for i in range (mask.nx):
#        for j in range (mask.nz):
#            e = final.v[i][j]**2
#            taper = mask.v[i][j]
#            energy += e*taper
#    
    # just for test
    if image_path != None:
        final.draw ('test', figure_name = image_path + '.png', cmap = 'gray', norm = 1e-7)        
        vel.draw ('', image_path +'_model.png', min_ = 1500, max_=5000)
        c.drawEnargyAtSource(image_path + '_energy_at_source.png')
        
        zoom = final.zoom(sx, sz, 250, 250)
        zoom.draw ('zoom', figure_name = image_path + '_zoom.png', cmap = 'gray', norm = 1e-7)        
    
            
    # detele scratch files
#    import shutil    
#    shutil.rmtree(test_name)
    
    return energy, entropy
        
def calcMisfitEnergy_FD (c, vel, g, image_path = None):
    assert (c.back == 0)
    assert (c.snap == 0)

    test_name = c.path + 'tests/'
    if not os.path.exists(test_name):
        os.makedirs(test_name)
            
    c.vp_file = test_name + 'vel.bin'
    c.wfl_file = test_name + 'snaps.bin'
    c.gather_file = test_name + 'gather.bin'
    param_name = test_name + 'param_bw.txt'
    
    # write model
    vel.writeToFile (c.vp_file)

    # TODO optimize me
    c.writeToFile(param_name)   
    
    # command line args
    args = ['/home/cloudera/TRM/acoustic_FD_TRM/bin',
            param_name]
    
    FNULL = open(os.devnull, 'w')
    subprocess.call(args, stdout=FNULL, stderr=subprocess.STDOUT)
        
    calc_g = c.readGather ()
    
    if image_path != None:
        calc_g.draw ('test', figure_name = image_path + '_calc.png', cmap = 'gray')        
  
    norm_energy = 0
    energy = 0
    for i in range (g.ntr):
        for j in range (g.nt):
            norm_energy += g.v[j][i] **2
            calc_g.v[j][i] -= g.v[j][i] 
            energy += calc_g.v[j][i]**2
    
    # just for test
    if image_path != None:
        calc_g.draw ('test', figure_name = image_path + '_misfit.png', cmap = 'gray')        
        
        vel.draw ('', image_path +'_model_FD.png')

            
    # detele scratch files
#    import shutil    
#    shutil.rmtree(test_name)
    
    return energy/norm_energy

def calc_energy_RT (g, tt):
    energy = 0
    for i in range (g.ntr):
        j = int (tt[i]/g.dt)
        if j < 0 or j >= g.nt:
            continue
        amp = g.v[j][i]
        energy += amp ** 2
        
    return energy
    
def calc_semb_RT (g, tt):
    sum_ = 0
    sum_sq = 0
    for i in range (g.ntr):
        j = int (tt[i]/g.dt)
        if j < 0 or j >= g.nt:
            continue
        amp = g.v[j][i]
        sum_ += amp
        sum_sq += amp ** 2
        
    if sum_sq == 0:
        return 0
        
    return (sum_**2)/(sum_sq*g.ntr)
    
def generateGeom_RT (xs, ys, nr, dx, dy):
    import geom_RT
    
    g = geom_RT.Geom()
    s = geom_RT.GeomSource(xs, ys, 0)
    ind = int(nr/2)
    for i in range (-ind,ind + 1):
        x = s.x + i*dx
        y = s.y + i*dy
        s.addRec(geom_RT.GeomRec(x, y, s.z, 0))
        
    g.addSource(s)
    return g 
    
def calcTT_RT (m, geom_name, figure_name = None):
    import geom_RT
    
    tmp_dir = 'tmp'
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
        
    rays_name = tmp_dir + '/rays.txt'
    tt_name = tmp_dir + '/tt.txt'
#    geom_name = tmp_dir + '/geom.txt'
    model_name = tmp_dir + '/model.txt'
    
    m.writeToFile (model_name)
#    g.writeToFile (geom_name)
       
    # command line args
    args = [
            '/media/sf_Google_Drive/Phd/Near surface/FAT/sources/tomo3D/tt_forward3d/tt_forward3d',
            '-g', # bending
            '-G'+geom_name,
            '-M'+model_name,
            '-R'+rays_name,
#            '-N4/4/4/0/8/1e-4/1e-7'
            ]
            
#    print(' '.join(args))

    with open(tt_name, 'w') as f:
        subprocess.call(args, stdout=f)
    
    g_tt = geom_RT.Geom ()
    g_tt.readFromFile (tt_name)
    
    # draw
    if figure_name != None:
        # dbg
        import forward_model
        import rays
        rs = rays.Rays ()
        rs.readFromFile(rays_name)
        
        forward_model.drawX (m, rs, m.ly()/2, 'Velocity and Ray Path', figure_name, 'jet', draw_points=True)
    
    tt = []
    s = g_tt.sources[0]
    for r in s.rec:
        tt.append(r.t)
    
    return tt
    

def weighted_choice(items):
  """
  Chooses a random element from items, where items is a list of tuples in
  the form (item, weight). weight determines the probability of choosing its
  respective item. Note: this function is borrowed from ActiveState Recipes.
  """
  weight_total = sum((item[1] for item in items))
  
#  l = []
#  for item, weight in items:
#    l.append(weight/weight_total*100)
#  print ('distrib', l)
  
  n = random.uniform(0, weight_total)
  for item, weight, fitness in items:
    if n < weight:
      return item
    n = n - weight
  return item

def bound_random(v, dv, v1, v2):
    rv = 100000
    if v == None:
        v = int((v1 + v2)/2)
    if dv == None:
        dv = int((v1 + v2)/2)
    while v + rv not in range(v1, v2) : 
        rv = random.randrange(-dv, dv, 1)
        
    return v+rv
  
def round_v (v1):
    dv1 = 10
    return round(v1/dv1)*dv1

def round_z (v1):
    dv1 = 5
    return round(v1/dv1)*dv1
    
def random_v1(v1 = None):
    return round_v(bound_random (v1, None, 1500, 5000))
  
def random_v2(v2 = None):
    return round_v(bound_random (v2, None, 1500, 5000))
    
def random_z1(z1 = None):
    return round_z(bound_random (z1, None, 50, 250)) # WARNINNG min depth 50m
    
def random_v(v1 = None):
    return round_v(bound_random (v1, None, 1500, 5000))    
    
class GA_helper ():   
    
    def init(self, c, g, m):
        self.c = c
        self.g = g
        self.m = m
                
        source_x = m.lx()/2
        source_y = m.ly()/2

        self.geom = generateGeom_RT(source_x, source_y, g.ntr, g.dh/1000., 0)
        self.geom_name = 'geom.txt'
        self.geom.writeToFile(self.geom_name) 
        
        self.fd_rt = -1   # FD or RT 
        self.max_min = -1 # maximization or minimization
        
#        self.rt_energy_semb = 1 # semblance or energy
#        self.fd_energy_entropy = 0 # entropy or energy or misfit_energy
        
    def define_FD_energy(self):
        self.fd_rt = 0
        self.fd_energy_entropy = 0

    def define_FD_entropy(self):
        self.fd_rt = 0
        self.fd_energy_entropy = 1       

    def define_FD_misfit_enerfy (self):
        self.fd_rt = 0
        self.fd_energy_entropy = 2     
        
    def define_RT_energy(self):
        self.fd_rt = 1
        self.rt_energy_semb = 0

    def define_RT_semb(self):
        self.fd_rt = 1
        self.rt_energy_semb = 1          
        
    def print_info (self):
        if self.fd_rt == 0:
            print ('Finite difference')
            if self.fd_energy_entropy == 0:
                print ('Energy')
                self.max_min = 0
            if self.fd_energy_entropy == 1:
                print ('Entropy')
                self.max_min = 1
            if self.fd_energy_entropy == 2:
                print ('Misfit energy')
                self.max_min = 1
                
        if self.fd_rt == 1:
            print ('Raytracing')
            if self.rt_energy_semb == 0:
                self.max_min = 0
                print ('Energy')
            if self.rt_energy_semb == 1:
                print ('Entropy')
                self.max_min = 0
                
        if self.max_min == 0:
            print ('Maximization')
        if self.max_min == 1:
            print ('Minimization')
            

    def random_dna(self):
        print ('random_dna not implemented')
        
    def random_population(self, pop_size):
        pop = []
        for i in range(pop_size):       
            dna = self.random_dna()
            pop.append(dna)
        return pop
    
    def getModel_RT(self, dna):
        print ('getModel_RT not implemented')
        
    def getTT_RT(self, dna, figure_name = None):
        dna_m = self.getModel_RT(dna)
        if dna_m == None:
            return None
        tt = calcTT_RT(dna_m, self.geom_name, figure_name)   
        return tt        

    def fitness(self, dna, image_path=None):
        if self.fd_rt == 0:
            return self.__fitness_FD(dna, image_path)
        if self.fd_rt == 1:
            return self.__fitness_RT(dna, image_path)
                
    def __fitness_RT(self, dna, image_path=None):
        tt = self.getTT_RT(dna);
        if self.rt_energy_semb == 0:
            return calc_energy_RT (self.g, tt)
        if self.rt_energy_semb == 1:
            return calc_semb_RT (self.g, tt)

    def __fitness_FD(self, dna, image_path=None):
        dna_m = self.getModel_FD(dna)
        if self.fd_energy_entropy == 0 or self.fd_energy_entropy == 1:
            energy, entropy = calcEnergy_FD (self.c, dna_m, image_path=image_path)
            if self.fd_energy_entropy == 0:
                return energy
            if self.fd_energy_entropy == 1:
                return entropy
        if self.fd_energy_entropy == 2:
            return calcMisfitEnergy_FD (self.c, dna_m, self.g, image_path=image_path)
            
    def getModel(self, dna):
        if self.fd_rt == 0:
            return self.getModel_FD(dna)
        if self.fd_rt == 0:
            return self.getModel_RT(dna)
            
    def getModel_FD(self, dna):
        print ('getModel not implemented')
    
    def mutate(self, dna, mutation_chance):        
        print ('mutate not implemented')
    
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
            fitness_val = self.fitness(individual)

            # Generate the (individual,fitness) pair, taking in account whether or
            # not we will accidently divide by zero.
            pair = (individual, fitness_val, fitness_val)

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
            fitness_val = self.fitness(individual)
            
            
            # not we will accidently divide by zero.
            if fitness_val == 0:
                pair = (individual, max_weight, fitness_val)
            else:
                assert (1.0/fitness_val<max_weight)
                pair = (individual, 1.0/fitness_val, fitness_val)
 
            weighted_population.append(pair)
        return weighted_population
    
        
    def getBest(self, weighted_population):
        (best_ind, maximum_weight, best_fitness) = weighted_population[0]
                
        for i in range(len(weighted_population)):
            if weighted_population[i][1] >= maximum_weight:
                (best_ind, maximum_weight, best_fitness) = weighted_population[i]
                
        return (best_ind, maximum_weight, best_fitness)
        
    def print_weight (self, weighted_population) :      
        for i in range(len(weighted_population)):
            print (weighted_population[i])
            
        (best_ind, maximum_weight, best_fitness) = self.getBest(weighted_population)
        
        weight_total = sum((item[1] for item in weighted_population))

        print ("Best individual", best_ind, maximum_weight, best_fitness)
        print ("Total fitness", weight_total)

    def draw (self, individual, images_path):
        #print ('individual',individual)
#        dna_m = self.getModel_FD(individual)
        self.fitness (individual, image_path=images_path)
        
        tt = self.getTT_RT(individual)
        if tt!= None:
            self.g.draw (tt = tt, figure_name = images_path +'_gather.png')
  
#     
    def crossover(self, dna1, dna2):
        return self.crossover2(dna1, dna2)

#    def crossover1(self, dna1, dna2):
#        child1 = []
#        child2 = []
#        for c in range(len(dna1)):
#            w = random.random()
#            v1 = None
#            v2 = None
#            if c == 0:
#                v1 = round_z(w*dna1[c] + (1-w)*dna2[c])
#                v2 = round_z((1-w)*dna1[c] + w*dna2[c])
#            if c == 0:
#                v1 = round_v(w*dna1[c] + (1-w)*dna2[c])
#                v2 = round_v((1-w)*dna1[c] + w*dna2[c])         
#            if c == 0:
#                v1 = round_v(w*dna1[c] + (1-w)*dna2[c])
#                v2 = round_v((1-w)*dna1[c] + w*dna2[c])
#            child1.append(v1)
#            child2.append(v2)
#            
##        print ('cross parents', dna1, dna2)
##        print ('cross childs', child1, child2)
#        return child1, child2
#        
    def crossover2(self, dna1, dna2):
        child1 = []
        child2 = []
        for c in range(len(dna1)):
            w = random.random()
            if w < 0.5:
                child1.append(dna1[c])
                child2.append(dna2[c])
            else:
                child1.append(dna2[c])
                child2.append(dna1[c])
            
#        print ('cross parents', dna1, dna2)
#        print ('cross childs', child1, child2)
        return child1, child2
           
class GA_helperI1 (GA_helper):   

    def __init__(self, c, g, m):
        self.init(c,g,m)    

    def random_dna(self):
        z1 = random_z1()
        v1 = random_v1()
        v2 = random_v2()
        
        dna = [z1, v1, v2]
        return dna
        
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
    
    @staticmethod
    def fillModel_RT (m, z1, v1, v2):
        for i in range (m.nx()):
            for j in range (m.ny()):
                for k in range (m.nz()):
                    z = k*m.dz 
                    v = v1
                    if (z > z1):
                        v = v2
                        
                    m.v[i][j][k] = v
        return m
    
    def getModel_RT(self, dna):
        m = copy.deepcopy(self.m)
        dna_m = self.fillModel_RT(m, dna[0]/1000., dna[1]/1000., dna[2]/1000.)
        return dna_m
    
    def getModel_FD(self, dna):
        import model_FD
        
        m = model_FD.model(self.c.nx, self.c.nz, self.c.dh, self.c.dh)
        dna_m = self.fillModel1(m, dna[0], dna[1], dna[2])
        return dna_m
            
    def mutate(self, dna, mutation_chance):        
        for c in range(len(dna)):
            if random.random() <= mutation_chance:
                if c == 0:
#                    prev = dna[c]
                    dna[c] = random_z1(dna[c])
#                    print ('mutate z1', prev, dna[c])
                if c == 1:
#                    prev = dna[c]
                    dna[c] = random_v1(dna[c])
#                    print ('mutate v1', prev, dna[c])
                if c == 2:
#                    prev = dna[c]
                    dna[c] = random_v2(dna[c])
#                    print ('mutate v2', prev, dna[c])
        
        return dna
    
class GA_helperI2 (GA_helper):   
    
    def __init__(self, c, g, m):
        self.init(c, g, m)  

    def empty_model (self):
        import model_FD
        lx = self.c.nx*self.c.dh
        lz = self.c.nz*self.c.dh
        new_dx = 10001
        new_dz = 25
        
        rand_m = model_FD.model(int(lx/new_dx)+1, int(lz/new_dz)+1, new_dx, new_dz)
        return rand_m
        
    def random_dna(self):
        rand_m = self.empty_model()
        for i in range (rand_m.nx):
            for j in range (rand_m.nz):
                rand_m.v[i][j] = self.random_v()
                    
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
            
    def mutate(self, dna, mutation_chance):        
        for c in range(len(dna)):
            if random.random() <= mutation_chance:
                dna[c] = random_v(dna[c])
                
        return dna


class GA_helperI3 (GA_helper):   
    
    def __init__(self, c, g, m):
        self.init(c, g, m)  
        self.layer_count = 3
        
    def random_z(self, z1 = None):
        lz = int(self.c.nz*self.c.dh)
        if z1 == None:
            z1 = int(lz/2)
        return round_z(bound_random (z1, None, 0, lz))    
        
    def sort_by_z(self, dna):
        for i in range(len(dna)):
            for j in range(i+1, len(dna)):
                if dna[i][0] > dna[j][0]:
                    tmp = dna[i]
                    dna[i] = dna[j]
                    dna[j] = tmp
        return dna
            
    def random_dna(self):
        dna = [[0, random_v()]]
    
        for l in range(self.layer_count-1):
            dna.append([self.random_z(),random_v()])
        return self.sort_by_z(dna)
    
    def getModel_FD(self, dna):
#        dna = self.sort_by_z(dna)
#        print (dna)
        dna_m = generate1DModel(self.c.nx, self.c.nz, self.c.dh, self.c.dh, dna)
        return dna_m
            
    def mutate(self, dna, mutation_chance):        
        for c in range(len(dna)):   
             if random.random() <= mutation_chance:
                 if c != 0: # WARNING special case
                    dna[c][0] = self.random_z(dna[c][0])
             if random.random() <= mutation_chance:
                 dna[c][1] = random_v(dna[c][1])                
               
        dna = self.sort_by_z(dna)
#        print ('mutate dna', dna)
        return dna
        
    def crossover(self, dna1, dna2):
        dna1, dna2 = self.crossover2(dna1, dna2)
        dna1 = self.sort_by_z(dna1)
        dna2 = self.sort_by_z(dna2)
#        print ('crossover dna1', dna1)
#        print ('crossover dna2', dna2)
        
        return dna1, dna2

#
# Main driver
# Generate a population and simulate GENERATIONS generations.
#

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
    weighted_population = helper.weight (population)

    best_generation = 0
    (global_best_ind, global_maximum_weight, global_best_fitness) = helper.getBest(weighted_population)
    
    print ('Init')
    helper.print_weight (weighted_population)
    print ("global best individual", global_best_ind, global_maximum_weight, global_best_fitness)

    # Simulate all of the generations.
    for generation in range(generatoin_count):
    
        population = []

        # Select two random individuals, based on their fitness probabilites, cross
        # their genes over at a random point, mutate them, and add them back to the
        # population for the next iteration.
        for _ in range(pop_size/2):
            # Selection
            ind1 = weighted_choice(weighted_population)
            ind2 = weighted_choice(weighted_population)

            # Crossover
            ind1, ind2 = helper.crossover(ind1, ind2)

            # Mutate and add back into the population.
            population.append(helper.mutate(ind1, mutation))
            population.append(helper.mutate(ind2, mutation))
            
        weighted_population = helper.weight (population)    
        (local_best_ind, local_maximum_weight, best_fitness) = helper.getBest(weighted_population)
        
        if local_maximum_weight > global_maximum_weight:
            global_best_ind = local_best_ind
            global_maximum_weight = local_maximum_weight
            global_best_fitness = best_fitness
            best_generation = generation
            print ("new global best!")
            helper.draw (local_best_ind, images_path + str(generation))

    
        print ('Generation', generation)
#        helper.print_weight (weighted_population)

        print ("global best individual", best_generation, global_best_ind, global_maximum_weight, global_best_fitness)
        print ("local best individual", local_best_ind, local_maximum_weight, best_fitness)
            
        weight_total = sum((item[1] for item in weighted_population))
        print ("Total weight", weight_total)

#        helper.fitness(local_best_ind,  image_path=images_path+'gen_'+str(generation))