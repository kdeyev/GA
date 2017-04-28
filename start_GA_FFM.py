#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 00:06:47 2017

@author: cloudera
"""
import os
import GA_FMM as GA
import copy
import subprocess

def modelingMultiGatherModel(model_path, vel):
    import model_FD
    
    c = model_FD.config(model_path)
    c.absPath()
    
    c.nz = 50
    c.nu0 = 20 #Hz
#    c.snap = 20
    c.snap = -1
    c.g_ns = 11
    
        
    # enable PML
#    c.free_surf = 0
#    c.npml=3
#    c.snap = c.nt

    
#    vel = GA.generate1DModel (c.nx, c.nz, c.dh, c.dh, model)
    
    vel.draw ('', model_path + 'vel_gen.png')    
#    exit ()
    vel.writeToFile (c.vp_file)

    rho = copy.deepcopy(vel)
    rho.resetValues (1)
    # use vel as roh
    rho.writeToFile (c.rho_file)
    
    gathers = []
    for shot in range(c.g_ns):
        c.generateGeomFiles(shot)
        
        c.wfl_file = param_name = c.path + 'snap_fw.bin'
        
        param_name = c.path + 'param_fw.txt'
        c.writeToFile(param_name)
        
        # command line args
        args = ['/home/cloudera/TRM/acoustic_FD_TRM/bin',
                param_name]
        
        log_name = c.path + 'log.txt'
        with open(log_name, 'w') as f:
            subprocess.call(args, stdout=f)

        g = c.readGather (shot)
        
        gathers.append(g)
    
#        c.back = 1
#        # enable PML
#    #    c.free_surf = 0
#    #    c.npml=3
#        
#        c.wfl_file = param_name = c.path + 'snap_bw.bin'
#    
#    #    # MUTE    
#        g = c.readGather (shot)
#    #    g.muteDirect (-0.1, 2000)
#    #    g.muteDirect (0.135, 3500, hyp = False)
#    #    g.muteOffset (0, 1000)
#    #    g.norm(1e+3)
#        
#        new_nather_name = c.gather_file + '_mute'
#        g.writeToFile (new_nather_name)    
#        c.gather_file = new_nather_name
#    
#        
    param_name = c.path + 'param_bw.txt'
    c.writeToFile(param_name)

def modelingOneModel(model_path, vel):
    import model_FD
    
    c = model_FD.config(model_path)
    c.absPath()
    
    c.nz = 50
    c.nu0 = 20 #Hz
    c.snap = 20
    c.g_ns = 1
    
        
    # enable PML
#    c.free_surf = 0
#    c.npml=3
#    c.snap = c.nt

    
#    vel = GA.generate1DModel (c.nx, c.nz, c.dh, c.dh, model)
    
    vel.draw ('', model_path + 'vel_gen.png', min_=1500, max_=5000)    
#    exit ()
    vel.writeToFile (c.vp_file)

    
    rho = copy.deepcopy(vel)
    rho.resetValues (1)
    g = c.readGather (0)
    rho.writeToFile (c.rho_file)
    
    c.generateGeomFiles(0)
    
    c.wfl_file = param_name = c.path + 'snap_fw.bin'
    
    param_name = c.path + 'param_fw.txt'
    c.writeToFile(param_name)
    
    # command line args
    args = ['/home/cloudera/TRM/acoustic_FD_TRM/bin',
            param_name]
    
    log_name = c.path + 'log.txt'
    with open(log_name, 'w') as f:
        subprocess.call(args, stdout=f)
        
    c.createImages ()
    
    c.back = 1
    # enable PML
#    c.free_surf = 0
#    c.npml=3
    
    c.wfl_file = param_name = c.path + 'snap_bw.bin'

#    # MUTE    
    g = c.readGather (0)
#    g.muteDirect (-0.1, 2000)
#    g.muteDirect (0.135, 3500, hyp = False)
#    g.muteOffset (0, 1000)
#    g.norm(1e+3)
    
    
    new_nather_name = c._gather_file + '_mute_0'
    g.writeToFile (new_nather_name)    
    c._gather_file = new_nather_name  
    c.g_gather_file = c.g_gather_file + '_0' + '_mute'
#    print (new_nather_name)
#    exit()

#    rho = GA.sourceWall_FD(g.sPos(), rho)
#    rho.writeToFile (c.rho_file)
#    rho.draw ('', model_path + 'rho_gen_bw.png') 
    
    vel.draw ('', model_path + 'vel_gen_bw.png', min_=1500, max_=5000)    
    vel.writeToFile (c.vp_file)

    
    param_name = c.path + 'param_bw.txt'
    c.writeToFile(param_name)
    
    # command line args
    args = ['/home/cloudera/TRM/acoustic_FD_TRM/bin',
            param_name]
    
    log_name = c.path + 'log.txt'
    with open(log_name, 'w') as f:
        subprocess.call(args, stdout=f)
        
    c.createImages ()

def run_SU (args, g):
    import subprocess
    import numpy
    import segypy
    
    p = subprocess.Popen(
        args, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
       
    v = numpy.transpose(g.v)
    in_data = bytearray()
    for i in range (g.ntr):
        header = bytearray(240)

        segypy.setDefaultSTHValue (g.nt, header, "ns")
        segypy.setDefaultSTHValue (int(g.dt*1000000), header, "dt")   
        in_data.extend(header)

        data = bytearray(g.nt*4) 
        for j in range (g.nt):
            segypy.setValue(v[i][j], data, j*4, 'f', segypy.endian, 1)
        in_data.extend(data)
        
    out, err = p.communicate(input=in_data)
    out_data_array = bytearray(out)
           
    trace_len = g.nt*4 + 240
    for i in range(g.ntr):
        v[i] = segypy.getValue(out_data_array, i*trace_len + 240, 'f', segypy.endian, g.nt)[0]
        
    g.v = numpy.transpose(v)
    
    return g
    
def testEntropy (helper, figure_name):
        
    correct_dna = [125., 2000., 3500.]
    fitness, info = helper.draw (correct_dna, images_path + 'correct')
    print ('correct', info)

    
    wrong_dna = [130., 2680., 3340.]   
    fitness, info = helper.draw (wrong_dna, images_path + 'wrong')
    print ('wrong', info)
    
    exit ()
 
#    
#def testObjective (helper, correct_dna, figure_name=None):
#    correct_fitness, correct_info = helper.fitness(correct_dna)
#    print ('correct', correct_dna, correct_fitness, correct_info)
#    
#    lz = 60
#    dz = 10
#    nz = lz/dz
#    
#    lv1 = 500
#    dv1 = 100
#    nv1 = lv1/dv1
#
#    lv2 = 500
#    dv2 = 100
#    nv2 = lv2/dv2
#    
#    import numpy
#    cube = numpy.zeros((2*nz+1,2*nv2+1,2*nv1+1))
#    
#    
#    import copy
#    info_cubes = {}
#    for k in correct_info.keys():
##        print (k)
#        info_cubes[k] = copy.deepcopy(cube)
#    
#    
#    for z in range(-nz,nz+1):
#        for v1 in range(-nv1,nv1+1):
#            for v2 in range(-nv2,nv2+1):
#                dna = [correct_dna[0]+z*dz, correct_dna[1]+v1*dv1, correct_dna[2]+v2*dv2]
##                fitness = [0]
#                individ = helper.fitness(Individ(dna))
#                cube[z+nz][v2+nv2][v1+nv1]=fitness/correct_fitness/2
#    
##                if fitness < correct_fitness:
##                    print ('wrong', dna, fitness, info)
#    
##                print(dna, info)
#                for k in correct_info.keys():
#                    correct_value = correct_info[k]
#                    value = info[k]
#                    c = info_cubes[k]
#                    c[z+nz][v2+nv2][v1+nv1]=value/correct_value/2
#
#    z = numpy.arange(correct_dna[0]-nz*dz,correct_dna[0]+nz*dz + dz,dz)
##    print (z)
#    v1 = numpy.arange(correct_dna[1]-nv1*dv1,correct_dna[1]+nv1*dv1 + dv1,dv1)
##    print (v1)
#    v2 = numpy.arange(correct_dna[2]-nv2*dv2,correct_dna[2]+nv2*dv2 + dv2,dv2)
##    print (v2)
#    GA.plotcube (cube,v1,v2,z,
#              x_label = 'V1',
#              y_label = 'V2',
#              z_label = 'Z',
#              figure_name = figure_name+'cube.png')
#    
#    for k in info_cubes.keys():
#        c = info_cubes[k]
#        GA.plotcube (c,v1,v2,z,
#              x_label = 'V1',
#              y_label = 'V2',
#              z_label = 'Z',
#              figure_name = figure_name + k + '.png')


def prepare_gather_mute_direct_offset(c, images_path ):
    gathers = []
    for shot in range(c.g_ns):
        g = c.readGather (shot)
        g.norm(1e+5)
    #    g.norm_ampl = 1e-03
        
        g.muteDirect (-0.1, 2000)
        g.muteDirect (0.135, 3500, hyp = False)
        g.muteOffset (0, 1000)
        g.draw ('', images_path + 'forward_gather.png')
        gathers.append(g)    
    return gathers      
        
def prepare_gather_phase(c, images_path):
    gathers = []
    for shot in range(c.g_ns):
        g = c.readGather (shot)
        g.norm_ampl = None # autonorm
    #    g.draw ('', images_path + 'orig.png')
        
        g.muteDirect (0.135, 3500, hyp = False)
        
        g.norm(1e+4)    
        g= run_SU(['/home/cloudera/cwp/bin/suaddnoise', 'sn=10000'], g)
        g= run_SU(['/home/cloudera/cwp/bin/suattributes', 'mode=phase'], g)        
        gathers.append(g)    
    return gathers      

def prepare_gather_agc(c, images_path):
    gathers = []
    for shot in range(c.g_ns):
        g = c.readGather (shot)
        
        # mute
        g.muteDirect (0.075, 5000, hyp = False)
        
        g.norm(1e+4)    
        g= run_SU(['/home/cloudera/cwp/bin/sugain', 'agc=1', 'wagc=0.05'], g)

        g.norm_ampl = 1 # autonorm    
        if images_path != None:
            g.draw ('', images_path + 'orig_' + str(shot) + '.png')
        gathers.append(g)
    
    return gathers   
    
def write_gathers (c, gathers):
#    c._gather_file = c._gather_file + '_proc'
    for shot in range(c.g_ns):
        g = gathers[shot]
        new_nather_name = c.g_gather_file + '_proc_' + str(shot)
        g.writeToFile (new_nather_name)    
        c._gather_file = new_nather_name  
        
def read_gathers (c):
    gathers = []
    c.g_gather_file = c.g_gather_file + '_proc'
        
    for shot in range(c.g_ns):
        g = c.readGather (shot)
        gathers.append(g)
        
    return gathers
    
def mute_offset(orig_gathers, offset, images_path):
    gathers = []
    for shot in range(len(orig_gathers)):
        g = copy.deepcopy(orig_gathers[shot])
        # mute
        g.muteOffset (offset, 5000000)
        
        if images_path != None:
            g.draw ('', images_path + 'mute_offset_' + str(shot) + '.png')
        gathers.append(g)
    
    return gathers   


def prepare_gather(c, images_path): 
    gathers = []
    for shot in range(c.g_ns):
        g = c.readGather (shot)
    
    #    # fd
    #    g.norm(1e+7)
    #    g.norm_ampl = 1e-03
    
    #    #rt
        # agc
        g.norm_ampl = None # autonorm
    #    g.draw ('', images_path + 'orig.png')
        
        g.norm(1e+4)    
        g= run_SU(['/home/cloudera/cwp/bin/suaddnoise', 'sn=10000'], g)
        g= run_SU(['/home/cloudera/cwp/bin/suattributes', 'mode=phase'], g)
    #    g= run_SU(['/home/cloudera/cwp/bin/sugain', 'agc=1', 'wagc=0.05'], g)
    
        # mute
    ###    g.muteDirect (-0.1, 2000)
        g.muteDirect (0.135, 3500, hyp = False)
    ###    g = model_FD.muteOffset (g, c.dr, 0, 1000)
    
        # mute reflection
    #    g.muteDirect (0.05, 2000, hyp = True)
    #    g.muteDirect (0.135, 3500, hyp = False)
    #    g.draw ('', images_path + 'forward_gather.png')
    #    exit()
        gathers.append(g)
    return gathers    
        
        
def prepare_helper (modelGeom, gathers, correct_dna, images_path):
    
    vel_constr = [[1500, 2500],
                  [2500, 3500],
                  [4500, 5500]]
    thick_constr = [[0, 100],
                   [50, 150],
                   [0, 200]]
                   
               
    for i in range(len(correct_dna)):
        for j in range(len(correct_dna[0])):
            th = correct_dna[i][j][0]
            v = correct_dna[i][j][1]
            if th < thick_constr[i][0] or th > thick_constr[i][1]:
                print ("wrong constriant")
                exit (0)
            if v < vel_constr[i][0] or v > vel_constr[i][1]:
                print ("wrong constriant")
                exit (0)
                
    helper = GA.GA_helperI4 (modelGeom, gathers, 0.01, True, len(correct_dna), len(correct_dna[0]), vel_constr, thick_constr)
    helper.define_FMM_energy_semb()
    
#    helper.putSpikeToGathers (correct_dna);
      
#    helper.addConstraint(GA.GA_helperI4_Constraint_Well (correct_dna))
    helper.addConstraint(GA.GA_helperI4_Constraint_V (correct_dna))
#    helper.addConstraint(GA.GA_helperI4_Constraint_Point(correct_dna, True, True))
        
    individ = helper.fitness(helper.createIndivid(correct_dna))
    print ('Correct answer:', individ.fitness)
 
    if images_path != None:
        helper.draw (individ, images_path + "correct")
    return helper
            
def test_correct_dna (modelGeom, orig_gathers, correct_dna,
        pop_size = 20, generatoin_count = 30, mutation = 0.1, images_path = None): 

    if not os.path.exists(images_path):
        os.makedirs(images_path)
        
    helper = prepare_helper(modelGeom, orig_gathers, correct_dna, None)
    print (helper)
        
    # Generate initial population. This will create a list of POP_SIZE strings,
    # each initialized to a sequence of random characters.
    population = helper.random_population(pop_size)

    for j in range(len(correct_dna[0])):
        for i in range(len(correct_dna)):
            for c in range(len(correct_dna[i][j])):
                population[j][i][j][c] = correct_dna[i][j][c]      
    
    return GA.GA_run_on_population(helper, images_path, population, generatoin_count, mutation)

        
if __name__ == "__main__":
#    model_path = '//home/cloudera/TRM/acoustic_FD_TRM/tests/evgeny/'
    model_path = 'C:/GA/tests/evgeny/'
#    model_path = "C:/Users/kostyad/Google Drive/Phd/Near surface/EAT/TRM/acoustic_FD_TRM/tests/evgeny/"
    if not os.path.exists(model_path):
        os.makedirs(model_path)
    
    import model_FD    
    c = model_FD.config(model_path)    
    param_name = c.path + 'param_bw.txt'
    c.readFromFile(param_name)   
    
    modelGeom  = model_FD.modelGeom (c.nx, c.nz, c.dh, c.dh)
    
    images_path = c.path + 'GA_images_evgeny_FMM/' 
    if not os.path.exists(images_path):
        os.makedirs(images_path)
     
#    gathers = prepare_gather_agc(c, None)
#    write_gathers (c, gathers)
    gathers = read_gathers (c)

    #    exit (0)
#    gathers = prepare_gather_mute_direct_offset(c, images_path)
    
#    new_nather_name = c.gather_file + '_mute'
#    g.writeToFile (new_nather_name)    
#    c.gather_file = new_nather_name
        
#    helper = GA.GA_helperI1 (c, g, m)
#    correct_dna = [125., 2000., 3500.]

#    helper = GA.GA_helperI3 (c, g, m)
#    correct_dna = [[0., 2000.], [125., 3500.], [200., 4000.]]


    correct_dna = [[[25., 2000.], [25., 2000.], [12.5, 2000.], [25., 2500.], [25., 2500.], [25., 2000.]],
                   [[100., 3000.], [100., 3000.], [120., 3500.], [100., 3500.], [120., 3500.], [100., 3500.]],
                   [[0., 5000.], [0., 5000.], [0., 5000.], [0., 5000.], [0., 5000.], [0., 5000.]]
                   ]

                   
#    test_correct_dna (c, m, gathers, correct_dna, 
#                      pop_size = 100, generatoin_count = 50, mutation = 0.1, images_path = images_path)
                   
    

               
    helper = prepare_helper (modelGeom, gathers, correct_dna, images_path)  

                   
#    modelingMultiGatherModel(model_path, helper.getModel_FD(correct_dna))
#    exit ()

#    testEntropy (helper, images_path)
#    exit ()
        
#    testObjective (helper, correct_dna, figure_name = images_path)
#    exit ()

#    GA.MonteCarlo(helper,correct_dna, 100000, mutation = 0)
#    exit()
    
#    GA.GA_test(helper,correct_dna, 1000000, mutation = 0.1)
#    exit()
        
    GA.GA_run (helper, images_path, correct_dna,
        pop_size = 1000, generatoin_count = 1000, mutation = 0.03)    