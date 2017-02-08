#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 00:06:47 2017

@author: cloudera
"""
import os
import GA
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
    c.g_ns = 3
    
        
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
#    c.snap = 20
    c.snap = -1
    
        
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
    
    c.generateGeomFiles()
    
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
    g = c.readGather ()
#    g.muteDirect (-0.1, 2000)
#    g.muteDirect (0.135, 3500, hyp = False)
#    g.muteOffset (0, 1000)
#    g.norm(1e+3)
    
    new_nather_name = c.gather_file + '_mute'
    g.writeToFile (new_nather_name)    
    c.gather_file = new_nather_name

    
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
 
    
def testObjective (helper, correct_dna, figure_name=None):
    correct_fitness, correct_info = helper.fitness(correct_dna)
    print ('correct', correct_dna, correct_fitness, correct_info)
    
    lz = 60
    dz = 10
    nz = lz/dz
    
    lv1 = 500
    dv1 = 100
    nv1 = lv1/dv1

    lv2 = 500
    dv2 = 100
    nv2 = lv2/dv2
    
    import numpy
    cube = numpy.zeros((2*nz+1,2*nv2+1,2*nv1+1))
    
    
    import copy
    info_cubes = {}
    for k in correct_info.keys():
#        print (k)
        info_cubes[k] = copy.deepcopy(cube)
    
    
    for z in range(-nz,nz+1):
        for v1 in range(-nv1,nv1+1):
            for v2 in range(-nv2,nv2+1):
                dna = [correct_dna[0]+z*dz, correct_dna[1]+v1*dv1, correct_dna[2]+v2*dv2]
#                fitness = [0]
                fitness, info = helper.fitness(dna)
                cube[z+nz][v2+nv2][v1+nv1]=fitness/correct_fitness/2
    
#                if fitness < correct_fitness:
#                    print ('wrong', dna, fitness, info)
    
#                print(dna, info)
                for k in correct_info.keys():
                    correct_value = correct_info[k]
                    value = info[k]
                    c = info_cubes[k]
                    c[z+nz][v2+nv2][v1+nv1]=value/correct_value/2

    z = numpy.arange(correct_dna[0]-nz*dz,correct_dna[0]+nz*dz + dz,dz)
#    print (z)
    v1 = numpy.arange(correct_dna[1]-nv1*dv1,correct_dna[1]+nv1*dv1 + dv1,dv1)
#    print (v1)
    v2 = numpy.arange(correct_dna[2]-nv2*dv2,correct_dna[2]+nv2*dv2 + dv2,dv2)
#    print (v2)
    GA.plotcube (cube,v1,v2,z,
              x_label = 'V1',
              y_label = 'V2',
              z_label = 'Z',
              figure_name = figure_name+'cube.png')
    
    for k in info_cubes.keys():
        c = info_cubes[k]
        GA.plotcube (c,v1,v2,z,
              x_label = 'V1',
              y_label = 'V2',
              z_label = 'Z',
              figure_name = figure_name + k + '.png')


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
        g.norm_ampl = None # autonorm
    #    g.draw ('', images_path + 'orig.png')
        
        # mute
        g.muteDirect (0.135, 3500, hyp = False)
        
        g.norm(1e+4)    
        g= run_SU(['/home/cloudera/cwp/bin/sugain', 'agc=1', 'wagc=0.05'], g)
        
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

    
class GA_helperI4_Constraint (GA.GA_constraint):
    def __init__(self, correct_dna):
        self.correct_dna = correct_dna
        self.nlayer = len(correct_dna)
        self.nx = len(correct_dna[0])
        
    def applyConstraint (self, dna):
        for i in range(self.nlayer):
            for j in range(self.nx):
                if i == 0: # first layer
                    # use correct velocity
                    dna[i][j][1] = correct_dna[i][j][1]
        
                
        return dna
        
    
if __name__ == "__main__":
    model_path = '//home/cloudera/TRM/acoustic_FD_TRM/tests/evgeny/'
    if not os.path.exists(model_path):
        os.makedirs(model_path)

#    helper = GA.GA_helperI4 (c, g, m, 2, 3)
#    correct_dna = [[[100., 2000.], [125., 2500.], [100., 2000.]],
#                   [[50., 3500.], [25., 3000.], [50., 3500.]]
#                   ]
#    modelingOneModel(model_path, helper.getModel_FD(correct_dna))
#    exit ()
    
    import model_FD    
    c = model_FD.config(model_path)    
    param_name = c.path + 'param_bw.txt'
    c.readFromFile(param_name)   

    import model_RT
    model_name = 'model.txt'
    m = model_RT.VelModel() 
    m.emptyModel(c.lx()/1000., 0.05, c.lz()/1000., 0.025, 0.025, 0.01)
    m.writeToFile(model_name)
    
    images_path = c.path + 'GA_images_evgeny_FMM/' 
    if not os.path.exists(images_path):
        os.makedirs(images_path)
     
    gathers = prepare_gather_agc(c, images_path)
#    gathers = prepare_gather_mute_direct_offset(c, images_path)
    
#    new_nather_name = c.gather_file + '_mute'
#    g.writeToFile (new_nather_name)    
#    c.gather_file = new_nather_name
        
#    helper = GA.GA_helperI1 (c, g, m)
#    correct_dna = [125., 2000., 3500.]

#    helper = GA.GA_helperI3 (c, g, m)
#    correct_dna = [[0., 2000.], [125., 3500.], [200., 4000.]]


    correct_dna = [[[50., 2000.], [50., 2000.], [25., 2000.], [50., 2500.], [50., 2500.], [50., 2000.]],
                   [[50., 3500.], [50., 3500.], [75., 3500.], [50., 3500.], [50., 3500.], [50., 3500.]]
                   ]
    helper = GA.GA_helperI4 (c, gathers, m, 0.01, True, len(correct_dna), len(correct_dna[0]))

                   
#    modelingMultiGatherModel(model_path, helper.getModel_FD(correct_dna))
#    exit ()

    helper.define_FMM_semb()
    
    constraint = GA_helperI4_Constraint (correct_dna)
    helper.addConstraint(constraint)

#    testEntropy (helper, images_path)
#    exit ()
    

    
#    testObjective (helper, correct_dna, figure_name = images_path)
#    exit ()

#    GA.MonteCarlo(helper,correct_dna, 100000, mutation = 0)
#    exit()
    
#    GA.GA_test(helper,correct_dna, 1000000, mutation = 0.1)
#    exit()
    
    GA.GA_run (helper, images_path, correct_dna,
        pop_size = 100, generatoin_count = 1000, mutation = 0.1)    