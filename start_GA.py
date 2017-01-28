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


def modelingOneModel(model_path, model):
    import model_FD
    
    c = model_FD.config(model_path)
    c.absPath()
    
    c.nz = 50
    c.nu0 = 20 #Hz
    c.snap = -1
    
        
    # enable PML
#    c.free_surf = 0
#    c.npml=3
#    c.snap = c.nt

    
    vel = GA.generate1DModel (c.nx, c.nz, c.dh, c.dh, model)
    
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
#    g.muteDirect (g, c.dr, -0.1, 2000)
#    g.muteDirect (g, c.dr, 0.16, 3500)
#    g.muteOffset (g, c.dr, 0, 1000)
#    g.norm(1e+3)
    
    
#    g.phase()
#    

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

#    
#def modelingOneModelEvgeny(model_path):
#        
#    c = model_FD.config(model_path)    
#    param_name = c.path + 'param_fw.txt'
#    c.readFromFile(param_name)   
##    c.absPath()
#
#
##    vel.writeToFile (c.vp_file)
##
#    m = model_FD.model (c.nx, c.nz, c.dh, c.dh)
#    m.readFromFile(c.vp_file)  
##    print (m.v)
##    exit ()
#    rho = copy.deepcopy(m)
#    rho.resetValues (1)
#    # use vel as roh
#    rho.writeToFile (c.rho_file)
#    #    
#    
##    c.generateGeomFiles()
#    
#    c.wfl_file = param_name = c.path + 'snap_fw.bin'
#    
#    param_name = c.path + 'param_fw.txt'
#    c.writeToFile(param_name)
#    
#    # command line args
#    args = ['/home/cloudera/TRM/acoustic_FD_TRM/bin',
#            param_name]
#    
##    log_name = c.path + 'log.txt'
##    with open(log_name, 'w') as f:
##        subprocess.call(args, stdout=f)
##    c.createImages ()
#    
##    # MUTE    
#    g = c.readGather ()
##    refl + refr
#    g.muteDirect (0.60, 9000, hyp=False, up = True)
#    g.muteDirect (0.80, 2000, hyp=True)
##    g.norm(1e+3)
#    
#    g.muteOffset (0, 1000)
#
##    refr
#    g.muteDirect (0.65, 9000, hyp=False, up = True)
#    g.muteDirect (0.60, 2000, hyp=True, up = False)
##    g.norm(1e+2)
#    g.draw ('', c.path + 'forward_gather_mute.png')
#    exit ()
#
#
#    new_nather_name = c.gather_file + '_mute'
#    g.writeToFile (new_nather_name)    
#    c.gather_file = new_nather_name
#
#    
#    param_name = c.path + 'param_bw.txt'
#    c.readFromFile(param_name)
##    c.absPath()
#    c.gather_file = new_nather_name
#    c.writeToFile(param_name)
#    
#    # command line args
#    args = ['/home/cloudera/TRM/acoustic_FD_TRM/bin',
#            param_name]
#    
#    log_name = c.path + 'log.txt'
#    with open(log_name, 'w') as f:
#        subprocess.call(args, stdout=f)
#        
#    c.createImages ()
#
#def modelingOneModelEvgenyFullTRM(model_path):
#        
#    c = model_FD.config(model_path)    
#    param_name = c.path + 'param_fw.txt'
#    c.readFromFile(param_name)   
#    c.absPath()
#
#
##    vel.writeToFile (c.vp_file)
##
#    m = model_FD.model (c.nx, c.nz, c.dh, c.dh)
#    m.readFromFile(c.vp_file)  
##    print (m.v)
##    exit ()
#    rho = copy.deepcopy(m)
#    rho.resetValues (1)
#    # use vel as roh
#    rho.writeToFile (c.rho_file)
#    #    
#    
#    c.generateGeomFilesRoundModel()
#    
#    c.wfl_file = param_name = c.path + 'snap_fw.bin'
#    
#    param_name = c.path + 'param_fw_full_TRM.txt'
#    c.writeToFile(param_name)
#    
#    # command line args
#    args = ['/home/cloudera/TRM/acoustic_FD_TRM/bin',
#            param_name]
#    
#    log_name = c.path + 'log.txt'
#    with open(log_name, 'w') as f:
#        subprocess.call(args, stdout=f)
#    c.createImages ()
#    
#
#    g = c.readGather ()
#    g.norm(1e+1)
#
#    new_nather_name = c.gather_file + '_mute'
#    g.writeToFile (new_nather_name)  
#    
##    param_name = c.path + 'param_bw.txt'
##    c.readFromFile(param_name)
##    c.absPath()
#    c.back = 1
#    c.gather_file = new_nather_name
#    param_name = c.path + 'param_bw_full_TRM.txt'        
#    c.writeToFile(param_name)
#    
#    # command line args
#    args = ['/home/cloudera/TRM/acoustic_FD_TRM/bin',
#            param_name]
#    
#    log_name = c.path + 'log.txt'
#    with open(log_name, 'w') as f:
#        subprocess.call(args, stdout=f)
#        
#    c.createImages ()

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
    
    
def testObjective (helper, correct_dna, figure_name=None):
    correct_fitness = helper.fitness(correct_dna)
    print ('correct', correct_dna, correct_fitness)
    
    lz = 20
    dz = 5
    nz = lz/dz
    
    lv1 = 100
    dv1 = 25
    nv1 = lv1/dv1

    lv2 = 50
    dv2 = 25
    nv2 = lv2/dv2
    
    import numpy
    cube = numpy.zeros((2*nz+1,2*nv2+1,2*nv1+1))
    for z in range(-nz,nz+1):
        for v1 in range(-nv1,nv1+1):
            for v2 in range(-nv2,nv2+1):
                dna = [correct_dna[0]+z*dz, correct_dna[1]+v1*dv1, correct_dna[2]+v2*dv2]
                fitness = [0]
                fitness = helper.fitness(dna)
                cube[z+nz][v2+nv2][v1+nv1]=fitness[0]/correct_fitness[0]/2
                if fitness > correct_fitness:
                    print ('wrong', dna, fitness)
                    
#    cube[nz*2][v2+nv2][0]=1
    
    z = numpy.arange(correct_dna[0]-nz*dz,correct_dna[0]+nz*dz + dz,dz)
    print (z)
    v1 = numpy.arange(correct_dna[1]-nv1*dv1,correct_dna[1]+nv1*dv1 + dv1,dv1)
    print (v1)
    v2 = numpy.arange(correct_dna[2]-nv2*dv2,correct_dna[2]+nv2*dv2 + dv2,dv2)
    print (v2)
    GA.plotcube (cube,v1,v2,z,
              x_label = 'V1',
              y_label = 'V2',
              z_label = 'Z',
              figure_name = figure_name)
    

def prepare_gather(c, images_path, ):
   
#    
    g = c.readGather ()
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
    g.muteDirect (0.2, 3500, hyp = False)
###    g = model_FD.muteOffset (g, c.dr, 0, 1000)

    # mute reflection
#    g.muteDirect (0.05, 2000, hyp = True)
#    g.muteDirect (0.15, 3500, hyp = False)
#    g.draw ('', images_path + 'forward_gather.png')
#    exit()
    return g    

if __name__ == "__main__":
    model_path = '//home/cloudera/TRM/acoustic_FD_TRM/tests/orhan/'
    if not os.path.exists(model_path):
        os.makedirs(model_path)
        
#    modelingOneModelEvgenyFullTRM(model_path)
#    exit ()

#    modelingOneModel(model_path, [[0, 2000], [125, 3500], [200, 4000]])
    
    import model_FD    
    c = model_FD.config(model_path)    
    param_name = c.path + 'param_bw.txt'
    c.readFromFile(param_name)   

    import model_RT
    model_name = 'model.txt'
    m = model_RT.VelModel() 
    m.emptyModel(c.nx*c.dh/1000., 0.05, c.nz*c.dh/1000., 0.025, 0.025, 0.01)
    m.writeToFile(model_name)
    
    images_path = c.path + 'GA_images_evgeny_phase_semb/' 
    if not os.path.exists(images_path):
        os.makedirs(images_path)
     
    g = prepare_gather(c, images_path)
    
        
    helper = GA.GA_helperI1 (c, g, m)
    correct_dna = [125., 2000., 3500.]
#    helper = GA.GA_helperI3 (c, g, m)
#    correct_dna = [[0., 2000.], [125., 3500.], [200., 4000.]]

    helper.define_RT_semb()

    
    testObjective (helper, correct_dna, figure_name = images_path + 'testObjective.png')
    exit ()
    
    GA.GA_run (helper, images_path, correct_dna,
        pop_size = 20, generatoin_count = 100, mutation = 0.1)    