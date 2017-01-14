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


def modelingOneModel(model_path):
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

    
    vel = GA.generate1DModel (c.nx, c.nz, c.dh, c.dh, [[0, 2000],[125, 3500],[200, 3000]])
    
#    vel.draw ('', model_path + 'vel_gen.png')    
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
    g.norm(1e+4)
    
    
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



if __name__ == "__main__":
    model_path = '//home/cloudera/TRM/acoustic_FD_TRM/tests/orhan/'
    if not os.path.exists(model_path):
        os.makedirs(model_path)
        
#    modelingOneModelEvgenyFullTRM(model_path)
#    exit ()

#    modelingOneModel(model_path)

    import model_FD    
    c = model_FD.config(model_path)    
    param_name = c.path + 'param_bw.txt'
    c.readFromFile(param_name)   

    import model_RT
    model_name = 'model.txt'
    m = model_RT.VelModel() 
    m.emptyModel(c.nx*c.dh/1000., 0.05, c.nz*c.dh/1000., 0.025, 0.025, 0.01)
    m.writeToFile(model_name)
    
#    #    # MUTE    
    g = c.readGather ()
#    g = model_FD.muteDirect (g, c.dr, -0.1, 2000)
#    g = model_FD.muteDirect (g, c.dr, 0.16, 3500)
#    g = model_FD.muteOffset (g, c.dr, 0, 1000)
#    g.draw ('', model_path + 'forward_gather.png')
#    exit()
 
    # 3
#    helper = GA.GA_helperI3 (c, g, m)
#    correct_dna = [[0, 2000],[125, 3500],[200, 3000]]

    # 1
    helper = GA.GA_helperI1 (c, g, m)
    correct_dna = [125, 2000, 3500]

    helper.define_FD_entropy()
    images_path = c.path + 'GA_images/' 
    
    GA.GA_run (helper, images_path, correct_dna,
        pop_size = 20, generatoin_count = 30, mutation = 0.1)    