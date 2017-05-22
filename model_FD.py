#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 21:49:46 2016

@author: cloudera
"""

import numpy
import matplotlib.pyplot as plt
import copy
import os
import math

font_size = 60
global_font = {'fontname':'serif', 'size':str(font_size)}

class modelGeom ():
    def __init__ (self, nx, nz, dx, dz, sx=0, sz=0):
        self.initGeom (nx, nz, dx, dz, sx, sz)
                        
    def initGeom (self, nx, nz, dx, dz, sx=0, sz=0):
        self.nx = nx
        self.nz = nz
        self.dx = dx
        self.dz = dz
        self.sx = sx
        self.sz = sz
#        self.x_nodes = [i*self.dx + self.sx for i in range (self.nx)]
#        self.z_nodes = [i*self.dz + self.sz for i in range (self.nz)]
                        
                        
#        print ('nx', self.nx)
#        print ('nx', self.nz) 
#        print ('dx', self.dx)
#        print ('dz', self.dz)
#        print (self.x_nodes)
#        print (self.z_nodes)
#    

    def lx (self):
        return self.sx + (self.nx-1)*self.dx
        
    def lz (self):
        return self.sz + (self.nz-1)*self.dz

    def getIndex (self, x, z):
        i = (x - self.sx)/self.dx
        j = (z - self.sz)/self.dz
        return [i, j]

    def getCoordByIndex (self, i,j):
#        assert (i.is_in1teger())
#        assert (j.is_integer())
        x = self.sx + i*self.dx
        z = self.sz + j*self.dz
        return [x, z]

               
class model (modelGeom):
    def __init__ (self, nx, nz, dx, dz, sx=0, sz=0):
        self.initGeom (nx, nz, dx, dz, sx, sz)
        
        self.v = numpy.zeros ((nx,nz))
        self.x_nodes = [i*self.dx + self.sx for i in range (self.nx)]
        self.z_nodes = [i*self.dz + self.sz for i in range (self.nz)]
                        
                        
#        print ('nx', self.nx)
#        print ('nx', self.nz) 
#        print ('dx', self.dx)
#        print ('dz', self.dz)
#        print (self.x_nodes)
#        print (self.z_nodes)
#    
    
    def readFromFile (self, filename): 
        self.v = numpy.fromfile(filename, dtype=numpy.float32, count=self.nx*self.nz)
        self.v.flags.writeable = True
        self.v = numpy.reshape(self.v, (self.nx, self.nz))    
    
    def readFromBuf (self, buf): 
        self.v = numpy.frombuffer(buf, dtype=numpy.float32, count=self.nx*self.nz)
        self.v = numpy.reshape(self.v, (self.nz, self.nx))
        self.v = numpy.transpose(self.v)                    

    def writeToFile (self, filename): 
        array = numpy.reshape(self.v, self.nx*self.nz)   
        array.astype('float32').tofile(filename)

    def resetValues (self, v):        
        for i in range (self.nx):
            for j in range (self.nz):
                self.v[i][j] = v

    def zoom (self, center_x, center_z, lx, lz):      
        center_x = int(center_x/self.dx)*self.dx
        center_z = int(center_z/self.dz)*self.dz

        lx = int(lx/self.dx)*self.dx
        lz = int(lz/self.dz)*self.dz
        
        start_x = center_x - lx
        start_x = max(start_x, self.sx)
        start_z = center_z - lz
        start_z = max(start_z, self.sz)
        
        start_x_index = int(start_x/self.dx)
        start_z_index = int(start_z/self.dz)
#        start_x = start_x_index*self.dx
#        start_z = start_z_index*self.dz
        
        end_x = center_x + lx
        end_x = min(end_x, self.lx())
        end_z = center_z + lz
        end_z = min(end_z, self.lz)
        
        nx = int((end_x - start_x)/self.dx)
        nz = int((end_z - start_z)/self.dz)

        out = model (nx, nz, self.dx, self.dz, start_x, start_z)        
        
        start_x_index = int(start_x/self.dx)
        start_z_index = int(start_z/self.dz)

        for i in range (out.nx):
            for j in range (out.nz):
                out.v[i][j] = self.v[i+start_x_index][j+start_z_index]

        return out
            
    def draw(self, label, figure_name = None, show=False, cmap = 'jet', norm = None, min_ = None, max_= None):    
        plt.rcParams['figure.figsize'] = 60, 15
#        fig = plt.figure (figsize=(60, 15))
        fig = plt.figure()
#        plt.rc('font', family='serif', size=60)
        ax = plt.subplot(111)
        ax.set_title(label, **global_font)
        plt.ylabel('Depth (m)', **global_font)
        plt.xlabel('Location (m)', **global_font)
             
        label_x = copy.deepcopy(self.x_nodes)
        label_z = copy.deepcopy(self.z_nodes)
        for i in range(len(label_x)):
            if i%200 != 0:
                label_x[i] = ""
    
        for i in range(len(label_z)):
            if i%15 != 0:
                label_z[i] = ""
                
        
#        print ('draw nx', self.nx)
#        print ('draw nz', self.nz)
        
        plt.xticks(range(len(label_x)), label_x, **global_font)
        plt.yticks(range(len(label_z)), label_z, **global_font)
        
        v = numpy.transpose(self.v)
#        a = self.dx/self.dz
#        a = 4
        a = 'auto'
        vmin = min_
        vmax = max_
        if norm != None:
            vmin = -norm
            vmax = norm
            
        cax = ax.imshow(v, aspect=a, vmin=vmin, vmax=vmax, cmap=cmap)
        
    #    fig.colorbar(cax, orientation='horizontal')
        fig.colorbar(cax).ax.tick_params(labelsize=font_size)
        
    #    plt.grid(True)
        if figure_name != None:
            plt.savefig(figure_name,bbox_inches='tight', dpi=100)
#            plt.savefig(figure_name,bbox_inches='tight')
            
        if show:
            plt.show ()
        plt.gcf().clear()
        plt.close('all')
        
    def setValue(self, x, z, v):
        [i, j] = self.getIndex (x, z)
        self.setValueByIndex(i, j, v)

    def setValueByIndex(self, i, j, v):
#        assert (i.is_integer())
#        assert (j.is_integer())
#
        i = int (i)
        j = int (j)
        
        if i in range (self.nx):
            if j in range (self.nz):
#                print ("setValueByIndex", i, j)
                self.v.flags.writeable = True
                self.v[i][j] = v
        
    def getValue(self, x, z):
        [i, j] = self.getIndex (x, z)
        return self.getValueByIndex(i, j)
        
    def getInterp(self): 
        from scipy.interpolate import RegularGridInterpolator
        
        fn = RegularGridInterpolator((self.x_nodes,self.z_nodes), self.v)
        return fn
         
#        x0 = round (i)
#        z0 = round (j)
#        
#        x = x - x0
#        z = z - z0
#        
#        x1 = x0 + 1
#        z1 = z0 + 1
#
#        output = (
#        self.v[x0][z0]*(1-x)*(1-z) +
#        self.v[x1][z0]*x*(1-z) +
#        self.v[x0][z0]*(1-x)*(1-z) +
#        self.v[x0][z1]*(1-x)*(1-y)*z +
#        self.v[x1][z1]*x*z +
#        self.v[x0][z1]*(1-x)*z +
#        self.v[x1][z0]*x*(1-z) +
#        self.v[x1][z1]*x*z)
#
#        return output
#        
    def getValueByIndex(self, i, j):
#        assert (i.is_integer())
#        assert (j.is_integer())
#
        i = int (i)
        j = int (j)
        
        if i in range (self.nx):
            if j in range (self.nz):
                return self.v[i][j]

#        k = round (k)
#        
#        if i < 0 or i >= self.nx ():
#            return None
#            
#        if j < 0 or j >= self.ny ():
#            return None
#        
#        if k < 0 or k >= self.nz ():
#            return None
#            
#            
#        return self.v[i][j][k]
#        
#        x0 = round (i)
#        y0 = round (j)
#        z0 = round (k)
#        
#        x = x - x0
#        y = y - y0
#        z = z - z0
#        
#        x1 = x0 + 1
#        y1 = y0 + 1
#        z1 = z0 + 1
#
#        output = (
#        self.v[x0][y0][z0]*(1-x)*(1-y)*(1-z) +
#        self.v[x1][y0][z0]*x*(1-y)*(1-z) +
#        self.v[x0][y1][z0]*(1-x)*y*(1-z) +
#        self.v[x0][y0][z1]*(1-x)*(1-y)*z +
#        self.v[x1][y0][z1]*x*(1-y)*z +
#        self.v[x0][y1][z1]*(1-x)*y*z +
#        self.v[x1][y1][z0]*x*y*(1-z) +
#        self.v[x1][y1][z1]*x*y*z)
#
#        return output
#
#        return self.v[i][j][k]

#
#class gather ():
#    def __init__ (self, ntr, nt, dt):
##        print ('ntr', ntr)
#        self.ntr = ntr
#        self.nt = nt
#        self.dt = dt
#        self.v = numpy.zeros ((ntr,nt))
#        
#    def readValues (self, filename): 
#        self.v = numpy.fromfile(filename, dtype=numpy.float32, count=self.ntr*self.nt)
#        self.v = numpy.reshape(self.v, (self.nt, self.ntr))
##        self.v = numpy.transpose(self.v)    
#
#    def writeToFile (self, filename): 
#        array = numpy.reshape(self.v, self.ntr*self.nt)   
#        array.astype('float32').tofile(filename)
#        
#    def draw(self, label, figure_name = None, show=False, norm = 1e-6, cmap = 'gray'):    
#        plt.rcParams['figure.figsize'] = 40, 30
##        fig = plt.figure (figsize=(60, 15))
#        fig = plt.figure()
#        plt.rc('font', family='serif', size=60)
#        ax = plt.subplot(111)
#        ax.set_title(label)
#        plt.ylabel('Time (s)')
#        plt.xlabel('Trace')
#        
##        print ('ntr', len(self.v), len(self.v[0]))
#             
#        label_t = [i*self.dt for i in range(self.nt)]  
#        
#        for i in range(len(label_t)):
#            if i%50 != 0:
#                label_t[i] = ""
#                
#        plt.yticks(range(len(label_t)), label_t)
#    
#        a = 'auto'
#        if norm != None:
#            cax = ax.imshow(self.v, aspect=a, vmin=-norm, vmax=norm, cmap=cmap)
#        else:
#            cax = ax.imshow(self.v, aspect=a, cmap=cmap)
#        
#    #    fig.colorbar(cax, orientation='horizontal')
#        fig.colorbar(cax)
#        
#    #    plt.grid(True)
#        if figure_name != None:
#            plt.savefig(figure_name,bbox_inches='tight', dpi=100)
##            plt.savefig(figure_name,bbox_inches='tight')
#            
#        if show:
#            plt.show ()
#        plt.gcf().clear()

class gather ():
    def __init__ (self, nt, dt, dh, s_pos, rec_pos):
#        print ('ntr', ntr)
        self.ntr = len(rec_pos)
        self.nt = nt
        self.dt = dt
        self.dh = dh
        self.norm_ampl = None
        self.v = numpy.zeros ((nt,self.ntr))
        self.s_pos = s_pos
        self.rec_pos = rec_pos
        
    def sPos(self):
        return self.s_pos
        
    def rPos(self):
        return self.rec_pos
        
    def readValues (self, filename): 
        self.v = numpy.fromfile(filename, dtype=numpy.float32, count=self.ntr*self.nt)
        self.v = numpy.reshape(self.v, (self.nt, self.ntr))
#        self.v = numpy.transpose(self.v)    

    def writeToFile (self, filename): 
        array = numpy.reshape(self.v, self.ntr*self.nt)   
        array.astype('float32').tofile(filename)
        
    def draw(self, label = '', figure_name = None, show=False, norm = None, cmap = 'gray', tt = None):    
        import matplotlib.pyplot as plt
        
        plt.rcParams['figure.figsize'] = 40, 30
#        fig = plt.figure (figsize=(60, 15))
        fig = plt.figure()
#        plt.rc('font', family='serif', size=60)
        ax = plt.subplot(111)
        ax.set_title(label, **global_font)
        plt.ylabel('Time (s)', **global_font)
        plt.xlabel('Trace', **global_font)
        
        if norm == None:
            norm = self.norm_ampl
#        print ('ntr', len(self.v), len(self.v[0]))
             
        label_t = [(i*self.dt*1000)/1000. for i in range(self.nt)]  
        
        for i in range(len(label_t)):
            if i%50 != 0:
                label_t[i] = ""
                
        plt.yticks(range(len(label_t)), label_t, **global_font)
    
        a = 'auto'
        if norm != None:
            cax = ax.imshow(self.v, aspect=a, vmin=-norm, vmax=norm, cmap=cmap)
        else:
            cax = ax.imshow(self.v, aspect=a, cmap=cmap)
            
        if tt != None:
            x = [v for v in range(len(tt))]
            y = [v/self.dt for v in tt]
            ax.plot(x, y, 'r-', linewidth=10, markersize=10)
        
    #    fig.colorbar(cax, orientation='horizontal')
        fig.colorbar(cax).ax.tick_params(labelsize=font_size)

        
    #    plt.grid(True)
        if figure_name != None:
            plt.savefig(figure_name,bbox_inches='tight', dpi=100)
#            plt.savefig(figure_name,bbox_inches='tight')
            
        if show:
            plt.show ()
        plt.gcf().clear()
        plt.close('all')
            
    def muteDirect (self, t0, v1, hyp = True, up=False) : 
        taper_samp = 20
        for i in range (self.ntr):
            offset = abs(self.rec_pos[i][0] - self.s_pos[0])
            if hyp:
                tmax = math.sqrt((offset / v1)**2 + t0**2)
            else:
                tmax = offset / v1 + t0
            
            sampmax = int(tmax/self.dt)
            if not up:
                for j in range (sampmax, self.nt):
                    f = 0
                    if j-sampmax < taper_samp:
                        f = 1-(j-sampmax)/float(taper_samp)
        #                print (f)
                    self.v[j][i] *= f
            else:
                for j in range (0, sampmax):
                    f = 0
#                    if sampmax-j < taper_samp:
#                        f = 1-(j-sampmax)/float(taper_samp)
        #                print (f)
                    self.v[j][i] *= f
#        return g
        

    def muteOffset (self, o1, o2) : 
        for i in range (self.ntr):
            offset = abs(self.rec_pos[i][0] - self.s_pos[0])
            f = 1
            if offset >= o1 and offset <= o2:
                f = 0
            for j in range (self.nt):
                self.v[j][i] *= f
#        return g

    def norm(self, v):
        for i in range (self.ntr):
            for j in range (self.nt):
                self.v[j][i] *= v

    def phase(self):
        from scipy.signal import hilbert
        
        v = numpy.transpose(self.v)
        
        for i in range (self.ntr):
            analytic_signal = hilbert(v[i])
            phase = numpy.angle(analytic_signal)
            phase = numpy.unwrap(phase)
#            print  (phase)
            for j in range (self.nt):
                self.v[j][i] = phase[j]

class snapShot ():
    def __init__(self,filename, m):
        self.spans = []
        with open(filename,'rb') as f:
            while True:
                buf=f.read(m.nx*m.nz*4)
                if not buf: break
                m = copy.deepcopy(m)
                m.readFromBuf(buf)
                self.spans.append(m)
    
    def draw(self, path, norm = None):
#        or i in range(len(self.spans)):
#            self.spans[i].amin
#        m = min ()
#        norm = 0.000001
        
        for i in range(len(self.spans)):
            filename = path + 'snap' + str (i+1) + '.png'
            self.spans[i].draw('', figure_name = filename, cmap = 'gray', norm=norm)
        
        
#
#def generateModel (nx, nz, dx, dz, modeltype):
#    m = model(nx, nz, dx, dz) 
#   
##    m1 = int(modeltype)
#    
#    factor = 1
##    factor = None
##    if m1 == 0:
##        factor = 1
##    if m1 == 1:
##        factor = 0.97
##    if m1 == 2:
##        factor = 1.03
#           
##    print ('generateModel nx', m.nx)
##    print ('generateModel nz', m.nz)
#        
#    if modeltype == 2:
#        for i in range (m.nx):
#            for j in range (m.nz):
#            
#                z = j*dz 
#                
#                if z <= 20:
#                    v = 2300
#                elif z <= 40:
#                    v = 2600
#                elif z <= 60:
#                    v = 2900
#                elif z <= 80:
#                    v = 3200
#                elif z <= 100:
#                    v = 3300
#                else:
#                    v = 3500
#                    
#                m.v[i][j] = v*factor
#        return m    
#        
#        
#    z1 = None
#    z2 = None
#    v1 = None
#    v2 = None
#    dv = None
#
#    if modeltype == 1:
#        z1 = 100
#        z2 = 125   
#        v1 = 2000
#        v2 = 3500
#
#    if modeltype == 3:
#        z1 = 00
#        z2 = 125
#        v1 = 2000
#        v2 = 3500
#        dv = (v2-v1)/(z2-z1)
#
#    if modeltype == 4:
#        z1 = 100
#        z2 = 125
#        v1 = 2000
#        v2 = 3500
#        dv = (v2-v1)/(z2-z1)
#
#    if modeltype == 5:
#        z1 = 0
#        z2 = 250
#        v1 = 2000
#        v2 = 2200
#        dv = (v2-v1)/(z2-z1)
#
#
#    for i in range (m.nx):
#        for j in range (m.nz):
#            z = j*dz 
#            v = v1
#            if dv != None and z >= z1 and z <= z2:
#                v = v1 + dv*(z-z1)
#            if (z > z2):
#                v = v2
#                
#            m.v[i][j] = v*factor
#    
#    return m    

class config ():
    def __init__(self, path):
        self.path = path
        self.ox=0          # Grid origin in X and Z directions
        self.oz=0    
        
        self.nx=2001       # Number of grid points in X and Z directions
        self.nz=2001 
        self.dh=5          # Grid step (equal step in both directions should be used)
        
        self.vp_file= 'vel.bin'     # Files with velocity and density models
        self.rho_file='rho.bin'    # size=nz*nx, single precision


        self.nt=1000         # Number of time samples in seismogram
        self.dt=0.002        # Time sampling in seismogram 
        self.nu0=20
        self._gather_file = 'gather.bin'    # File to save the snapshots. size=nr*nt*n_snaps

        self._ns = 1
        self._nr=1001          # Number of geophones

        self._src_file= 'src.txt'     # Files that contain sources and geophones positions  
        self._rcv_file= 'rcv.txt'     # (see below their format)

        self.npml=50    # Size of PML layer in grid points
                        # (increase it to reduce artificial reflections, 
                        #  decrease it to reduce calculation time )
        self.damp=0.1   # Damping factor in PML
                        # (decrease it in case of instability)
        
        self.free_surf=0  # Use 1 for free-surface and 0 for PML at the top of the model 
        
        self.apprx=1    # Use 1 for 4th order scheme, 0 for 2nd order scheme.
                        # Note: model sampling should be >15 points per central wavelength for 4th order
                        #  and >30 points per central wavelength for 2nd order scheme
        
        self.snap=10    # Number of snapshots to save. Use 0 to skip saving the snapshots
        
        self.wfl_file = 'snap_fw.bin'    # File to save the snapshots. size=nr*nt*n_snaps

        self.verb=1          # Verbosity flag
        
        self.expl=0          # Exploiding reflector modeling. If 1 then half-velocity model will be used.
                        # Shots should be located inside the medium at the positions of reflectors.
                        # Note that half-velocity may lead to dispersion due to coarse grid step in 
                        # the model.
        
        self.back=0     # Use 0 for forward propagation, 1 for backward-propagation.
                        # In forward propagation mode the Ricker wavelet is used.
                        # In backward propagation mode the seismogram from file defined in "out="
                        # is propagated in reverse-time. In that case the zero-offset image
                        # is stored in the file defined by "wfl=" at the last position (zero-time).
                        
        self.dr = 50
        self.offset = 2500
        
        self.g_ns = 1
        self.g_gather_file = 'gather.bin'    # File to save the snapshots. size=nr*nt*n_snaps
        self.g_src_file= 'src.txt'     # Files that contain sources and geophones positions  
        self.g_rcv_file= 'rcv.txt'     # (see below their format)

        
#        self.absPath()
    
        
    def absPath (self):
#        print (self.vp_file)
#        print (self.rho_file)
#        print (self._gather_file)
#        print (self._src_file)
#        print (self._rcv_file)
#        print (self.wfl_file)  
#        
        self.vp_file= self.path + self.vp_file
        self.rho_file=self.path + self.rho_file
        self.g_gather_file = self.path + self.g_gather_file   
        self.g_src_file= self.path + self.g_src_file
        self.g_rcv_file= self.path + self.g_rcv_file 
        self.wfl_file = self.path + self.wfl_file
#
#        print (self.vp_file)
#        print (self.rho_file)
#        print (self._gather_file)
#        print (self._src_file)
#        print (self._rcv_file)
#        print (self.wfl_file)    

    
    def writeToFile (self, filename):
        f = open(filename, 'w')
        f.write ('ox=%s\n' % (self.ox))    
        f.write ('oz=%s\n' % (self.oz))    
        
        f.write ('nx=%s\n' % (self.nx))    
        f.write ('nz=%s\n' % (self.nz))    
        
        f.write ('dh=%s\n' % (self.dh))
        
        f.write ('vp=%s\n' % (self.vp_file))
        f.write ('rho=%s\n' % (self.rho_file))
        
        f.write ('nt=%s\n' % (self.nt))
        f.write ('dt=%s\n' % (self.dt))
        f.write ('nu0=%s\n' % (self.nu0))
        f.write ('out=%s\n' % (self._gather_file))
        
        f.write ('ns=%s\n' % (self._ns))
        f.write ('nr=%s\n' % (self._nr))
        f.write ('src=%s\n' % (self._src_file))
        f.write ('rcv=%s\n' % (self._rcv_file))
        
        f.write ('npml=%s\n' % (self.npml))
        f.write ('damp=%s\n' % (self.damp))
        f.write ('free_surf=%s\n' % (self.free_surf))
        f.write ('apprx=%s\n' % (self.apprx))
        
        f.write ('snap=%s\n' % (self.snap))
        f.write ('wfl=%s\n' % (self.wfl_file))
        
        f.write ('verb=%s\n' % (self.verb))
        f.write ('expl=%s\n' % (self.expl))
        f.write ('back=%s\n' % (self.back))        
        
        f.write ('dr=%s\n' % (self.dr)) 
        f.write ('offset=%s\n' % (self.offset)) 
        
        f.write ('g_ns=%s\n' % (self.g_ns)) 
        f.write ('g_gather_file=%s\n' % (self.g_gather_file)) 
        f.write ('g_src_file=%s\n' % (self.g_src_file)) 
        f.write ('g_rcv_file=%s\n' % (self.g_rcv_file)) 
        

    def readFromFile (self, filename):
        f = open(filename, 'r')
        d = {}
        for line in f:
            k = line.split('=')
            v = k[1]
            k = k[0]
            d[k]=v  
    
        self.ox = float(d['ox'])
        self.oz = float(d['oz'])

        self.nx = int(d['nx'])
        self.nz = int(d['nz'])

        self.dh = float(d['dh'])
        self.vp_file = d['vp'].rstrip()
        self.rho_file = d['rho'].rstrip()
        
        self.nt = int(d['nt'])
        self.dt = float(d['dt'])
        self.nu0 = float(d['nu0'])
        self._gather_file = d['out'].rstrip()


        self._ns = int(d['ns'])
        self._nr = int(d['nr'])
        self._src_file = d['src'].rstrip()
        self._rcv_file = d['rcv'].rstrip()

        self.npml = float(d['npml'])
        self.damp = float(d['damp'])
        self.free_surf = int(d['free_surf'])
        self.apprx = int(d['apprx'])

        self.snap = int(d['snap'])
        self.wfl_file = d['wfl'].rstrip()

        self.verb = int(d['verb'])
        self.expl = int(d['expl'])
        self.back = int(d['back'])

        if 'dr' in d:
            self.dr = float(d['dr'])
            
        if 'offset' in d:
            self.offset = float(d['offset'])
            
        if 'g_ns' in d:
            self.g_ns = int(d['g_ns'])
      
        if 'g_gather_file' in d:
            self.g_gather_file = d['g_gather_file'].rstrip()
            
        if 'g_src_file' in d:
            self.g_src_file = d['g_src_file'].rstrip()
            
        if 'g_rcv_file' in d:
            self.g_rcv_file = d['g_rcv_file'].rstrip()
            
#        print (self.vp_file)
#        print (self.rho_file)
#        print (self._gather_file)
#        print (self._src_file)
#        print (self._rcv_file)
#        print (self.wfl_file)  

    def lx(self):
        return (self.nx-1)*self.dh

    def lz(self):
        return (self.nz-1)*self.dh
        
    def setCurrentShot(self, shot):
        self._src_file = self.g_src_file + '_' + str(shot)
        self._rcv_file = self.g_rcv_file + '_' + str(shot)
        self._gather_file = self.g_gather_file + '_' + str(shot)
        
    def generateGeomFiles(self, shot) : 
        self.setCurrentShot(shot)
        
        f = open(self._src_file, 'w')
        lx = self.lx()
#        print ('lx', lx)
        xs = lx/2
        if (self.g_ns > 1):
            ds = lx/(self.g_ns-1)
#            print ('ds', ds)
            xs = self.ox + ds*shot;
        f.write ('%s %s\n' % (xs, 0))
        
        f = open(self._rcv_file, 'w')
        x_start = xs - self.offset
        x_end = xs + self.offset
        x_start = max(x_start,self.ox)
        x_end = min(x_end,lx)
        self._nr = (x_end - x_start)/self.dr + 1
        for i in range(self._nr):
            x = x_start + i * self.dr
            f.write ('%s %s\n' % (x, 0))
            
    def generateGeomFilesRoundModel(self) :
#        self._src_file = self._src_file
        self._rcv_file = self._rcv_file + '_round'
        lx = self.lx()
        lz = self.lz()
        sx = 0
        sz = 0
        
        f = open(self.src_file, 'w')
        f.write ('%s %s\n' % (lx/2, 0))
        
        self._nr = 0
        
        f = open(self._rcv_file, 'w')
        
#        for i in range(self.nx-1):
#            x = sx + i * self.dh
#            f.write ('%s %s\n' % (x, sz))
#            self._nr  += 1
            
#        for i in range(self.nx-1):
#            x = sx + i * self.dh    
#            f.write ('%s %s\n' % (x, lz))
#            self._nr  += 1
#            
        for j in range(self.nz-1):
            z = sz + j * self.dh
            f.write ('%s %s\n' % (sx, z))
            self._nr  += 1
            
        for j in range(self.nz-1):
            z = sz + j * self.dh
            f.write ('%s %s\n' % (lx, z))
            self._nr  += 1
            
    def generateGeomFilesEvgeny(self) : 
        f = open(self._src_file, 'w')
        f.write ('%s %s\n' % (self.nx*self.dh/2, 0))
        
        f = open(self._rcv_file, 'w')
        self._nr = self.nx
        for i in range(self._nr):
            x = i * self.dh
            f.write ('%s %s\n' % (x, 0))

    def readGather (self, shot, read = True):
        self.setCurrentShot(shot)
        
            
        g = gather(self.nt, self.dt, self.dr, self._getSourcePosition(), self._getRecPosition())
        if read:
            g.readValues (self._gather_file)
            
        return g
        
    def drawEnargyAtSource (self, figure_name):
        en_on_source = numpy.fromfile(self._gather_file + '_s', dtype=numpy.float32)

        plt.rcParams['figure.figsize'] = 40, 30
        plt.figure (figsize=(40, 30))
#        plt.rc('font', family='serif', size=60)
        ax = plt.subplot(111)
#        ax.set_title("Plan")
        plt.xlabel('Time (s)', **global_font)
        
        en_on_source = [abs(v) for v in en_on_source]
        x = [(self.nt-v)*self.dt for v in range(len(en_on_source))]
        
        ax.plot(x, en_on_source)
       
        plt.grid(True)
        plt.savefig(figure_name,bbox_inches='tight')
#        if show:
#            plt.show ()
        plt.gcf().clear()
        plt.close('all')
            
    def createImages (self):
        
        images_dir = self.path + 'images/'
        if not os.path.exists(images_dir):
            os.makedirs(images_dir)

        m  = model (self.nx, self.nz, self.dh, self.dh)
        m.readFromFile(self.vp_file)            
        if self.back == 0:
            m.draw ('', images_dir + 'forward_model.png', min_=1500, max_=5000)
        else:
            m.draw ('', images_dir + 'backward_model.png', min_=1500, max_=5000)
        
        norm = 1e-8
        g = self.readGather (0)
        if self.back == 0:
            g.draw ('', images_dir + 'forward_gather.png', norm = norm)
        else:
            g.draw ('', images_dir + 'backward_gather.png', norm = norm)
            
#        if self.back == 1:
#            self.drawEnargyAtSource(images_dir + 'energy_at_source.png')
##            g.draw ('', images_dir + 'backward_gather_on_s.png',norm=None)
        
        snaps = snapShot(self.wfl_file, m)
        
        if self.back == 0:
            snaps.draw (images_dir + 'forward_', norm = norm)
            
        if self.back == 1:
            snaps.draw (images_dir + 'backward_', norm = norm)
            
        final = snaps.spans[len(snaps.spans)-1]
        
        # take source position
#        [sx, sz] = self._getSourcePosition()
        
        if self.back == 1:
            final.draw ('test', figure_name = images_dir + 'final.png', cmap = 'gray', norm = norm)               
            zoom = final.zoom(g.s_pos[0], g.s_pos[1], 250, 250)
            zoom.draw ('zoom', figure_name = images_dir + 'final_zoom.png', cmap = 'gray', norm = norm)        
        
            
    def _getSourcePosition (self):
        assert(self._ns == 1)
#        print (self._src_file)
        f = open(self._src_file, 'r')
        l = f.readline()
        l = l.split (" ")
        sx = float (l[0])
        sz = float (l[1])
        
        return [sx, sz]

    def _getRecPosition (self):
        f = open(self._rcv_file, 'r')
        rec = []
        for line in f:
            k = line.split(' ')
            x = float(k[0])
            y = float(k[1])
            rec.append([x,y])
        return rec
        
     
def draw_convergence (func, xlabel, ylabel, label=None, figure_name = None, show=False, last = True):    
    import matplotlib.pyplot as plt
    
    plt.rcParams['figure.figsize'] = 40, 30
#        fig = plt.figure (figsize=(60, 15))
#    fig = plt.figure()
    plt.rc('font', family='serif', size=30)
    ax = plt.subplot(111)
    if label != None:
        ax.set_title(label)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    
    for i in range(len(func)):
        l = 5
        if i == len(func)-1 and last:
            l = 10
        ax.plot(func[i], linewidth=l)
    
#    plt.grid(True)
    if figure_name != None:
        plt.savefig(figure_name,bbox_inches='tight', dpi=100)
#            plt.savefig(figure_name,bbox_inches='tight')
        
    if show:
        plt.show ()
    plt.gcf().clear()
    plt.close('all')
   
            
#def modelingOneModel(c, model_type):
#        
##    c = config(model_path)
##    c.nz = 50
#        
#    # enable PML
##    c.free_surf = 0
##    c.npml=3
##    c.snap = c.nt
#
#    
#    vel = generateModel (c.nx, c.nz, c.dh, c.dh, model_type)
##    vel.draw ('', model_path + 'vel_gen.png')
#    
#    vel.writeToFile (c.vp_file)
#
#    rho = copy.deepcopy(vel)
#    rho.resetValues (1)
#    # use vel as roh
#    rho.writeToFile (c.rho_file)
#
#    c.generateGeomFiles()
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
#    log_name = c.path + 'log.txt'
#    with open(log_name, 'w') as f:
#        subprocess.call(args, stdout=f)
#        
#    c.createImages ()
#    
#    c.back = 1
#    # enable PML
##    c.free_surf = 0
##    c.npml=3
#    
#    c.wfl_file = param_name = c.path + 'snap_bw.bin'
#
##    # MUTE    
#    g = c.readGather ()
#    g = muteDirect (g, c.dh, -0.1, 2000)
#    g = muteDirect (g, c.dh, 0.16, 3500)
#    g = muteOffset (g, c.dh, 0, 1000)
##    
#
#    new_nather_name = c._gather_file + '_mute'
#    g.writeToFile (new_nather_name)    
#    c._gather_file = new_nather_name
#
#    
#    param_name = c.path + 'param_bw.txt'
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
#
##path = '/home/cloudera/TRM/acoustic_FD_TRM/tests/'
#    
#
#     
#def drawMod1(v1, v2, sl, figure_name):    
#    import matplotlib.pyplot as plt
#
#    fig = plt.figure (figsize=(40, 30))
#    plt.rc('font', family='serif', size=60)
#    ax = plt.subplot(111)
#    ax.set_title('Energy')
#    plt.ylabel('V2 (m/s)')
#    plt.xlabel('V1 (m/s)')
#        
#    cax = ax.imshow(sl)
#    
#    plt.xticks(range(len(v1)), v1)
#    plt.yticks(range(len(v2)), v2)
#        
#    fig.colorbar(cax)
#    plt.grid(True)
#    plt.savefig(figure_name,bbox_inches='tight', dpi=100)
##    plt.show ()
#    plt.gcf().clear()

#
#def modeling():
#
#    if not os.path.exists(path):
#        os.makedirs(path)
#    
#    for model_type in range (1, 6):
#    
#        model_path = path + 'model' + str (model_type) +'/'
#        if not os.path.exists(model_path):
#            os.makedirs(model_path)
#            
#        modelingOneModel(model_path, model_type)
#
#def inverse (c):
#    param_name = c.path + 'param_bw.txt'
#    c.readFromFile(param_name)    
#    
#    # only first and last snapshots are needed
#    c.snap = -1
#
#    m = model(c.nx, c.nz, c.dh, c.dh)
#    
#        
#    test_num = 0
#    
#    v1start = 1450
#    v2start = 3450
#    v1stop = 2550
#    v2stop = 4050
#    v1step = 10
#    v2step = 10
#    v1num = int((v1stop-v1start)/v1step) + 1
#    v2num = int((v2stop-v2start)/v2step) + 1
#
#    v1a = [v1start + i*v1step for i in range (v1num)]
#    v2a = [v2start + i*v2step for i in range (v2num)]
#    
#    matrix = numpy.zeros((v1num, v2num))
#    
#    if not os.path.exists(c.path + 'tests/' ):
#        os.makedirs(c.path + 'tests/' )
#    drawMod1(v1a, v2a, matrix, c.path + 'tests/' + 'energy.png')
#    
#    for i in range(v1num):
#        for j in range(v2num):
#            v1 = v1a[i]
#            v2 = v2a[j]
#
#
#            test_num += 1       
#            m = fillModel1(m, 125, v1, v2)
#            e = calcEnergy(c, m, test_num)    
#            print (v1, v2, e)
#
#            matrix [i][j]=e
#
#            drawMod1(v1a, v2a, matrix, c.path + 'tests/' + 'energy.png')           
#        
#def calcEnergy (c, vel, test_num = 0, image_path = None):
#    assert (c.back == 1)
#    assert (c.snap == -1)
#
#    test_name = c.path + 'tests/' + str (test_num) + '/' 
#    if not os.path.exists(test_name):
#        os.makedirs(test_name)
#            
#    c.vp_file = test_name + 'vel.bin'
#    c.wfl_file = test_name + 'snaps.bin'
#    param_name = test_name + 'param_bw.txt'
#    
#    # write model
#    vel.writeToFile (c.vp_file)
#
#    c.writeToFile(param_name)   
#    
#    # command line args
#    args = ['/home/cloudera/TRM/acoustic_FD_TRM/bin',
#            param_name]
#    
#    FNULL = open(os.devnull, 'w')
#    subprocess.call(args, stdout=FNULL, stderr=subprocess.STDOUT)
#        
#    # read snaps
#    snaps = snapShot(c.wfl_file, vel)
#    # only 1 snaps should be red
#    assert(len(snaps.spans) == 1)
#    
#    # take last snap
#    final = snaps.spans[0]
#    
#    # take source position
#    [sx, sz] = c._getSourcePosition()
##    print ('sx', sx, 'sz', sz)
#    
#    
#    [six, siz] = final.getIndex(sx, sz)
#    six = int (six)
#    siz = int (siz)
#
#    total_power = 0
#    for i in range (final.nx):
#        for j in range (final.nz):
#            total_power += final.v[i][j]**2
#
#    entropy = 0
#    for i in range (final.nx):
#        for j in range (final.nz):
#            p_norm = final.v[i][j]**2/total_power
#            entropy += p_norm*math.log(p_norm)
#            
#    entropy = -entropy
#    
#    energy = final.v[six][siz]**2
#    
##    for i in range (mask.nx):
##        for j in range (mask.nz):
##            [x,z] = final.getCoordByIndex (i, j)
##            dist = math.sqrt((sx- x)**2 + (sz-z)**2)
##            # dist 0 - taper = 1
##            # dist = area - taper = 0.5
##            taper = 1/(dist/area + 1)
##            mask.v[i][j] = taper**4
#
##    mask.draw ('', 'mask.png')
#    
#    #take energy in source point
##    energy = 0
##    for i in range (mask.nx):
##        for j in range (mask.nz):
##            e = final.v[i][j]**2
##            taper = mask.v[i][j]
##            energy += e*taper
##    
#    # just for test
#    if image_path != None:
#        final.draw ('test' + str (test_num), figure_name = image_path + '.png', cmap = 'gray', norm = 1e-7)        
#        vel.draw ('', image_path +'_model.png', min_ = 1500, max_=5000)
#        c.drawEnargyAtSource(image_path + '_energy_at_source.png')
#        
#        zoom = final.zoom(sx, sz, 250, 250)
#        zoom.draw ('zoom', figure_name = image_path + '_zoom.png', cmap = 'gray', norm = 1e-7)        
#    
#
#            
#    # detele scratch files
#    import shutil    
##    shutil.rmtree(test_name)
#    
#    return energy, entropy
#    
#    
#if __name__ == "__main__":
#    model_type = 1
##    model_path = path + 'model' + str (model_type) +'/'
#    model_path = '/home/cloudera/TRM/acoustic_FD_TRM/test/'
#    c = config(model_path)
#    c.readFromFile('/home/cloudera/TRM/acoustic_FD_TRM/test/1.txt')
#    c.generateGeomFilesEvgeny()
#    
#    modelingOneModel(c, model_type)
#
##    path = '/home/cloudera/TRM/acoustic_FD_TRM/tests/model1/'
##    c = config(path)
#        
##    inverse (c)
    