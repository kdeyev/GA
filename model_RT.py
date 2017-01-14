#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 10:14:44 2016

@author: cloudera
"""

import numpy
import matplotlib.pyplot as plt

eps_dist = 0.001 

def printArray (a):
    a = [str (v) for v in a]
    return ' '.join(a)
    
def toIndex (v, start, delta):
    return (v-start)/delta
        
def toValue (v, start, delta):
    return v*delta+start
        
def getDelta (array):
    delta = array[1]-array[0]
    for i in range(1, len(array)):
        cd = array[i]-array[i-1]
        assert(abs(delta - cd) < eps_dist)
    return delta
    
   
class VelModel ():
    def __init__ (self):
        return
        
    def loadFromFile (self, filename):
        rows = 0
        f = open ( filename , 'r')
        line = f.readline().split()
        rows = rows + 1
        nx = int (line[0])
        ny = int (line[1])
        nz = int (line[2])
        Vw = float (line[3])
        Va = float (line[4])
        self.x_nodes = f.readline().split()
        rows = rows + 1
        self.x_nodes = [float(v) for v in self.x_nodes]
        self.y_nodes = f.readline().split()
        rows = rows + 1
        self.y_nodes = [float(v) for v in self.y_nodes]

        assert (nx == self.nx ())
        assert (ny == self.ny ())
        
        self.topo = []
        for i in range(self.nx ()):
            nodes = f.readline().split()
            rows = rows + 1
            nodes = [float(v) for v in nodes]
            self.topo.append(nodes)
        
        self.z_nodes = f.readline().split()
        rows = rows + 1
        self.z_nodes = [float(v) for v in self.z_nodes]
        assert (nz == self.nz ())
        
        self.v = numpy.loadtxt(filename, skiprows=rows )
        self.v = numpy.reshape(self.v, (nx, ny, nz))
#        self.v = numpy.transpose(self.v, (0,2,1))

        self.dx = getDelta (self.x_nodes) 
        self.dy = getDelta (self.y_nodes) 
        self.dz = getDelta (self.z_nodes)

    def x (self, index) :
        return self.x_nodes[index]
    def y (self, index) :
        return self.y_nodes[index]
    def z (self, index) :
        return self.z_nodes[index]
        
    def nx (self) :
        return len(self.x_nodes)
    def ny (self) :
        return len(self.y_nodes)
    def nz (self) :
        return len(self.z_nodes)
        
    def lx (self) :
        return self.x_nodes[len(self.x_nodes)-1]
    def ly (self) :
        return self.y_nodes[len(self.y_nodes)-1]
    def lz (self) :
        return self.z_nodes[len(self.z_nodes)-1]
        
    def resetValues (self, v):
      for i in range (self.nx()):
        for j in range (self.ny()):
            for k in range (self.nz()):
                self.v[i][j][k] = v      

    def moveModel (self, dx, dy) :
        self.x_nodes = [v + dx for v in self.x_nodes]
        self.y_nodes = [v + dy for v in self.y_nodes]
        

    def setValue(self, x, y, z, v):
        i = self.indexX (x)
        j = self.indexY (y)
        k = self.indexZ (z)
#        assert (i.is_integer())
#        assert (j.is_integer())
#        assert (k.is_integer())
        i = round (i)
        j = round (j)
        k = round (k)
        self.v[i][j][k] = v

    def getValue(self, x, y, z):
#        return numpy.interp3(x,y,z,self.v,self.x_nodes,self.y_nodes,self.z_nodes)
        i = self.indexX (x)
        j = self.indexY (y)
        k = self.indexZ (z)
#        assert (i.is_integer())
#        assert (j.is_integer())
#        assert (k.is_integer())
        i = round (i)
        j = round (j)
        k = round (k)

        # hack
        if i == 1 and self.nx () == 1:
            i = 0

        # hack
        if j == 1 and self.ny () == 1:
            j = 0

        if i < 0 or i >= self.nx ():
            return None
            
        if j < 0 or j >= self.ny ():
            return None
        
        if k < 0 or k >= self.nz ():
            return None
            
        return self.v[i][j][k]

    def getInterpValue(self, x, y, z):
        if x < 0 or x > self.lx ():
            return None
            
        if y < 0 or y > self.ly ():
            return None
        
        if z < 0 or z > self.lz ():
            return None

        from scipy.interpolate import RegularGridInterpolator
        fn = RegularGridInterpolator((self.x_nodes,self.y_nodes,self.z_nodes), self.v)
        return fn([[x,y,z]])[0]
        
        i = self.indexX (x)
        j = self.indexY (y)
        k = self.indexZ (z)

            
        x0 = round (i)
        y0 = round (j)
        z0 = round (k)
        
        x = x - x0
        y = y - y0
        z = z - z0
        
        x1 = x0 + 1
        y1 = y0 + 1
        z1 = z0 + 1

        output = (
        self.v[x0][y0][z0]*(1-x)*(1-y)*(1-z) +
        self.v[x1][y0][z0]*x*(1-y)*(1-z) +
        self.v[x0][y1][z0]*(1-x)*y*(1-z) +
        self.v[x0][y0][z1]*(1-x)*(1-y)*z +
        self.v[x1][y0][z1]*x*(1-y)*z +
        self.v[x0][y1][z1]*(1-x)*y*z +
        self.v[x1][y1][z0]*x*y*(1-z) +
        self.v[x1][y1][z1]*x*y*z)

        return output

        
    def indexX (self, v):
        return toIndex (v, self.x_nodes[0], self.dx)
    def indexY (self, v):
        return toIndex (v, self.y_nodes[0], self.dy)
    def indexZ (self, v):
        return toIndex (v, self.z_nodes[0], self.dz)
        
    def valueX (self, v):
        return toValue (v, self.x_nodes[0], self.dx)
    def valueY (self, v):
        return toValue (v, self.y_nodes[0], self.dy)
    def valueZ (self, v):
        return toValue (v, self.z_nodes[0], self.dz)
        
    def emptyModel (self, lx, ly, lz, dx, dy, dz):
        nx = int (lx/dx) + 1
        ny = int (ly/dy) + 1
        nz = int (lz/dz) + 1
        
#        assert (lx == nx*dx)
#        assert (ly == ny*dy)
#        assert (lz == nz*dz)
        
        
        self.x_nodes = [toValue(i, 0, dx) for i in range (nx)]
        self.y_nodes = [toValue(i, 0, dy) for i in range (ny)]
        self.z_nodes = [toValue(i, 0, dz) for i in range (nz)]

        self.topo = numpy.zeros((nx, ny))
        self.v = numpy.zeros((nx, ny, nz))
        
        self.dx = dx#getDelta (self.x_nodes) 
        self.dy = dy#getDelta (self.y_nodes) 
        self.dz = dz#getDelta (self.z_nodes)

        
    def writeToFile (self, filename):
        f = open(filename, 'w')
        f.write ('%s %s %s %s %s\n' % (len (self.x_nodes), len (self.y_nodes), len (self.z_nodes), 1.5, 0.33))
        f.write (printArray (self.x_nodes) + str ('\n'))
        f.write (printArray (self.y_nodes) + str ('\n'))
        for t in self.topo:
            f.write (printArray(t) + str ('\n'))
        f.write (printArray(self.z_nodes) + str ('\n'))
        
        v = self.v
        v = numpy.reshape(v, (self.nx () * self.ny (), self.nz ()))
        for t in v:
            f.write (printArray(t) + str ('\n'))
        
        
    def print_ (self):
        print (self.v)
        print (self.x_nodes)
        print (self.z_nodes)
        
    
def drawProjX (model, showY, ax, fig, cmap, norm = None):
    proj_x = numpy.zeros((model.nx(), model.nz ()))
    y_index = int(model.indexY (showY))
    
    for i in range (model.nx()):
        for k in range (model.nz()):             
            proj_x[i][k] = model.v[i][y_index][k]

    proj_x = numpy.transpose(proj_x)
    
    a = model.dz/model.dz
#    a = 1
    a = 'auto'
    if norm != None:
        cax = ax.imshow(proj_x, aspect=a, vmin=-norm, vmax=norm, cmap=cmap)
    else:
        cax = ax.imshow(proj_x, aspect=a, cmap=cmap)
    
#    fig.colorbar(cax, orientation='horizontal')
    fig.colorbar(cax)
    
                    
   
def drawProjZ (model, showX, showY, figure_name, show = False):
                 
    import matplotlib.pyplot as plt

#    plt.rcParams['figure.figsize'] = 40, 30
    plt.figure (figsize=(40, 30))
    plt.rc('font', family='serif', size=60)
    ax = plt.subplot(111)
#    ax.set_title(title)
    
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Depth (m)')
    plt.gca().invert_yaxis() 
    
    x_index = int(model.indexX (showX))
    y_index = int(model.indexY (showY))

    Z = []
    V = []
        
    for k in range (model.nz()):
        z = model.z(k)*1000
        v = model.v[x_index][y_index][k]*1000
        Z.append(z)
        V.append(v)

    Z1 = []
    V1 = []
    Z1.append(0)
    V1.append(1500)       
    Z1.append(0)
    V1.append(4000)
        
        
    ax.plot(V, Z)
    ax.plot(V1, Z1)
    
    
    if figure_name != None:
        plt.savefig(figure_name,bbox_inches='tight', dpi=100)

    if show:
        plt.show()
    plt.gcf().clear()
