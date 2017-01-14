#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 11:40:57 2016

@author: cloudera
"""
import numpy
import math
import matplotlib.pyplot as plt

class Ray ():
    def __init__(self):
        self.vals = []

    def readFromFile (self, f):
        for line in f:
            if line[0] == ">":
                return
            line = line.split()
            line = [float(v) for v in line]
            self.vals.append(line)
            
    def getX(self):
        return numpy.transpose(self.vals)[0]

    def getY(self):
        return numpy.transpose(self.vals)[1]

    def getZ(self):
        return numpy.transpose(self.vals)[2]

    def calcLen (self):
        l = 0
        for i in range (len(self.vals)-1):
            l += math.sqrt ((self.vals[i][0]-self.vals[i+1][0])**2 +
                             (self.vals[i][1]-self.vals[i+1][1])**2 +
                              (self.vals[i][2]-self.vals[i+1][2])**2
            )
        return l
            
class Rays (): 
    def __init__(self):
        self.rays = []

    def readFromFile (self, filename):
        f = open(filename)
        line = f.readline()
        for line in f:
            r = Ray()
            r.readFromFile(f)
            x = r.getX ()
            y = r.getY ()    
            z = r.getZ ()
            if z[0] != 0 or z[len(z)-1] != 0:
                continue
            self.rays.append(r) 
    

def drawRaysY(m, rs, ytoShow, ax, draw_points = False):
    for ray in rs.rays:
        x = ray.getX ()
        y = ray.getY ()    
        z = ray.getZ ()
        X = []
#        Y = []
        Z = []
        for i in range (len (x)):
#            if abs(y[i] - ytoShow) < 0.025 : 
                X.append(x[i])
#                Y.append(y[i])
                Z.append(z[i])
        
#        if (Z[0] > 0.010 or Z[len(Z)-1] > 0.010) and draw_points == True:
#            continue
        
        X = [m.indexX(v) for v in X]
        Z = [m.indexZ(v) for v in Z]
        ax.plot(X, Z)
        if draw_points:
            ax.plot(X[len(X)-1], Z[len(Z)-1], 'rv', markersize=20)
            ax.plot(X[0], Z[0], 'b^', markersize=20)
        
def drawRaysXY(m, rs, figure_name, show = False, draw_points = False):
    plt.rcParams['figure.figsize'] = 40, 30
    plt.figure (figsize=(40, 30))
    plt.rc('font', family='serif', size=60)
    ax = plt.subplot(111)
    ax.set_title("Plan")
    plt.ylabel('Location (m)')
    plt.xlabel('Location (m)')
    
#    label_x = [int (v*1000) for v in m.x_nodes]
#    label_z = [int (v*1000) for v in m.z_nodes]
#    plt.xticks(range(len(label_x)), label_x, fontsize=12)
#    plt.yticks(range(len(label_z)), label_z, fontsize=12)
    
    for ray in rs.rays:
        x = ray.getX ()
        y = ray.getY ()          
#        z = ray[2]

        x = [int(v*1000) for v in x]
        y = [int(v*1000) for v in y]
        ax.plot(x, y)
        if draw_points:
            ax.plot(x[len[x]-1], y[len[y]-1], 'b^', markersize=20)
            ax.plot(x[0], y[0], 'rv', markersize=20)
            
    plt.savefig(figure_name,bbox_inches='tight')
    if show:
        plt.show ()
    plt.gcf().clear()