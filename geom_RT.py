#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 23:07:41 2016

@author: cloudera
"""


class GeomRec ():
    def __init__ (self, x=0, y=0, z = 0, c = 1):
        # c = 0 refraction
        # c = 1 reflection
        self.x = x
        self.y = y
        self.z = z
        self.c = c
        self.t = 0
        self.dt = 0
        return
    
    def writeToFile (self, f):
        f.write('r %s %s %s %s %s %s\n' % (self.x, self.y, self.z, self.c, self.t, self.dt))
     
    def readFromFile (self, f):
        line = f.readline().split()
        self.x = float(line [1])
        self.y = float(line [2])
        self.z = float(line [3])
        self.c = float(line [4])
        self.t = float(line [5])
        self.dt = float(line [6])
        
class GeomSource ():
    def __init__ (self, x=0, y =0, z = 0):
        self.x = x
        self.y = y
        self.z = z
        self.rec = []
        
    def addRec(self, rec):
        # c = 0 refraction
        # c = 1 reflection
        self.rec.append(rec)
        
    def generateRec (self, appX, appY, dx, dy, c = 0) :
        indX = int(appX/dx)
        indY = int(appY/dy) 
        for i in range (-indX,indX + 1):
            for j in range (-indY,indY + 1):
                x = self.x + i*dx
                y = self.y + j*dy
                self.addRec(GeomRec(x, y, self.z, c))
                
    def moveSource (self, x, y) :
        dx = x - self.x
        dy = y - self.y
        self.x = x
        self.y = y
        for r in self.rec:
            r.x = r.x + dx
            r.y = r.y + dy
        
    
    def writeToFile (self, f):
        f.write('s %s %s %s %s\n' % (self.x, self.y, self.z, len(self.rec)))
        for r in self.rec:
            r.writeToFile(f)
            
    def readFromFile (self, f):
        
        line = f.readline().split()
        self.x = float(line [1])
        self.y = float(line [2])
        self.z = float(line [3])
        
        num_rec = int (line [4])
        for i in range(num_rec):
            rec = GeomRec ()
            rec.readFromFile (f)
            self.addRec (rec)
            
class Geom ():
    def __init__ (self):
        self.sources = []
        return
        
    def addSource (self,s):
        self.sources.append(s)
        
    def writeToFile (self, filename):
        f = open(filename, 'w')
        f.write(str (len(self.sources))+str ('\n'))
        for s in self.sources:
            s.writeToFile(f)
    def readFromFile (self, filename):
        f = open(filename, 'r')
        line = f.readline().split()
        num_source = int (line[0])
        for i in range (num_source):
            source = GeomSource ()
            source.readFromFile (f)
            self.addSource(source)
            
def offset (s, r):
    import math
    dx = r.x - s.x
    dy = r.y - s.y
    return math.sqrt (dx**2 + dy**2)