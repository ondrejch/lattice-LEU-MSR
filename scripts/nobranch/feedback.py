#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 17:47:14 2019

@author: agoluogl
"""

import converge
from lattice import Lattice
import mapscan
from salts import Salt
import os
import time
import matplotlib.pyplot as plt


class Feedback(object):
    'Determine feeback effects for a particular lattice'
    def __init__(self, salt:str='flibe', iters:int=11, enr:float=0.15, temp:bool=True, density:bool=False, sf:float=0.1, l:float=20.0):

        if sum([temp,density]) > 1: # Right now can only examine one variable at a time for feedback effects
            print("Only one variable can be True")
            sys.exit("Too many variables to examine")

        'Constructor with default values'
        self.l:float         = l          # Hex lattice size [cm]
        self.sf:float        = sf         # Fuel salt fraction
        self.salt_name:str   = salt       # Salt identifier
        self.case:bool       = density    # Doppler only (False) or Doppler and Density (True)
        self.dotemp:bool     = temp       # find feedback effects of temperature (True) or do not (False)
        self.tempstart:float = 900.0      # Default start of the range at 900K
#add option to choose temperature start
        self.tempstep:float  = 20.0       # Default step size between temperatures to examine feeback effects
        self.iters:int       = iters      # Number of files (steps) to compare
        self.tempend:float   = self.tempstart+self.tempstep*iters # End of Temp Range
        self.mainpath:str    ="~/git/lattice-LEU-MSR/scripts/ashley/" #change so that it gets whatever current directory

        self.k:float         = None          # ANALYTICAL_KEFF
        self.kerr:float      = None          # ANALYTICAL_KEFF Error
        self.cr:float        = None          # CONVERSION_RATIO
        self.crerr:float     = None          # CONVERSION_RATIO Error

        self.tempsK:list     = [0]*iters     # Salt temperatures to compare [K]
        self.deckdirs:list   = [0]*iters     # directory names
        self.qsubnames:list  = [0]*iters     # run.sh file names
        self.lattice:obj     = [0]*iters     # Lattice()
        self.keffs:list      = [0]*iters     # k-effective values
        self.kerrs:list      = [0]*iters     # errors in k-effective values
        self.crs:list        = [0]*iters     # conversion ratios
        self.crerrs:list     = [0]*iters     # errors in conversion ratios
        self.datapath:str    = [0]           # paths for data files based on options (T/F) & values (900K)

    def run(self):
        if self.dotemp:
            self.datapath="tempfeedbacks{:.0f}K_{:.0f}K.csv".format(self.tempstart,self.tempend)
            os.system("rm "+self.mainpath+"data/"+self.datapath)
            with open(self.datapath,"w") as df:
                df.write("Salt: {},Temp Range: {}K - {}K\n".format(self.salt_name,self.tempstart,self.tempend))
                df.write("K-Effective, Error, Conversation Ratio, Error\n")
            for i in range(self.iters):
                self.tempsK[i]=self.tempstart+self.tempstep*i
                self.deckdirs[i]="lat{}".format(self.tempsK[i])
                self.qsubnames[i]="run{}.sh".format(i)
                self.lattice[i]=Lattice()
                self.lattice[i].tempK=self.tempsK[i]
                self.lattice[i].deck_path=os.path.expanduser(self.mainpath+self.deckdirs[i])
                self.lattice[i].qsub_path=os.path.expanduser(self.mainpath+self.qsubnames[i])
                self.lattice[i].for_feedbacks() #writes and saves all of the necessary files, submits to cluster

    def getdata(self):
            for i in range(self.iters):
                self.lattice[i].get_calculated_values()
                self.keffs[i]=self.lattice[i].k
                self.kerrs[i]=self.lattice[i].kerr
                self.crs[i]=self.lattice[i].cr
                self.crerrs[i]=self.lattice[i].crerr
                with open(self.datapath,"a+") as df:
                    df.write("{},{},{},{}\n".format(self.keffs[i],self.kerrs[i],self.crs[i],self.crerrs[i]))
                self.lattice[i].cleanup()
            os.system("mkdir data; mv "+self.datapath+" "+self.mainpath+"data/"+self.datapath)




#     def plotdata(self):
#            for i in range(self.iters):





#add: if salt, if enr, if density, if l, if sf, etc...



#-----------------------------------------------------
if __name__ == '__main__':
#    print("This module handles feedback effects.")
    testing = Feedback()
    testing.run()
    testing.getdata()
#    print(testing.keffs)
