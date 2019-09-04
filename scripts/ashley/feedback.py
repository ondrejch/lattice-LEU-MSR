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
"""

#Questions
#where is the lattice pitch being used? its not a hexagonal array of fuel things
#but is it chemical? how


test=Lattice() #makes default lattice
#print(test.s) #prints the salt that is in the lattice
#print(test.get_deck()) #prints the cell cards for the lattice in serpent

test.all() #makes input files, saves them, runs them, and gets data

print(test.k)
keff=test.k
kefferr=test.kerr
conversion_ratio=test.cr
crerr=test.crerr


test.cleanup()

"""

#%%
class Feedback(object):
    'Determine feeback effects for a particular lattice'
    def __init__(self, salt:str='flibe', iters:int=11, enr:float=0.15, temp:bool=True, density:bool=False, sf:float=0.1, l:float=20.0):

        if sum([temp,density]) > 1: # Right now can only examine one variable at a time for feedback effects
            print("Only one variable can be True")
            sys.exit("Too many variables to examine")

        'Constructor with default values'
        self.l:float         = l          # Hex lattice size [cm]
        self.sf:float        = sf         # Fuel salt fraction
#        self.r:float         = self.r()   # Diameter of the fuel salt channel [cm]
        self.salt_name:str   = salt       # Salt identifier
#        self.s               = Salt(self.salt_formula, enr) # Salt used
        self.case:bool       = density    # Doppler only (False) or Doppler and Density (True)
        self.dotemp:bool     = temp       # find feedback effects of temperature (True) or do not (False)
        self.tempbase:float  = 900.0      # Default middle of range at 900K
        self.iters:int       = iters      # Number of files (steps) to compare
        self.tempstep:float  = 20.0       # Default step size between temperatures to examine feeback effects

        self.k:float       = None       # ANALYTICAL_KEFF
        self.kerr:float    = None       # ANALYTICAL_KEFF Error
        self.cr:float      = None       # CONVERSION_RATIO
        self.crerr:float   = None       # CONVERSION_RATIO Error
        self.tempsK:list   = [0]*iters  # Salt temperatures to compare [K]
        self.deckdirs:list = [0]*iters  # directory names
        self.tests:list    = [0]*iters  # tests

    def do_things(self):
        if self.dotemp:
            self.tempstart=self.tempbase - self.iters*self.tempstep/2
            for i in range(self.iters):
                self.tempsK[i]=self.tempstart+self.tempstep*i
                self.deckdirs[i]="lat{}".format(self.tempsK[i])
                os.system("mkdir " + self.deckdirs[i])
                self.tests[i]=Lattice()
                self.tests[i].tempK=self.tempsK[i]
                self.tests[i].deck_path=os.path.expanduser("~/git/salt-management-DMSR/ashley/"+self.deckdirs[i])

"""

#add: if salt, if enr, if density, if l, if sf, etc...

    def call_Lattice(self):
        print(Lattice(self).k)          #check correctness
        self.k=Lattice(self).k          #k-effective value
        self.kerr=Lattice(self).kerr    #error in k-effective
        self.cr=Lattice(self).cr        #conversion ratio
        self.crerr=Lattice(self).crerr  #error in conversion ratio
        Lattice(self).cleanup()

"""



#-----------------------------------------------------
if __name__ == '__main__':
    print("This module handles feedback effects.")
    testing = Feedback()
    testing.do_things()
