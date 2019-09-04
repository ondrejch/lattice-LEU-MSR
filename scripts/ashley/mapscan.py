#!/usr/bin/python3
#
# Class to find critical enrichment of a lattice
#
# Ondrej Chvala, ochvala@utk.edu
# 2019-08-17
# GNU/GPL

import math
import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from concurrent import futures
import threading
import time
import lattice
import converge

my_debug:int = 5

#SALT_FRACTIONS  = [0.07,0.08,0.09,0.10,0.11,0.12,0.13]
#LATTICE_PITCHES = [10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0]

SALT_KEYS = ['flibe', 'lif', 'naf', 'nafbe12', 'nafbe30', 'nafrbf2', 'nafkf'] # list(lattice.SALTS.keys())
SALT_FRACTIONS  = [0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.18,0.20,0.225,0.25,0.275,0.30,0.325,0.35,0.375,0.40,0.425,0.45,0.475,0.50,0.525,0.55]
LATTICE_PITCHES = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0,25.0,26.0,28.0,30.0,32.0,34.0,36.0,38.0,40.0,45.0,50.0,55.0,60.0]

class LatticeData(object):
    'Class holding the calculated data'
    def __init__(self, salt:str='flibe', sf:float=0.1, l:float=20.0):
        self.salt:str       = salt      # Salt key
        self.sf:float       = sf        # Salt fraction in the lattice
        self.l:float        = l         # Lattice pitch [cm]
        self.enr:float      = None      # Uranium enrichment
        self.rho:float      = None      # Reactivity [pcm]
        self.rhoerr:float   = None      # sigma_{rho} [pcm]

    def __repr__(self):
        result = "Salt: %s, sf: %5.3f, l: %5.3f [cm]" % (self.salt, self.sf, self.l)
        if self.enr:
            result += "  enr: %6.4f, rho: %6.2f +- %6.2f pcm" % (self.enr, self.rho, self.rhoerr)
        return result

    def testfunction(self) -> float:
        'Development testing function'
        time.sleep(1)
        print("* test * ", threading.current_thread(), self)
        return self.sf*self.l

class MyTest(object):
    'Development testing class'
    def __init__(self):
        self.data = []

    def runtest(self):
        with futures.ThreadPoolExecutor(max_workers = 3) as executor:
            to_do = []
            for pitch in LATTICE_PITCHES:
                lat = LatticeData()
                lat.l = pitch
                future = executor.submit(lat.testfunction)
                to_do.append(future)

            for future in futures.as_completed(to_do):
                res = future.result()
                msg = '{} result: {!r}'
                print(msg.format(future, res))
                self.data.append(res)
        print(self.data)

class ScanConverge(object):
    'Go over sf/l phase space and coverge enrichments'
    def __init__(self, salt='flibe', sf_list=None, l_list=None):
        try:
           salt_formula = lattice.SALTS[salt]
        except:
            ValueError("Salt "+salt+" is undefined.")
        self.salt = salt
        self.sf_list = []
        if sf_list:
            self.sf_list = sf_list
        else:
            self.sf_list = SALT_FRACTIONS
        self.l_list = []
        if l_list:
            self.l_list = l_list
        else:
            self.l_list = LATTICE_PITCHES
        self.data = []

        # Increase convergence by using old values to start regula falsi
        # https://stackoverflow.com/questions/29974122/interpolating-data-from-a-look-up-table#30057858
        old_LUT = np.genfromtxt("/home/ondrejch/L/old_"+salt+".dat", delimiter=" ")
        self.LUTxy = old_LUT[:, :2]
        self.LUTval= old_LUT[:, 2]
        del old_LUT
        self.old_tree = cKDTree(self.LUTxy)

    def doconverge(self, sf, pitch) -> LatticeData:
        'Converge one lattice'
        c = converge.Converge(self.salt, sf, pitch)
        if not c.read_rhos_if_done():    # Was the enrichment not found already?
            xy = (sf, pitch)
            dist, ind = self.old_tree.query(xy, k=2) # Find nearest old enrichments
            d1, d2 = dist.T                     # Distance from our point
            v1, v2 = self.LUTval[ind].T         # Value - enrichment
            v = (d1)/(d1 + d2)*(v2 - v1) + v1   # Linear interpolation
            c.enr_min = v *0.6                  # Set regula falsi min
            c.enr_max = v *1.5                  #                  max
            c.iterate_rho()                     # Start iterations
            c.save_iters()
        res = LatticeData(self.salt, sf, pitch)
        res.enr    = c.conv_enr
        res.rho    = c.conv_rho
        res.rhoerr = c.conv_rhoerr
        if my_debug:
            print("* DBG: ", repr(res))
        return res

    def runscan(self):
        'Threaded convergence scan for sf x pitch phase space'
        with futures.ThreadPoolExecutor(max_workers=100) as executor:
            to_do = []
            for sf in self.sf_list:
                for pitch in self.l_list:
                    future = executor.submit(self.doconverge, sf, pitch)
                    to_do.append(future)

            for future in futures.as_completed(to_do):  # <7>
                res = future.result()  # <8>
                if my_debug:
                    msg = '{} result: {!r}'
                    print(msg.format(future, res)) # <9>
                self.data.append(res)
        if my_debug:
            print(self.data)

# ------------------------------------------------------------
if __name__ == '__main__':
    print("This module handles phase space scanning.")
#    input("Press Ctrl+C to quit, or enter else to test it.")
#    m = MyTest()
#    m.runtest()
    m = ScanConverge()
    m.runscan()


