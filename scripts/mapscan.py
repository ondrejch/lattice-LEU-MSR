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

my_debug:int = 1

SALT_FRACTIONS  = [0.07,0.08]
LATTICE_PITCHES = [10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0]

SALT_KEYS = ['flibe', 'lif', 'naf', 'nafbe12', 'nafbe30', 'nafrbf2', "nafzrf", 'nafkf'] # list(lattice.SALTS.keys())
SALT_FRACTIONS  = [0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.18,0.20,0.225,0.25,0.275,0.30,0.325,0.35,0.375,0.40,0.425,0.45,0.475,0.50,0.525,0.55]
LATTICE_PITCHES = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0,25.0,26.0,28.0,30.0,32.0,34.0,36.0,38.0,40.0,45.0,50.0,55.0,60.0]

class ConvergedPoint(object):
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

class ScanConverge(object):
    'Go over sf/l phase space and coverge enrichments'
    def __init__(self, salt='flibe', sf_list=None, l_list=None):
        try:
           salt_formula = lattice.SALTS[salt]
        except:
            raise ValueError("Salt "+salt+" is undefined.")
        self.salt      = salt       # Salt key
        self.conv_list = []         # List of convergence objects, 1 per thread
        self.data      = []         # List of ConvergedPoint results 
        self.sf_list   = []         # Salt fractions to scan
        if sf_list:
            self.sf_list = sf_list
        else:
            self.sf_list = SALT_FRACTIONS
        self.l_list    = []         # Channel pitches to scan
        if l_list:
            self.l_list = l_list
        else:
            self.l_list = LATTICE_PITCHES

        # Increase convergence by using old values to start regula falsi
        # https://stackoverflow.com/questions/29974122/interpolating-data-from-a-look-up-table#30057858
        old_LUT = np.genfromtxt("/home/ondrejch/L/old_"+salt+".dat", delimiter=" ")
        self.LUTxy = old_LUT[:, :2]
        self.LUTval= old_LUT[:, 2]
        del old_LUT
        self.old_tree = cKDTree(self.LUTxy)

    def doconverge(self, c) -> ConvergedPoint:
        'Converge one lattice'
        tl = threading.local()          # Prevent threads overwriting each others data
        tl.res = ConvergedPoint(self.salt, c.sf, c.l)
        tl.is_converged:bool = False
        if not c.read_rhos_if_done():   # Was the enrichment not found already?
        #    c.cleanup_force_all()       # Wipe the directory
            tl.xy = (c.sf, c.l)
            tl.dist, tl.ind = self.old_tree.query(tl.xy, k=2) # Find nearest old enrichments
            tl.d1, tl.d2 = tl.dist.T                     # Distance from our point
            tl.v1, tl.v2 = self.LUTval[tl.ind].T         # Value - enrichment
            tl.v = (tl.d1)/(tl.d1 + tl.d2)*(tl.v2 - tl.v1) + tl.v1   # Linear interpolation
            c.enr_min = tl.v *0.5               # Set regula-falsi min
            c.enr_max = tl.v *1.8               #                  max
            if c.enr_max > 0.99:                # Sanity check
                c.enr_max = 0.99
            tl.is_converged = c.iterate_rho()   # Run iterations
            if tl.is_converged:
                c.save_iters()
        else:
            tl.is_converged = True
        if tl.is_converged:
            tl.res.enr    = c.conv_enr
            tl.res.rho    = c.conv_rho
            tl.res.rhoerr = c.conv_rhoerr
        else:
            tl.res.enr    = -1.0
            tl.res.rho    = -1.0
            tl.res.rhoerr = -1.0

        if my_debug > 3:
            print("* DBG: ", repr(tl.res))
        return tl.res

    def runscan(self):
        'Threaded convergence scan for sf x pitch phase space'
        with futures.ThreadPoolExecutor(max_workers=50) as executor:
            to_do = []
            for sf in self.sf_list:
                for l in self.l_list:
                    if not self.is_converged(sf, l):
                        self.conv_list.append(converge.Converge(self.salt, sf, l))   # Each point is a class on a list
                        future = executor.submit(self.doconverge,self.conv_list[-1]) # -1: last on the list 
                        to_do.append(future)
                        time.sleep(0.5)

            for future in futures.as_completed(to_do):
                res = future.result()
                if my_debug:
                    msg = '{} result: {!r}'
                    print(msg.format(future, res))
                self.data.append(res)                   # Record all results into a list
        if my_debug:
            print(self.data)

    def is_converged(self, sf, l):
        'Is this lattice already converged?'
        if not self.data:
            return False
        for d in self.data:
            if d.sf == sf and d.l == l:
                return True
        return False

    def read_data(self, savefile=None) -> bool:
        'Read converged lattices from savefile'
        if not savefile:
            c = converge.Converge(self.salt)    # Temp class just to get the main_path
            savefile = c.main_path + "_converged.dat"
        try:
            fh = open(savefile, 'r')
            for myline in fh.readlines():
                myline  = myline.strip().split()
                sf     = float(myline[0])
                l      = float(myline[1])
                cpt        = ConvergedPoint(self.salt, sf, l)
                cpt.enr    = float(myline[2])
                cpt.rho    = float(myline[3])
                cpt.rhoerr = float(myline[4])
                self.data.append(cpt)
            fh.close()
            return True
        except:
            return False

    def save_data(self, savefile=None) -> bool:
        'Save converged lattices to savefile'
        if not self.data:
            print("ERROR: No data to save!")
            return False
        if not savefile:
            if self.conv_list:
                savefile = self.conv_list[0].main_path + "_converged.dat"
            else:
                c = converge.Converge(self.salt)    # Get the path
                savefile = c.main_path + "_converged.dat"
        try:
            fh = open(savefile, 'w')
            for d in self.data:
                fh.write("%8.6f\t%8.5f\t%14.12f\t%10.2f %6.1f\n" % \
                    (d.sf, d.l, d.enr, d.rho, d.rhoerr) )
            fh.close()
        except IOError as e:
            print("[ERROR] Unable to write to file: ", \
                   self.deck_path + '/' + self.deck_name)
            print(e)
            return False
        
# ------------------------------------------------------------
if __name__ == '__main__':
    print("This module handles phase space scanning.")
#    input("Press Ctrl+C to quit, or enter else to test it.")
    m = ScanConverge()
    m.read_data()
    m.runscan()
    m.save_data()



