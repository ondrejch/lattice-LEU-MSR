#!/usr/bin/python3
#
# Class to find critical enrichment of a lattice
#
# Ondrej Chvala, ochvala@utk.edu
# 2019-08-16
# GNU/GPL

import math
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
import threading
import sys
import os
import shutil
import time

from lattice import Lattice

my_debug:int = 0

def rho(k:float) -> float:
    'K to rho [pcm]'
    return 1e5*(k-1.0)/k

class Converge(object):
    'Converge enrichment for a particular lattice'
    def __init__(self, salt:str='flibe', sf:float=0.1, l:float=20.0):
        'Constructor with default values'
        self.sf:float       = sf        # Fuel salt fraction
        self.l:float        = l         # Hex lattice size [cm]
        self.salt:str       = salt      # Salt key
        self.rho_tgt:float  = 200.0     # Target reactivity [pcm]
        self.eps_rho:float  = 200.0     # Reactivity epsilon [pcm]
        self.eps_enr:float  = 1e-9      # Epsilon enrichment
        self.main_path:str = os.path.expanduser('~/L/')+salt # Main path
        self.iter_path:str = self.main_path + "/"+"%08.6f"%self.sf + \
                "/%08.5f"%self.l        # Where the lattice decks are

        self.RhoData = namedtuple("RhoData" ,"enr rho rho_err")
        self.rholist        = []        # List of results
        self.sleep_sec:int  = 30        # Sleep timer between results read attempts [s]
        # Constants for iteration boundaries
        self.enr_min:float  = 0.007     # LEU enrichment boundaries
        self.enr_max:float  = 0.21      #  for iterative search
        self.iter_max:int   = 20        # Maximum # of iterations
        self.conv_enr:float = None      # Converged value of enrichment
        self.conv_rho:float = None      # Converged value of rho [pcm]
        self.conv_rhoerr:float = None   # Converged value of sigma_rho [pcm]
        self.force_recalc:bool = False  # Force recalculation of existing data

    def __repr__(self):
        'Pretty printing'
        result = '''Convergence for lattice with %s, sf: %5.3f, l: %5.3f cm
Iteration boudaries %5.4f %5.4f, max iters: %d''' % \
        (self.salt, self.sf, self.l, self.enr_min, self.enr_max, self.iter_max)
        if self.rholist:
            result += '\nList of rho(enr) found during convergece:'
            for r in self.rholist:
                result += '\n%10.8f   ->  %10.2f  +- %6.1f' % (r)
        if self.conv_enr:
            result += "\nRESULT: enr %10.8f, rho %10.2f +- %6.1f" % \
                        (self.conv_enr, self.conv_rho, self.conv_rhoerr)
        return result

    def iterate_rho(self) -> bool:
        'Execute the convergence search https://en.wikipedia.org/wiki/Regula_falsi#Example_code'
        # Create and run the edge points
        rho0:float = +1.0
        rho1:float = -1.0
        enr0:float = 0.0
        enr1:float = 0.0
        rho0err:float = 0.0
        rho0err:float = 0.0

        while rho0 > 0.0 or rho1 < 0.0:
            lat0 = Lattice(self.salt, self.sf, self.l, self.enr_min)
            lat0.set_path_from_geometry()
            lat1 = Lattice(self.salt, self.sf, self.l, self.enr_max)
            lat1.set_path_from_geometry()

            lat0.save_qsub_file() # As long as we are using defaults, shared qsub file will work
            if self.force_recalc or not lat0.get_calculated_values():
                lat0.cleanup()
                lat0.save_deck()
                lat0.run_deck()
            if self.force_recalc or not lat1.get_calculated_values():
                lat1.cleanup()
                lat1.save_deck()
                lat1.run_deck()
            is_done:bool = False            # Wait for Serpent
            while not is_done:
                if lat0.get_calculated_values() and lat1.get_calculated_values():
                    is_done = True
                else:
                    if my_debug:
                        print("[DEBUG CONV] sleeping ...")
                    time.sleep(self.sleep_sec)  # Wait a minute for Serpent ...

            rho0 = rho(lat0.k)  # [pcm]
            rho1 = rho(lat1.k)  # [pcm]           
            if self.enr_max > 0.98 and rho1 < 0.0:
                print("ERROR: lattice " + repr(lat1) +" cannot get critical.")
                return False
            enr0 = lat0.s.enr
            enr1 = lat1.s.enr
            rho0err = 1e5*lat0.kerr
            rho1err = 1e5*lat1.kerr
            if rho0 > 0.0:
                self.enr_min /= 1.5
            if rho1 < 0.0:
                self.enr_max *= 1.5
                if self.enr_max > 0.99:     # Sanity check
                    self.enr_max = 0.99

        self.rholist.append( self.RhoData(enr0, rho0, rho0err) )
        self.rholist.append( self.RhoData(enr1, rho1, rho1err) )
        
        if my_debug > 2:
            print(repr(self.rholist))

        # Regula Falsi root search, Illinois algorithm
        n_iter:int      = 0
        side:int        = 0
        enri:float      = None
        rhoi:float      = None
        rhoierr:float   = None
        while n_iter < self.iter_max:
            n_iter += 1
            if my_debug > 2:
                print("[DEBUG RF] ", rho0, enr0 , rho1, enr1)
            # Enrichment value for this iteration
            d_rho = rho0 - rho1
            if d_rho == 0.0:
                print("ERROR: div by zero", self.salt, self.sf, self.l, enri)
                return False
            enri = ((rho0-self.rho_tgt)*enr1 - (rho1-self.rho_tgt)*enr0) / d_rho
            if my_debug:
                print("[DEBUG RF] new enr: ", enri)
            if abs(enr1 - enr0) < self.eps_enr*abs(enr0+enr1):
                break   # Enrichment values close to each other, done
            # New lattice run
            l = Lattice(self.salt, self.sf, self.l, enri)
            l.set_path_from_geometry()
            if self.force_recalc or not l.get_calculated_values():
                l.cleanup()
                l.save_deck()
                l.run_deck()
                while not l.get_calculated_values():
                    if my_debug:
                        print("[DEBUG RF] sleeping ...")
                    time.sleep(self.sleep_sec)  # Wait a minute for Serpent ...
            rhoi   = rho(l.k)       # [pcm]
            rhoierr = 1e5*l.kerr     # [pcm]
            self.rholist.append( self.RhoData(enri, rhoi, rhoierr) )
            if (rhoi-self.rho_tgt)*(rho1-self.rho_tgt) > 0.0: # Same sign as enr1
                enr1 = enri
                rho1 = rhoi
                if side == -1:
                    rho0 = (rho0-self.rho_tgt)/2.0 + self.rho_tgt
                side = -1
            if (rho0-self.rho_tgt)*(rhoi-self.rho_tgt) > 0.0: # Same sign as enr0
                enr0 = enri
                rho0 = rhoi
                if side == 1:
                    rho1 = (rho1-self.rho_tgt)/2.0 + self.rho_tgt
                side = 1
#            if abs((rho0-self.rho_tgt)*(rho1-self.rho_tgt)) < self.eps_rho**2:
            if abs(rhoi-self.rho_tgt) < self.eps_rho:
                break   # Reactivities close, done
        self.conv_enr    = enri
        self.conv_rho    = rhoi
        self.conv_rhoerr = rhoierr
        if my_debug > 4:
            print('DONE: ',self)
        return True

    def plot_iters(self, plot_file:str='rf_enr_iter.png'):
        'Plot how the iteration went'
        if not self.rholist:
            print("Warning: No iterations to plot!")
            return
        xvals = [x[0] for x in self.rholist]
        yvals = [x[1] for x in self.rholist]
        yerrs = [x[2] for x in self.rholist]
        fig = plt.figure()
        plt.errorbar(x=xvals, y=yvals, yerr=yerrs)
        plt.title("RF, %s, %5.3f %5.3f" % (self.salt, self.sf, self.l) )
        plt.xlabel("Enrichment")
        plt.ylabel("Reactivity [pcm]")
        plt.grid(True)
        if plot_file == None:
            plt.show()
        else:
            plt.savefig(self.iter_path+'/'+plot_file, bbox_inches='tight')
        plt.close()

    def read_rhos_if_done(self, save_file:str='converge_data.txt') -> bool:
        'Try to load previous search file'
        if os.path.exists(self.iter_path + '/' + save_file) and os.path.getsize(self.iter_path + '/' + save_file) > 50:
            fh = open(self.iter_path + '/' + save_file, 'r')
        else:
            if my_debug:
                print("Results not available in: ", \
                   self.iter_path + '/' + save_file)
            return False
        myline  = fh.readline().strip()
        mysalt  = myline.split()[5]          
        mysf    = float(myline.split()[7])
        myl     = float(myline.split()[9])
        if not (mysalt==self.salt and mysf==self.sf and myl==self.l):
            print("ERROR: Lattice parameters do not match!")
            return False

        for myline in fh.readlines():
            myline  = myline.strip().split()
            myenr   = float(myline[0])
            myrho   = float(myline[1])
            myrhoerr= float(myline[2])
            self.rholist.append( self.RhoData(myenr, myrho, myrhoerr) )


        if len(self.rholist) < 3:
            return False
        found_enr1   = self.rholist[-2][0]
        found_enr0   = self.rholist[-1][0]
        found_rho1   = self.rholist[-2][1]
        found_rho0   = self.rholist[-1][1]
        found_rhoerr = self.rholist[-1][2]
        if abs((found_rho0 - self.rho_tgt)*(found_rho1-self.rho_tgt)) < self.eps_rho**2 or \
            abs(found_enr1 - found_enr0) < self.eps_enr*abs(found_enr0+found_enr1):
            self.conv_enr    = self.rholist[-1][0]
            self.conv_rho    = self.rholist[-1][1]
            self.conv_rhoerr = self.rholist[-1][2]
            if my_debug:
                print("*** READ :",self)
            return True
        else:
            return False

    def save_iters(self, save_file:str='converge_data.txt'):
        'Save history of the iterative search'
        if not self.rholist:
            print("Warning: No iterations to save!")
            return
        result = '# enr, rho, sig_rho for %s sf %5.3f l %5.3f\n' \
                                % (self.salt, self.sf, self.l)
        for r in self.rholist:
            result += '%10.8f\t%10.2f\t%6.1f\n' % (r)

        try:
            fh = open(self.iter_path + '/' + save_file, 'w')
            fh.write(result)
            fh.close()
        except IOError as e:
            print("[ERROR] Unable to write to file: ", \
                   self.iter_path + '/' + save_file)
            print(e)

    def cleanup_force_all(self):
        'Nuke base dir'
        if os.path.isdir(self.iter_path):
            shutil.rmtree(self.iter_path)

    def cleanup(self, preserve_last:bool=False):
        'Delete all run directories'
        n_lats = len(self.rholist) -1
        if preserve_last:
            n_last -= 1         # Save the last run
        for i in range(n_lats):
            l = Lattice(self.salt, self.sf, self.l, rholist[i][1])
            l.cleanup()


# ------------------------------------------------------------
if __name__ == '__main__':
    print("This module find critical enrichment of a lattice.")
#    input("Press Ctrl+C to quit, or enter else to test it.")
    c = Converge('flibe',0.07,11.0)
    c.enr_min = 0.006149271528
    c.enr_max = 0.016040111061
    c.read_rhos_if_done()
    c.iterate_rho()
    c.save_iters()
    print("Enrichment for %s sf %5.3f l %5.3f -> %7.5f" % \
            (c.salt, c.sf, c.l, c.conv_enr) )

