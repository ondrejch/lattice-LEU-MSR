#!/usr/bin/python3
#
# Class to find reactivity feedback coefficients of a lattice
#
# Ondrej Chvala, ochvala@utk.edu
# 2020-03-20
# GNU/GPL

import time
import numpy as np
from lattice import Lattice
from converge import rho

my_debug:int = 1

GRAPHITE_CTE:float = 3.5e-6   # Graphite linear thermal expansion coefficient (CTE) [m/m per K]
GRAPHITE_RHO:float = 1.80     # Graphite density at 950 K [g/cm3]
FB_BASE_TEMP:float = 900.0    # Base temperature for feedbacks [K]
FB_TEMPS      = [800.0, 850.0, 900.0, 950.0, 1000.0]  # Temps [K] for feedback calculations
SLEEP_SEC:int = 30        # Sleep timer between results read attempts [s]
FEEDBACKS = ['fs.dopp', 'fs.void', 'fs.both', 'gr.dopp']

def graphite_linear_expansion(l0:float=100, tempK:float=950.0):
    'Return new lenght based on graphite thermal expansion'
    l_f = l_0 * (1.0 + GRAPHITE_CTE * (tempK - 950.0))
    return l_f

def graphite_density_expansion(tempK:float=950.0):
    'Return new density based on graphite thermal expansion'
    unit_f  = (1.0 + GRAPHITE_CTE * (tempK - 950.0))
    rho_f   = GRAPHITE_RHO / unit_f**3
    return rho_f


class Feedbacks(object):
    '''Calcualte reactivity coefficients of a lattice'''
    def __init__(self, salt:str='flibe', sf:float=0.1, l:float=20.0, e:float=0.015):
        'Constructor with default values'
        self.salt:str  = salt      # Salt key
        self.sf:float  = sf        # Fuel salt fraction
        self.l:float   = l         # Hex lattice size [cm]
        self.e:float   = e         # Fuel salt enrichment
        self.fb_lats   = {}        # Lattice objects for feedback calculations
        self.fb_rhos   = {}        # Reactivities
#        self.fb_rhos_err = {}        # Reactivity sigmas - irrelevant since we run the same statistiscs
        self.alpha     = {}        # Temperature-reactivity feedbacks
        self.force_recalc:bool = False  # Force recalculation of existing data

    def run_feedback(self, feedback:str="fs.dopp"):
        '''Run all salt feedback cases:
            fs.dopp - Fuel salt Doppler
            fs.void - Fuel salt Void
            fs.both - Fuel salt Doppler + void
            gr.dopp - Graphite Doppler
            fs.gr
        '''
        for t in FB_TEMPS:
            fb_lat_name = feedback + "."  + str("%04.0f" % t)
            self.fb_lats[fb_lat_name] = Lattice(self.salt, self.sf, self.l, self.e)
            mylat           = self.fb_lats[fb_lat_name]   # Shorthand
            if   feedback == "fs.dopp":      # Salt Doppler
                mylat.tempK     = FB_BASE_TEMP
                mylat.mat_tempK = t
            elif feedback == "fs.void":      # Salt Void
                mylat.tempK     = t
                mylat.mat_tempK = FB_BASE_TEMP
            elif feedback == "fs.both":      # Salt Void + Doppler
                mylat.tempK     = t
                mylat.mat_tempK = t
            elif feedback == "gr.dopp":
                mylat.tempK     = FB_BASE_TEMP
                mylat.mat_tempK = FB_BASE_TEMP
                mylat.grtempK   = t + 50.0
            else:
                raise ValueError("Feedback " + feedback + " not implemented!")
            if t == FB_BASE_TEMP:   # Base case
                mylat.set_path_from_geometry()
            else:                   # Feedback cases
                mylat.set_path_from_geometry(fb_lat_name)
            if mylat.mat_tempK < 900.0:     # TODO this should be fixed if we generalize
                mylat.lib = '06c'           # nuclear data libraries
            if mylat.grtempK < 900.0:
                mylat.grtempK = '06c'
            if my_debug:
                print(mylat.deck_path)
            if self.force_recalc or not mylat.get_calculated_values():
                mylat.cleanup()
                mylat.save_deck()
                mylat.run_deck()

    def read_feedback(self, feedback:str="fs.dopp"):
        '''Read feedback data and fit reactivity coefficients'''
        while True:     # Wait for all cases to finish
            is_done = True
            for t in FB_TEMPS:
                fb_lat_name = feedback + "."  + str("%04.0f" % t)
                if not self.fb_lats[fb_lat_name].get_calculated_values():
                    is_done = False
            if is_done:     # All done
                break
            if my_debug:
                print("[DEBUG RF] sleeping ...")
            time.sleep(SLEEP_SEC)           # Wait a minute for Serpent ...

        self.fb_rhos[feedback] = []         # Add reactivities as a list
#        self.fb_rhos_err[feedback] = []
        for t in FB_TEMPS:
            fb_lat_name = feedback + "."  + str("%04.0f" % t)
            self.fb_rhos[feedback].append(rho(self.fb_lats[fb_lat_name].k))
#            self.fb_rhos_err[feedback].append(self.fb_lats[fb_lat_name].kerr) # We ignore sigmas, same statistics for all cases
#        w=[1./abs(r*e) for r, e in zip(f.fb_rhos[feedback],f.fb_rhos_err[feedback])]
        (slope, intercept) = np.polyfit(FB_TEMPS, f.fb_rhos[feedback], 1)
        self.alpha[feedback] = slope


# ------------------------------------------------------------
if __name__ == '__main__':
    print("This module finds lattice feedbacks.")
#    input("Press Ctrl+C to quit, or enter else to test it.")
    f = Feedbacks()
    for fb in FEEDBACKS:
        f.run_feedback(fb)
    for fb in FEEDBACKS:
        f.read_feedback(fb)

