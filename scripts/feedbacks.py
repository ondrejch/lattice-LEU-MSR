#!/usr/bin/python3
#
# Class to find reactivity feedback coefficients of a lattice
#
# Ondrej Chvala, ochvala@utk.edu
# 2020-03-20
# GNU/GPL

import time
from lattice import Lattice
from converge import rho

my_debug:int = 0

GRAPHITE_CTE:float = 3.5e-6   # Graphite linear thermal expansion coefficient (CTE) [m/m per K]
FB_BASE_TEMP:float = 900.0    # Base temperature for feedbacks [K]
FB_TEMPS      = [800.0, 850.0, 900.0, 950.0, 1000.0]  # Temperatures [K] for feedback calculations
SLEEP_SEC:int = 30        # Sleep timer between results read attempts [s]


class Feedbacks(object):
    '''Calcualte reactivity coefficients of a lattice'''
    def __init__(self, salt:str='flibe', sf:float=0.1, l:float=20.0, e:float=0.015):
        'Constructor with default values'
        self.salt:str       = salt      # Salt key
        self.sf:float       = sf        # Fuel salt fraction
        self.l:float        = l         # Hex lattice size [cm]
        self.e:float        = e         # Fuel salt enrichment
        self.fb_lats        = {}        # Lattice objects for feedback calculations
#        self.fs_dopp        = None
#        self.fs
        self.force_recalc:bool = False  # Force recalculation of existing data

    def calculate_salt_doppler(self):
        '''Run all salt doppler cases'''
        for t in FB_TEMPS:
            print(t)
            fb_lat_name = "fs.dopp." + str(t)
            if my_debug:
                print(fb_lat_name)
            self.fb_lats[fb_lat_name] = Lattice(self.salt, self.sf, self.l, self.e)
            self.fb_lats[fb_lat_name].tempK = FB_BASE_TEMP
            self.fb_lats[fb_lat_name].mat_tempK = t
            if t == FB_BASE_TEMP:
                self.fb_lats[fb_lat_name].set_path_from_geometry()
            else:
                self.fb_lats[fb_lat_name].set_path_from_geometry(fb_lat_name)
            print(self.fb_lats[fb_lat_name].deck_path)
            if self.force_recalc or not self.fb_lats[fb_lat_name].get_calculated_values():
                self.fb_lats[fb_lat_name].cleanup()
                self.fb_lats[fb_lat_name].save_deck()
                self.fb_lats[fb_lat_name].run_deck()
                while not self.fb_lats[fb_lat_name].get_calculated_values():
                    if my_debug:
                        print("[DEBUG RF] sleeping ...")
                    time.sleep(SLEEP_SEC)  # Wait a minute for Serpent ...



# ------------------------------------------------------------
if __name__ == '__main__':
    print("This module finds lattice feedbacks.")
#    input("Press Ctrl+C to quit, or enter else to test it.")
    f = Feedbacks()
    f.calculate_salt_doppler()
