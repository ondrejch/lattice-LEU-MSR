#!/usr/bin/python3
#
# Ondrej Chvala, ochvala@utk.edu
# 2019-08-12
# GNU/GPL

import math
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
import time
import molmass
import serpentTools
serpentTools.settings.rc['verbosity']='error'

from salts import Salt

SALTS = {
    "flibe"  : "72%LiF + 16%BeF2 + 12%UF4",
    "lif"    : "73%LiF + 27%UF4",
    "naf"    : "78%NaF + 22%UF4",
    "nafbe12": "74%NaF + 12%BeF2 + 14%UF4",
    "nafbe30": "58%NaF + 30%BeF2 + 12%UF4",
    "nafrbf2": "46%NaF + 33%RbF + 21%UF4",
    "nafzrf" : "49%NaF + 38%ZrF4 + 13%UF4",
    "nafkf"  : "50.5%NaF + 21.5%KF + 28%UF4" }

NUCLEAR_DATA_LIBS = {'ENDF7', 'ENDF8', 'JEFF33' }

do_plots = False
my_debug = False

class Lattice(object):
    'Infinite 2D lattice of graphite hexagonal blocks with a cylindrical fuel channel'
    def __init__(self, salt:str='flibe', sf:float=0.1, l:float=20.0, e:float=0.015):
        try:                            # Check if salt is known
            self.salt_formula = SALTS[salt]
        except:
            ValueError("Salt "+salt+" is undefined.")
        if sf>0.6 or sf<1e-4 or l<0.05 or l>100.0 or e>1.0 or e<0.0: # Reject bad input
            raise ValueError("Bad parameters: ", sf, l, e)

        # Lattice parameters
        self.l:float       = l          # Hex lattice size [cm]
        self.sf:float      = sf         # Fuel salt fraction
        self.r:float       = self.r()   # Diameter of the fuel salt channel [cm]
        self.salt_name:str = salt       # Salt identifier
        self.s             = Salt(self.salt_formula, e) # Salt used
        self.tempK:float   = 900.0      # Salt temperature [K]

        self.k:float       = None       # ANALYTICAL_KEFF
        self.kerr:float    = None       # ANALYTICAL_KEFF Error
        self.cr:float      = None       # CONVERSION_RATIO
        self.crerr:float   = None       # CONVERSION_RATIO Error

        self.nuc_libs:str  = 'ENDF7'    # Nuclear data library
        self.lib:str       = '09c'      # CE xsection temp selection
        self.queue:str     = 'fill'     # NEcluster torque queue
        self.histories:int = 5000       # Neutron histories per cycle
        self.ompcores:int  = 8          # OMP core count
        self.deck_name:str = 'lat'      # Serpent input file name
        self.deck_path:str = '.'        # Where to run the lattice deck
        self.qsub_path:str = os.path.expanduser('~/run.sh')  # Full path to the qsub script
        self.main_path:str = os.path.expanduser('~/L/')+salt # Main path
        self.boron_graphite:float = 2e-06     # 2ppm boron in graphite
        if my_debug:
            print("DEBUG LATTICE ", self.salt_formula, self.sf, self.l, self.s.enr)

    def __repr__(self):
        return "Lattice: "+repr(self.s)+" sf: "+repr(self.sf)+" l: "+repr(self.l)

    def set_path_from_geometry(self):
        'Sets path to directory to run cases based on geometry'
        self.deck_path = self.main_path + "/" + "%08.6f"%self.sf + \
            "/%08.5f"%self.l + "/%014.12f"%self.s.enr

    def hexarea(self) -> float:
        'Area of the lattice [cm2]'
        return 2.0 * math.sqrt(3.0) * self.l**2

    def r(self) -> float:
        'Radius of inner channel [cm]'
        return math.sqrt(self.sf * self.l**2 *math.sqrt(3.0) / (2.0 * math.pi))

    def get_cells(self) -> str:
        'Cell cards for Serpent input deck'
        cells = '''
%______________cell definitions_____________________________________
cell 11  0  fuelsalt   -1      % fuel channel
cell 12  0  graphite    1 -2   % graphite
cell 99  0  outside     2      % graveyard
'''
        return cells.format(**locals())

    def get_surfaces(self) -> str:
        'Surface cards for Serpent input deck'
        surfaces = '''
%______________surface definitions__________________________________
surf 1   cyl    0.0 0.0 {self.r}  % fuel channel radius
surf 2   hexxc  0.0 0.0 {self.l}  % reflective unit cell boundary
'''
        return surfaces.format(**locals())

    def get_graphite(self) -> str:
        'Graphite material cards'
        gr_lib   = self.lib
        gr_frac  = 1.0 - self.boron_graphite
        b10_frac = 0.2 * self.boron_graphite
        b11_frac = 0.8 * self.boron_graphite
        graphite = '''
%  NUCLEAR GRAPHITE: Natural concentration of carbon
%  DENSITY: 1.80 G/CC
mat graphite -1.80 moder graph 6000 tms 950
rgb 130 130 130
6000.{gr_lib} {gr_frac}
5010.{gr_lib} {b10_frac} % boron impirity eq.
5011.{gr_lib} {b11_frac}
%  THERMAL SCATTERING LIBRARY FOR GRAPHITE
therm graph 0 gre7.18t gre7.20t gre7.22t
'''
        return graphite.format(**locals())

    def get_data_cards(self) -> str:
        'Data cards for the lattice'
        data_cards = '''
% Boundary condition
set bc 2

% Neutron population and criticality cycles
set pop {self.histories} 100 40 % {self.histories} neutrons, 100 active, 40 inactive cycles

% Analog reaction rate
% set arr 2
'''
        if self.nuc_libs == "ENDF7":
            data_cards += '''
% Data Libraries
set acelib "sss_endfb7u.sssdir"
set declib "sss_endfb7.dec"
set nfylib "sss_endfb7.nfy"
'''
        if do_plots:
            data_cards += '''
% Plots
plot 3 1500 1500
% mesh 3 1500 1500
'''
        return data_cards.format(**locals())

    def get_deck(self) -> str:
        'Serpent deck for the lattice'
        deck = '''\
set title "MSR lattice cell, l {self.l}, sf {self.sf}, salt {self.salt_formula}, Uenr {self.s.enr} "
'''
        deck += self.get_surfaces()
        deck += self.get_cells()
        deck += self.s.serpent_mat(self.tempK)
        deck += self.get_graphite()
        deck += self.get_data_cards()
        return deck.format(**locals())

    def save_deck(self):
        'Saves Serpent deck into an input file'
        try:
            os.makedirs(self.deck_path)
            fh = open(self.deck_path + '/' + self.deck_name, 'w')
            fh.write(self.get_deck())
            fh.close()
        except IOError as e:
            print("[ERROR] Unable to write to file: ", \
                   self.deck_path + '/' + self.deck_name)
            print(e)

    def save_qsub_file(self):
        'Writes run file for TORQUE.'
        qsub_content ='''#!/bin/bash
#PBS -V
#PBS -N Serp2MSR_lat
#PBS -q {self.queue}
#PBS -l nodes=1:ppn={self.ompcores}

hostname
rm -f done.dat
cd ${{PBS_O_WORKDIR}}
module load mpi
module load serpent

sss2 -omp {self.ompcores} {self.deck_name} > myout.out
awk 'BEGIN{{ORS="\\t"}} /ANA_KEFF/ || /CONVERSION/ {{print $7" "$8;}}' {self.deck_name}_res.m > done.out
rm {self.deck_name}.out
'''.format(**locals())
        try:                # Write the deck
            f = open(self.qsub_path, 'w')
            f.write(qsub_content)
            f.close()
        except IOError as e:
            print("Unable to write to file", fname)
            print(e)

    def run_deck(self):
        'Runs the deck using qsub_path script'
        if self.queue is 'local':    # Run the deck locally
            os.chdir(self.deck_path)
            os.system(self.qsub_path)
        else:               # Submit the job on the cluster
            os.system('cd '+ self.deck_path + ' && qsub ' + self.qsub_path)

    def get_calculated_values(self) -> bool:
        'Fill k and cr for lattice if calculated'
        if os.path.exists(self.deck_path+'/done.out') and os.path.getsize(self.deck_path+'/done.out') > 30:
            pass
        else:                   # Calculation not done yet
            return False

        results = serpentTools.read(self.deck_path + '/' + self.deck_name + "_res.m")
        self.k      = results.resdata["anaKeff"][0]
        self.kerr   = results.resdata["anaKeff"][1]
        self.cr     = results.resdata["conversionRatio"][0]
        self.crerr  = results.resdata["conversionRatio"][1]

        if my_debug:
            print("[DEBUG Lat] ---> k = {self.k}, CR = {self.cr}".format(**locals()))
        return True

    def cleanup(self):
        'Delete the run directory'
        if os.path.isdir(self.deck_path):
            shutil.rmtree(self.deck_path)

# ------------------------------------------------------------
if __name__ == '__main__':
    print("This module handles a simple lattice.")
    input("Press Ctrl+C to quit, or enter else to test it.")
    l = Lattice()
    print("***** Serpent deck: \n" +  l.get_deck() + "\n***** ")
#    l.ompcores=32
    l.deck_path = os.path.expanduser('~/tmp/lat_test')
    l.save_qsub_file()
    l.save_deck()
    l.run_deck()
    while not l.get_calculated_values():
        print("Wating for Serpent ...")
        time.sleep(5.0)
    print(l.k, l.cr)



