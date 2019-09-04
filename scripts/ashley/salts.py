#!/usr/bin/python3
#
# Module that handles molten salt properties for neutronic simulations
#
# Ondrej Chvala, ochvala@utk.edu
# 2019-08-06
# GNU/GPL

import math
from collections import namedtuple
import molmass          # https://pypi.org/project/molmass/
import copy

my_debug = False

MOLARVOLUMES = { # Melt composition molar volumes at 600 and 800 degC
    'LiF' : (13.411, 14.142),
    'BeF2': (23.6, 24.4),
    'NaF' : (19.08, 20.2),
    'KF'  : (28.1, 30.0),
    'ZrF4': (47.0,  50.0),
    'ThF4': (46.43, 47.59),
    'UF4' : (46.43, 47.59)}

class MeltPart(object):
    'Storage for salt density fit calculation'
    def __init__(self, f:str, molf:float, enr:float):
        try:
            self.molar_vols = MOLARVOLUMES[f]
        except:
            raise ValueError("Molar volumes of "+f+" undefined!")
        self.formula:str      = f
        self.molar_frac:float = molf
        self.s = Salt("100%"+f, enr)
    def __repr__(self):
        return "%s, %s" % (repr(self.formula), repr(self.s))

class IsoWeightFraction(object):
    '''Class for salts isotopic weight fractions.
       Ntuples are immutable in Python, use a class instead'''
    def __init__(self, Z:int, A:int, wf:float):
        self.Z:int      = Z
        self.A:int      = A
        self.wf:float   = wf
    def __repr__(self):
        return "%2i %3i  %10.8f" % (self.Z, self.A, self.wf)

class Salt(object):
    'Class for salt parsing, based on salt formula and enrichment'
    def __init__(self, f:str="72%LiF + 16%BeF2 + 12%UF4", e:float=0.02):
        'Constructor using salt formula, uranium enrichment, and salt density'
        try:
            f = f.strip().replace(" ", "")
        except:
            raise ValueError("Formula " + self.formula + " error")

        self.formula:str    = f         # Chemical formula for a salt
        self.enr:float      = e         # Uranium enrichment
        self.Li7dep:float   = 0.99995   # Li-7 depletion level
        self.mol_mass:foat  = None      # Molar mass of the salt
        # Salt isotopic composition - isotopes repeat per melt parts
        self.isolist = []   # For internal processing use only
        self.SaltIso = namedtuple("SaltIso", "Z A atoms amass wfrac molefract")
        # Salt isotopic weight fractions, each isotope is unique
        self.wflist = []

        # Update database if isotopes for our MSR enrichments
        self.ELEMENTS = copy.deepcopy(molmass.ELEMENTS) # Local copy to allow different enrichments
        self.ELEMENTS['Li'].isotopes[6].abundance  = 1.0 - self.Li7dep
        self.ELEMENTS['Li'].isotopes[7].abundance  = self.Li7dep
        wf_u234:float = 0.0089 * self.enr
        wf_u236:float = 0.0046 * self.enr
        wf_u238:float = 1.0 - (wf_u234 + self.enr + wf_u236)

        self.ELEMENTS['U'].isotopes[234].abundance = wf_u234
        self.ELEMENTS['U'].isotopes[235].abundance = self.enr
        self.ELEMENTS['U'].isotopes[236]=molmass.Isotope(236.0455611, wf_u236, 236) # Add to dbase
        self.ELEMENTS['U'].isotopes[238].abundance = wf_u238

        # Density calculation
        self.melt_parts = []        # List of , enr:floatMeltPart objects
        self.density_a:float = None # Linear density interpolation slope
        self.density_b:float = None # Intercept

        if my_debug:
            print(self)

    def __repr__(self):
        result = "Salt: %s, Uenr= %f " % (self.formula, 100.0*self.enr) + "%"
        if self.isolist:
            for i in self.isolist:
                result += "\n"+repr(i)
        if self.mol_mass:
            result += "\nMolar mass %f g/mole" % (self.mol_mass)
        if self.wflist:
            result += "\nIsotopic Weight fractions:"
            twf = 0.0
            for w in self.wflist:
                twf += w.wf
                result += "\n"+repr(w)
            result += "\n---Sum: %10.8f" % twf
        return result

    def _formula_parse_iso(self):
        'Parse chemical formula of the salt and get list of all isotopes'
        tot_moles:float = 0.0                   # Total molar fraction, should add to 1
        for meltpart in self.formula.split('+'):# Separate melt components
            mfract, comp = meltpart.split('%')  # Separate component pct. fractions
            mfract = float(mfract)/100.0        # Molar % -> fraction
            tot_moles += mfract                 # Add molar fractions of compositions
            comp_f = molmass.Formula(comp)      # Turn component into a molmass formula
            for symbol in comp_f._elements:     # Elements in a component
                ele = self.ELEMENTS[symbol]  # Get the element object
                for m, n_atoms in comp_f._elements[symbol].items(): # Number of atoms in a component
                    if(m != 0):                 # This should be zero per molmass
                        raise Exception("Error in parsing")
                    for A in ele.isotopes.keys():   # Get data for each isotope
                        Z     = ele.protons
                        amass = ele.isotopes[A].mass
                        wfrac = ele.isotopes[A].abundance
                        if wfrac > 0.0:         
                            isotuple = self.SaltIso(Z, A, n_atoms, amass, wfrac, mfract)
                            self.isolist.append(isotuple)
        self.isolist.sort()                     # Looks nicer sorted
        if abs(tot_moles - 1.0) > 1e-5:         # Sanity check
            raise ValueError("User Error: Formula "+self.formula+" molar fractions do not add to 1.0!")

    def _molar_mass(self):
        'Establish molar mass of the salt'
        if not self.isolist:        # Generate list of isotopes
            self._formula_parse_iso()
        self.mol_mass = 0.0
        for i in self.isolist:      # Add molar weights from all isotopes
            self.mol_mass += i.molefract * i.atoms * i.amass * i.wfrac

    def get_molar_mass(self) ->float:
        'Returns molar weight [g/mole]'
        if not self.mol_mass:   # Establish molar mass of the salt
            self._molar_mass()
        return self.mol_mass

    def _isotopic_fractions(self):
        'Establish isotopic fractions'
        if not self.isolist:    # Generate list of isotopes
            self._formula_parse_iso()
        if not self.mol_mass:   # Establish molar mass of the salt
            self._molar_mass()
        for i in self.isolist:  # Process all isotopes in the isolist
            w_l = [x for x in self.wflist if x.Z==i.Z and x.A==i.A]
            if not w_l:         # Isotope not in self.wflist, add new one
                iwf = IsoWeightFraction(i.Z, i.A, i.molefract * i.atoms * i.amass * i.wfrac)
                self.wflist.append(iwf)
            else:               # Isotope is in self.wflist, add to the mass
                if len(w_l) > 1:    # We should get list of length 1, each isotope is unique
                    raise Exception("Error: there should only be one ", w_l)
                w = w_l[0]
                w.wf += i.molefract * i.atoms * i.amass * i.wfrac
        twf = 0.0               # Total weight fraction, should add to 1
        for w in self.wflist:   # Normalize each isotope by molar mass of the salt
            w.wf /= self.mol_mass   # Fraction of the mass of that isotope
            twf += w.wf         # Add each weight fraction
        if abs(twf - 1.0) > 1e-12:          # Sanity check
            raise ValueError("Error: weight fractions do not add to 1.0!")

    def _fit_density(self):
        'Uses molar counting method to get density fit coefficients'
        for m in self.formula.split('+'):   # Separate melt components
            mfract, comp = m.split('%')     # Separate component pct. fractions
            mfract = float(mfract)/100.0    # Molar % -> fraction
            self.melt_parts.append( MeltPart(comp, mfract, self.enr) )
        weight_600C = 0.0
        weight_800C = 0.0
        volume_600C = 0.0
        volume_800C = 0.0
        if self.melt_parts == 1:    # Single component melt
            self.melt_parts[0].s.mol_mass = self.get_molar_mass()
        for mp in self.melt_parts:
            weight_600C += mp.molar_frac*mp.s.get_molar_mass()
            weight_800C += mp.molar_frac*mp.s.get_molar_mass()
            volume_600C += mp.molar_frac*MOLARVOLUMES[mp.formula][0]
            volume_800C += mp.molar_frac*MOLARVOLUMES[mp.formula][1]
        density_600C = weight_600C / volume_600C
        density_800C = weight_800C / volume_800C
        self.density_a = (density_800C - density_600C) / (800.0 - 600.0)
        self.density_b = density_600C - self.density_a*600.0
        if my_debug:
            print("  Density at 600 and 800C:", density_600C, density_800C)
            print("  Fit a, b:", self.density_a, self.density_b)

    def densityK(self, tempK:float) -> float:
        'Returns density [g/cm3] based on temperature in Kelvin'
        return self.densityC(tempK - 273.15)

    def densityC(self, tempC:float) -> float:
        'Returns density [g/cm3] based on temperature in degC'
        if tempC < 600 or tempC > 800:
            print("Warning: temperature data is interpolated between 600 and 800C.")
        if not self.density_a or not self.density_b:
            self._fit_density()     # Necessary to prevent infinite recursion..
        return self.density_a * tempC + self.density_b

    def nice_name(self)->str:
        'Return salt name with spaces around + sign'
        return self.formula.replace('+',' + ')        

    def serpent_mat(self, tempK:float=900, lib="09c", rgb:str="240 30 30"):
        'Returns Serpent deck for the salt material'
        if not self.wflist:         # Generate list of isotopic weight fractions
            self._isotopic_fractions()
        if my_debug:                # Check uranium enrichment 
            u= 0.0
            for w in self.wflist: 
                if w.Z == 92: 
                    u += w.wf
            for w in self.wflist:
                if w.Z == 92:
                    print("DEBUG SALT: %d -> %8.3f" % (w.A, 100.0*w.wf/u) )
        mat  = "% Fuel salt: " + self.nice_name() + ", U enrichment " + str(self.enr)
        mat += "\nmat fuelsalt %12.8f rgb %s tmp %8.3f\n" % (-1.0*self.densityK(tempK),rgb,tempK)
        for w in self.wflist:
            mat += "%3d%03d.%s  %14.12f" % (w.Z, w.A, lib, -1.0*w.wf)
            mat += "    %  "+ self.ELEMENTS[w.Z].symbol +"-"+ str(w.A) +"\n"
        return mat

# This executes if someone tries to run the module
if __name__ == '__main__':
    print("This is a salt processing module.")
    input("Press Ctrl+C to quit, or enter else to test it. ")
    s = Salt()
    print(s)
    print("\n\n--> Serpent deck:\n",s.serpent_mat())
    print("--> Density [g/cm3] at 700C: ",s.densityC(700))
    print("--> Density [g/cm3] at 900K: ",s.densityK(900))

