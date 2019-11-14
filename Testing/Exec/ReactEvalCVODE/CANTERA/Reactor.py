###############################################################
#
# Cantera UPF : freely-propagating premixed flames
# Description : perform calculations of UPF for a range of eq. ratio
#               and compute global indices (integrated fuel, CO, NO consumption, ...)
# Author : A. Felden
# Last modified : 02/2018
#
###############################################################

#import :
from cantera import *
import numpy as np
import csv
import sys,os

#################################################################
# Prepare your run
#################################################################
#Mechanism used for the process
cas         = 'dodecane_reduced_mech_jetsurf.cti' 
#dodecane_reduced_mech_jetsurf.cti' 

gas         = Solution(cas,'gas')

########
#STORAGE
########

#################
#Run parameters:

#################
#Stoechiometry :

m            = gas.n_species

fuel_species = 'NC12H26'
#phi          = 1.0
#stoich_O2    = gas.n_atoms(fuel_species,'C') + 0.25*gas.n_atoms(fuel_species,'H')
#stoich_O2    =  0.25*gas.n_atoms(fuel_species,'H')
#air_N2_O2_ratio = 3.76

ifuel = gas.species_index(fuel_species)
io2 = gas.species_index('O2')
in2 = gas.species_index('N2')

x        = np.zeros(m,'d')
x[ifuel] = 0.1 #phi                        # 4.34783e-002
x[io2]   = 0.2 #stoich_O2                  # 2.00866e-001
x[in2]   = 0.7 #stoich_O2*air_N2_O2_ratio  # 7.55656e-001

#General parameter values :
tin        	= 2000.0        # unburned gas temperature
p_min 	        = 101325.0

#Gas def
gas.TPX         = tin, p_min, x

csv_file = "datafromSC.dat_"+str(cas)+"-"+str(fuel_species)
with open(csv_file, 'w') as outfile:
    writer = csv.writer(outfile,delimiter=' ',quotechar=' ', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(['1'])
    writer.writerow([gas.P,gas.T, 0.0,]+list(gas.X))


#################
#Simulation

#################
#r = IdealGasConstPressureReactor(gas)
r = IdealGasReactor(gas)
sim = ReactorNet([r])
time_init = 5.0e-05
nbdt = 500
dt_react_incr = time_init / nbdt
#print("dt rect ?", dt_react_incr)
times = np.zeros(nbdt)
data = np.zeros((nbdt,4))

print('#Time', 'temp', 'e', 'h', 'Press', 'rho')
time = 1.0e-10
time = 0.0
for n in range(nbdt):
    time += dt_react_incr 
    sim.advance(time)
    times[n] = time 
    data[n,0] = r.T
    data[n,1:] = r.thermo[fuel_species,'O2','H2'].Y * r.density
    #print('%10.3e %10.3f %10.3f %14.6e ')% (sim.time, r.T, r.thermo.P, r.thermo.u)
    print(sim.time, r.T)#, r.thermo.u, r.thermo.h, r.thermo.P, r.density)
#print(r.thermo.Y)
#print(data[nbdt-1,1:])
