# exampleAirBrayton.py
"""
This example problem is adapted from example 9.8 of Fundamentals of Engineering Thermodynamics,
Moran, Shapiro, et al. 7th edition.

Consider the following Brayton cycle:
Air enters the compressor at 100 kPa, 27 C and is compressed to 1000 kPa.  The temperature at
inlet of the first turbine stage is 1127 C.  The expansion takes place isentropically in two
stages, with reheat to 1127 C between the stages at a constant pressure of 300 kPa.
A regenerator having an effectiveness of 100% is also incorporated into the cycle.

Find:
a) The temperature at each state point
b) Thermal efficiency

"""

import numpy as np

import sys
sys.path.insert(0,'..') #<-- assumes EasyProp is one folder up

import EasyProp as EP

myFluid = EP.EasyProp('Air','SI')
numStatePoints = 8;
regen_eff = 1.0;
eta_turbine = 1.0;
eta_compressor = 1.0;

# declare fluid property arrays
h = np.zeros((numStatePoints,),dtype=np.float)
h_s = np.zeros_like(h);
s = np.zeros_like(h);
s_s = np.zeros_like(h);
P = np.zeros_like(h);
T = np.zeros_like(h);
T_s = np.zeros_like(h);

# state point 0:
T[0] = 27 # C - given
P[0] = 100 # kPa - given
h[0] = myFluid.h_pT(P[0],T[0])
s[0] = myFluid.s_pT(P[0],T[0])

# compression to state point 1:
P[1] = 1000 # kPa - given
s_s[1] = s[0]
h_s[1] = myFluid.h_ps(P[1],s_s[1])
h[1] = h[0] - (h[0] - h_s[1])/eta_compressor
T[1] = myFluid.T_ph(P[1],h[1])

w_comp = h[0] - h[1]

print("Compressor specific work = %g kJ/kg"%w_comp)

# skip to combustor outlet / turbine #1 inlet
P[3] = 1000 #kPa (regeneration and combustion stages are assumed isobaric)
T[3] = 1127 # C given
h[3] = myFluid.h_pT(P[3],T[3])
s[3] = myFluid.s_pT(P[3],T[3])

# expansion through turbine #1
P[4] = 300 #kPa given
s_s[4] = s[3]
h_s[4] = myFluid.h_ps(P[4],s_s[4])
h[4] = h[3] - (h[3] - h_s[4])*eta_turbine
T[4] = myFluid.T_ph(P[4],h[4])

w_turb1 = h[3] - h[4]

# isobaric reheat
T[5] = 1127 # C given
P[5] = P[4]
h[5] = myFluid.h_pT(P[5],T[5])
s[5] = myFluid.s_pT(P[5],T[5])

q_rh = h[5] - h[4]

# expansion through turbine #2
P[6] = 100 # kPa assume regeneration and exhaust are done isobaric processes
s_s[6] = s[5]
h_s[6] = myFluid.h_ps(P[6],s_s[6])
h[6] = h[5] - (h[5] - h_s[6])*eta_turbine
T[6] = myFluid.T_ph(P[6],h[6])

w_turb2 = h[5] - h[6]

# use regen effectiveness to find enthalpy change across regenerator
h[2] = h[1] + regen_eff*(h[6] - h[1])
P[2] = P[1] # regen isobaric on low temp side
T[2] = myFluid.T_ph(P[2],h[2])

# state point 7 - regenerator exhaust
h[7] = h[6] - (h[2] - h[1]) # energy balance assuming equal mass flow rate on both sides of regen.
P[7] = P[6] # again, assume isobaric
T[7] = myFluid.T_ph(P[7],h[7])

# first law balance

q_s = h[3] - h[2]
q_r = h[0] - h[7]

q_net = q_s + q_r + q_rh
w_net = w_comp + w_turb1 + w_turb2
eta_th = w_net/(q_s+q_rh)

# report the results:
for i in range(numStatePoints):
    print("T[%d] = %g C "%(i,T[i]))
 
print("Net specific heat transferred = %g kJ/kg"%q_net)
print("Net specific work = %g kJ/kg"%w_net)
print("Cycle Thermal Efficiency = %5.2f percent"%(eta_th*100.))
    

