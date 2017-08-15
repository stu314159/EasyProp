# simpleRankine.py 
"""
Example problem using EasyProp for fluid properties.
Problem is adapted from Example 8-5 of Introduction to Thermal and Fluids Engineering,
1st Edition, D. Kaminski and M. Jensen, John Wiley & Sons

A Rankine cycle uses water as the working fluid.  The boiler operates at 8 MPa
and produces saturated vapor.  Saturated liquid exits the condenser at 7.5 kPa.
The pump has an isentropic efficiency of 82 percent, and the turbine has an 
isentropic efficiency of 88 percent.  The steam flow rate is 2.8x10**4 kg/h.
Cooling water enters the condenser at 20 C and leaves at 38 C.

Find:
a) the net power produced (kW)
b) the net heat transfer rate into the boiler (kW)
c) cycle thermal efficiency (%)
d) cooling water flow rate (kg/h)

"""
import numpy as np

import sys
sys.path.insert(0,'..') #<-- assumes EasyProp is one folder up

import EasyProp as EP

myFluid = EP.EasyProp('Water','SI')
numStatePoints = 4;
eta_turbine = 0.88;
eta_pump = 0.82;

# declare fluid property arrays
h = np.zeros((numStatePoints,),dtype=np.float)
h_s = np.zeros_like(h);
s = np.zeros_like(h);
s_s = np.zeros_like(h);
P = np.zeros_like(h);
T = np.zeros_like(h);
T_s = np.zeros_like(h);

# state point 1: condensate pump inlet
P[0] = 7.5; # given pressure at condenser outlet / condensate pump inlet
h[0] = myFluid.hL_p(P[0])
s[0] = myFluid.sL_p(P[0])

