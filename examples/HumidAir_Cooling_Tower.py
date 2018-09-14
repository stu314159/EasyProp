import numpy as np

import sys
sys.path.insert(0,'..') #<-- assumes EasyProp is one folder up

import EasyProp as EP

water = EP.simpleFluid('Water','SI')
humidAir = EP.HumidAir('SI')
air = EP.simpleFluid('Air','SI')
numStatePoints = 5;



# declare fluid property arrays
h = np.zeros((numStatePoints,),dtype=np.float) # enthalpy of water, keep these in kJ/kg or BTU/lbm
ha = np.zeros((numStatePoints,),dtype=np.float) # mixture enthalpy of dry air husing HumidAir obejct
ha_d = np.zeros_like(h);
T = np.zeros_like(h);
R = np.zeros_like(h); # relative humidity (for air properties)
W = np.zeros_like(h);

P = 101.325; # Pa, atmospheric pressure (constant)
m_dot_water = 4.5e7; # kg/hr of water in/out


# problem parameters
# State point #1 - warm water in:
T[0] = 38; # C
h[0] = water.h_pT(P,T[0])

# State point #2 - cooled water out:
T[1] = 30; #C
h[1] = water.h_pT(P,T[1])

# State point #3 - atmospheric air in
T[2] = 25; # C
R[2] = 0.35;
h[2] = water.h_pT(P,T[2])
ha[2] = humidAir.h_PTR(P,T[2],R[2]) # figure out the syntax, okay?
ha_d[2] = air.h_pT(P,T[2])
W[2] = humidAir.w_PTR(P,T[2],R[2]); 

# state point #4 - moist air out
T[3] = 35;
R[3] = 0.9;
h[3] = water.h_pT(P,T[3])
ha[3] = humidAir.h_PTR(P,T[3],R[3]);
ha_d[3] = air.h_pT(P,T[3])
W[3] = humidAir.w_PTR(P,T[3],R[3]);

# state point #5 - make-up water in
T[4] = 20;
h[4] = water.h_pT(P,T[4])

print "h = "
print h

print "ha = "
print ha

print "ha_d = "
print ha_d

print "W = "
print W

m_dot_air = m_dot_water*(h[1]-h[0])/(ha[2]-ha[3]+(W[3]-W[2])*h[4]);
print "air mass flow rate = %g kg/hr"%m_dot_air

m_dot_makeup = m_dot_air*(W[3]-W[2])
print "makeup water flow rate = %g kg/hr"%m_dot_makeupf

