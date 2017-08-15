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
m_dot_steam = 2.8e4; # kg/h

# declare fluid property arrays
h = np.zeros((numStatePoints,),dtype=np.float)
h_s = np.zeros_like(h);
s = np.zeros_like(h);
s_s = np.zeros_like(h);
P = np.zeros_like(h);
T = np.zeros_like(h);
T_s = np.zeros_like(h);

# state point 0: condensate pump inlet
P[0] = 7.5; # given pressure at condenser outlet / condensate pump inlet
h[0] = myFluid.hL_p(P[0])
s[0] = myFluid.sL_p(P[0])

# compression from state point 0 to 1
P[1] = 8000.; # boiler pressure of 8 MPa given
s_s[1] = s[0]
h_s[1] = myFluid.h_ps(P[1],s_s[1]) # isentropic compression
h[1] = h[0] - (h[0] - h_s[1])/eta_pump; #<-- correct for pump efficiency
w_pump = h[0] - h[1] # real pump specific work (kJ/kg)

# isobaric heat addition in the boiler from state point 1 to 2
P[2] = P[1]
h[2] = myFluid.hV_p(P[2]) #<-- saturated vapor leaves the boiler (given)
s[2] = myFluid.sV_p(P[2])

q_s = h[2] - h[1] #<-- specific heat supplied (kJ/kg)

# expansion from state point 2 to 3 in the turbine
s_s[3] = s[2]
P[3] = P[0]
h_s[3] = myFluid.h_ps(P[3],s_s[3])
h[3] = h[2] - (h[2] - h_s[3])*eta_turbine
w_turbine = h[2] - h[3] # real turbine specific work (kJ/kg)

# isobaric heat rejection in the condenser from state point 3 to 0
q_r = h[0] - h[3] #<-- specific heat rejected (kJ/kg)

w_net = w_pump + w_turbine;
q_net = q_s + q_r;

# check energy balance; w_net should equal q_net
print "Net specific work = %g kJ/kg "%w_net
print "Net specific heat = %g kJ/kg "%q_net


# answer some questions:

#a) net power produced = w_net*m_dot_stm
net_power = w_net*m_dot_steam*(1./3600.) #(kW = kJ/s)

print "Net power = %g kW"%net_power

#b) net heat transfer rate into the boiler
net_Qdot_supplied = q_s*m_dot_steam*(1./3600.) #kW =  kJ/s

print "Net heat transfer rate into the boiler = %g kW"%net_Qdot_supplied

#c) cycle thermal efficiency
eta_th = w_net/q_s 
print "Thermal Efficiency = %4.1f percent "%(eta_th*100.)

#d) cooling water flow rate.

# state point 4 = condenser cooling water inlet condition
Po = 101.325 # kPa assume standard atmospheric pressure for condenser cooling water
Tcond_in = 20. #C - given
Tcond_out = 38. #C - given
h_cond_in = myFluid.h_pT(Po,Tcond_in)
h_cond_out = myFluid.h_pT(Po,Tcond_out)

# energy balance on condenser for cooling water flow rate
m_dot_cool = m_dot_steam*(h[3]-h[0])/(h_cond_out - h_cond_in); #kg/h

print "Condenser cooling water flow rate = %g kg/h"%m_dot_cool



