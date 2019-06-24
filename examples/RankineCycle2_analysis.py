# RankineCycle2_analysis.py
"""
Python-based analysis of Rankine Cycle2.

State point conditions:

state point 1: P = 1.5 psia, x = 0.0
state point 2: P = 164 psia, MCP efficiency = 84 percent
state point 3: P = 164 psia, x = 0.0
state point 4: P = 838 psia, MFP efficiency = 89 percent
state point 5: P = 838 psia, x = 1.0
state point 6: P = 164 psia, HP Turbine efficiency = 94%
state point 7: P = 164 psia, x = 1.0
state point 8: P = 164 psia, T = 490 F
state point 9: P = 1.5 psia, LP Turbine efficiency = 92%
state point 10: P = 838 psia, x = 0.0
state point 11: P = 164 psia, isenthalpic expansion in trap T1
state point 12: P = 164 psia, x = 0.0, liquid drain from moisture
separator
state point 13: P = 2190 psia, T = 610 F
state point 14: P = 2170 psia, T = 535 F
state point 15: P = 14.7 psia, T = 91 F
state point 16: P = 14.7 psia, T = 114.9 F

Calculate the following:

The enthalpy at all state points and the flow fraction f1 and f2.
The net specific work of all turbines and pumps (BTU/lbm)
The net specific heat supplied at the steam generator and condenser
(BTU/lbm)
Cycle thermal efficiency
Flow rate of steam from the Steam Generator and the rate of cooling
water flow through the condenser.
 

"""

import numpy as np
import scipy.optimize as opt

import sys
sys.path.insert(0,'..') #<-- assumes EasyProp is one folder up

import EasyProp as EP

myFluid = EP.simpleFluid('Water','USCS')

numStatePoints = 16+1;#<-- so state point arrays can be treated as 1-based.
# need to decide how to treat the state point numbering hackery.

# some given parameters
eta_mcp = 0.84 # main condensate pump isentropic efficiency
eta_mfp = 0.89 # main feed pump isentropic efficiency
eta_hpt = 0.94 # high pressure turbine isentropic efficiency
eta_lpt = 0.92 # low pressure turbine isentropic efficiency

m_dot_p = 113.5e6 # lbm/hr, primary coolant system flow rate

# declare fluid property arrays
h = np.zeros((numStatePoints,),dtype=np.float)
h_s = np.zeros_like(h);
x = np.zeros_like(h);
s = np.zeros_like(h);
s_s = np.zeros_like(h);
P = np.zeros_like(h);
T = np.zeros_like(h);
T_s = np.zeros_like(h);

# state point 1
P[1] = 1.5 # psia given
h[1] = myFluid.hL_p(P[1])
s[1] = myFluid.sL_p(P[1])

# from sp 1 -> 2: compression in MCP

s_s[2] = s[1]
P[2] = 164 # psi, given
h_s[2] = myFluid.h_ps(P[2],s_s[2])
h[2] = h[1] - (h[1] - h_s[2])/eta_mcp

# state point 3
P[3] = P[2] # isobaric mixing in the OFWH
h[3] = myFluid.hL_p(P[3])
s[3] = myFluid.sL_p(P[3])

# from sp3 -> 4: compression in MFP
P[4] = 838 # psia, given
s_s[4] = s[3]
h_s[4] = myFluid.h_ps(P[4],s_s[4])
h[4] = h[3] - (h[3] - h_s[4])/eta_mfp

# from sp 4 -> 5: isobaric heat addition in the SG
P[5] = P[4]
h[5] = myFluid.hV_p(P[5])
s[5] = myFluid.sV_p(P[5])

# from sp 5 -> 6: expansion in the HP turbine
P[6] = 164. # psia, given
s_s[6] = s[5]
h_s[6] = myFluid.h_ps(P[6],s_s[6])
h[6] = h[5] - eta_hpt*(h[5] - h_s[6])
x[6] = myFluid.x_ph(P[6],h[6]) # this is needed for the energy balance

# from sp 6 -> 7 isobaric, moisture separation
P[7] = P[6]
h[7] = myFluid.hV_p(P[7])
s[7] = myFluid.sV_p(P[7])

# from sp 7 -> 8: isobaric heat addition in the re-heater
P[8] = P[7]
T[8] = 490 #F, given
h[8] = myFluid.h_pT(P[8],T[8])
s[8] = myFluid.s_pT(P[8],T[8])

# from sp 8 -> 9: expansion in the LP turbine
P[9] = 1.5 # psia, given
s_s[9] = s[8]
h_s[9] = myFluid.h_ps(P[9],s_s[9])
h[9] = h[8] - eta_lpt*(h[8] - h_s[9])


# from sp 9 -> 1: isobaric heat rejection in the condenser

# state point 10 high pressure outlet from the Re-heater
P[10] = P[5]
h[10] = myFluid.hL_p(P[10])

# state point 10 -> 11: isenthalpic expansion in trap 1
P[11] = P[2]
h[11] = h[10]

# state point 12: sat liquid drain from moisture separator
P[12] = P[6]
h[12] = myFluid.hL_p(P[12])

# flow fractions f1 and f2 will be determined
# using energy balance relations on the reheater
# and the open feedwater heater

def reheater_balance(f):
    bal = f[0]*h[5] + \
    (1-f[0])*x[6]*(1-f[1])*h[7] \
    -((1-f[0])*x[6]*(1-f[1])*h[8]+f[0]*h[10]);
    
    return bal

def ofwh_balance(f):
    ff2 = (1-f[0])*x[6]*(1-f[1])
    ff11 = f[0]
    ff12 = (1-f[0])*(1-x[6])
    ff7 = (1-f[0])*x[6]*f[1]
    e_in = ff2*h[2]+ff11*h[11]+ff7*h[7]+ff12*h[12]
    e_out = h[3]
    bal = e_in - e_out
    return bal

def wrapped_function(f):
    return [np.abs(reheater_balance(f)), np.abs(ofwh_balance(f))]

f =opt.broyden1(wrapped_function,[0.1,0.1])

# f now holds a vector of flow fractions.
# f[0] corresponds to f1 on the given schematic
# f[1] corresponds to f2 on the given schematic

print("Flow fractions [f1,f2] = ",f)

# calculate specific work of the turbines and pumps

w_mcp = (1-f[0])*x[6]*(1-f[1])*(h[1] - h[2])
w_mfp = h[3] - h[4]
w_hpt = (1-f[0])*(h[5] - h[6])
w_lpt = (1-f[0])*x[6]*(1-f[1])*(h[8] - h[9])

w_net = w_mcp + w_mfp + w_hpt + w_lpt

# report results
print("MCP work = %g BTU/lbm"%w_mcp)
print("MFP work = %g BTU/lbm"%w_mfp)
print("HP Turbine work = %g BTU/lbm"%w_hpt)
print("LP Turbine work = %g BTU/lbm"%w_lpt)
print("Net work = %g BTU/lbm" %w_net)

# net heat in the SG and condenser
q_sg = h[5] - h[4]
q_cond = (1-f[0])*x[6]*(1-f[1])*(h[1]-h[9])

q_net = q_sg + q_cond

print("q in through SG = %g BTU/lbm"%q_sg)
print("q through Condenser = %g BTU/lbm"%q_cond)
print("Net heat = %g BTU/lbm"%q_net)

# compute thermal efficiency
eta_th = w_net / q_sg
print("Cycle thermal efficiency = %5.2f percent."%(eta_th*100.))

# to get flow rate through the steam system and condenser cooling system
# get state point properties for sp 13-16

# state point 13 - primary core outlet
P[13] = 2190 # psia, given
T[13] = 610 # F, given
h[13] = myFluid.h_pT(P[13],T[13])

# state point 14 - primary SG outlet
P[14] = 2170 # psia, given
T[14] = 535 # F, given
h[14] = myFluid.h_pT(P[14],T[14])

# use primary/secondary energy balance for secondary flow rate
m_dot_s = m_dot_p*(h[13] - h[14])/(h[5] - h[4])

print("Mass flow rate of secondary = %g lbm/hr " % m_dot_s)

# state point 15 - condenser cooling inlet
P[15] = 14.7 # psia, given
T[15] = 91 # F, given
h[15] = myFluid.h_pT(P[15],T[15])

# state point 16 - condenser cooling outlet
P[16] = 14.7 # psia, given
T[16] = 114.9 # F, given
h[16] = myFluid.h_pT(P[16],T[16])

# use energy balance on the condenser to get flow rate
m_dot_c = m_dot_s*(-1.*q_cond)/(h[16] - h[15])

print("Mass flow rate of condenser cooling water = %g lbm/hr"%m_dot_c)
