#simpleOtto.py
"""
This example is adapted from Example 9.1 from Fundamentals of Engineering Thermodynamics
Moran, Shapiro, et al. 7th edition.

The temperature at the beginning of the compression process of an air-standard Otto cycle with a compression ratio of 8 is 540 R, the pressure is 1 atm an the cylinder volume is
0.02 ft**3.  The maximum temperature during the cycle is 3600 R.  Determine:

a) the temperature and pressure at the end of each process of the cycle; 
b) the thermal efficiency, and 
c) the mean effective pressure

Instead of an air-standard analysis; we will assume variable specific heats

"""

import numpy as np
import sys
sys.path.insert(0,'..')

import EasyProp as EP

air = EP.simpleFluid('Air','USCS');

numSP = 5; # reserve SP0 for the dead state
# declare fluid property arrays
u = np.zeros((numSP,),dtype=np.float) # BTU/lbm , internal energy
u_s = np.zeros_like(u);
s = np.zeros_like(u); # BTU/lbm-R, entropy
s_s = np.zeros_like(u);
P = np.zeros_like(u); # psia, pressure
T = np.zeros_like(u); # F, temperature
T_s = np.zeros_like(u); 
v = np.zeros_like(u); # specific volume


r_v = 8.0 # compression ratio
V = 0.02 # ft**3, volume of cylinder, given

# establish dead state condition
T[0] = 80.
P[0] = 14.7


# state point 1, beginning of compression stroke
T[1] = 80.
P[1] = 14.7
u[1] = air.u_pT(P[0],T[0])
s[1] = air.s_pT(P[0],T[0])
v[1] = air.v_pT(P[0],T[0])

# process 1 -> 2, isentropic compression
s[2] = s[1]
v[2] = v[1]/r_v
T[2] = air.T_sv(s[2],v[2])
P[2] = air.P_vT(v[2],T[2])
u[2] = air.u_pT(P[2],T[2])

w_comp = u[1] - u[2]

# process 2 -> 3, isometric heat addition to 3600 R (3140 F)
v[3] = v[2]
T[3] = 3140. # F, given
P[3] = air.P_vT(v[3],T[3])
u[3] = air.u_pT(P[3],T[3])
s[3] = air.s_pu(P[3],u[3])

q_s = u[3] - u[2]

# process 3 --> 4, isentropic expansion to v[4] = v[1]
s[4] = s[3]
v[4] = v[1]
T[4] = air.T_sv(s[4],v[4])
P[4] = air.P_vT(v[4],T[4])
u[4] = air.u_pT(P[4],T[4])

w_exp = u[3] - u[4]

# process 4 --> 1, isometric heat rejection
q_r = u[1] - u[4]

q_net = q_r + q_s
w_net = w_comp + w_exp
eta_th = w_net/q_s

m_air = V/v[1]; # lbm, based on cylinder volume (ft**3) and initial air specific volume

W_net = w_net*m_air

print("T = ", T)
print("u = ", u)
print("v = ", v)
print("Net specific heat transferred = %g BTU/lbm"%q_net)
print("Net specific work = %g BTU/lbm"%w_net)
print("Net total work = %g BTU"%W_net)
print("Thermal efficiency = %g percent"%(eta_th*100.))
