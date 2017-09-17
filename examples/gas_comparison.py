# gas comparison based on analysis of section 7.4 of "Nuclear Engery Conversion" by El-Wakil

import numpy as np

import sys
sys.path.insert(0,'..') #<-- assumes EasyProp is one folder up

import EasyProp as EP

numFluids = 6;

physProp = np.zeros((numFluids,2),dtype=np.float)

def getPhysProp(fluid,p,T):
    M = fluid.M()
    g = fluid.gamma_pT(p,T)
    Pr = fluid.Prandtl_pT(p,T)
    value = (Pr**0.6)*(M)*(((g-1.)/g)**3.)
    return value
    

# fluid 1 - Hydrogen
# fluid 2 - Helium
# fluid 3 - CO2
# fluid 4 - Nitrogen
# fluid 5 - Air
# fluid 6 - Argon


# construct fluids
h2 = EP.simpleFluid('Hydrogen','USCS')
he = EP.simpleFluid('Helium','USCS')
co2 = EP.simpleFluid('CO2','USCS')
n2 = EP.simpleFluid('Nitrogen','USCS')
air = EP.simpleFluid('Air','USCS')
argon = EP.simpleFluid('argon','USCS')
#xenon = EP.simpleFluid('xenon','USCS')

# make fluid list
fluidList = [h2,he,co2,n2,air,argon]

T = 80.;
p = 14.7;
for i in range(len(fluidList)):
    physProp[i][0] = getPhysProp(fluidList[i],p,T);

T = 600.;
p = 14.7;
for i in range(len(fluidList)):
    physProp[i][1] = getPhysProp(fluidList[i],p,T);
    

# get normalization constant
normC = 1./physProp[0][0]; # so that normC*physProp = 1.0 for Hydrogen at 80 F

physProp *= normC;

print physProp




