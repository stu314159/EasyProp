#EasyProp.py
"""
Implementation Module for EasyProp.

EasyProp is another wrapper for the CoolProp Python Wrapper.

The goal is to provide an interface that is:
a) even simpler than the CoolProp Python Wrapper albeit with a reduced API; and
b) callable from a MATLAB environment with a more "MATLAB-ish" syntax.

"""

import CoolProp.CoolProp as CP


class PropertyConverter(object):
    def __init__(self):
        self.P_fact = 6.89476; # 1 psi = 6.89476 kPa
        self.h_fact = 2.326; # 1 BTU/lbm = 2.326 kJ/kg
        self.s_fact =  4.1868; # 1 BTU/lbm*R = 4.1868 kJ/kg*K
        self.c_fact = self.s_fact; # why not have two names?
        self.rho_fact = 16.0185; # 1 lbm/ft^3 = 16.0185 kg/m^3
        
        self.mu_fact = 0.020885434224573; # 1 N-s/m^2 = 0.02088... lbf-s/ft^2
        
    def mu_toUS(self,mu_SI):
        """
        convert viscosity to US units
        """
        return (mu_SI*self.mu_fact);
    
    def F_toK(self,T_F):
        """
        convert temperature in F to temperature in K
        
        """
        return ((T_F-32)*(5./9.) + 273.15)
    
    def K_toF(self,T_K):
        """
        convert temperature in K to F
        """
        
        return ((T_K-273.15)*(9./5.)+32.)
        
    def P_toUS(self,P_SI):
        """
        convert kPa to psi
        """
        return P_SI/self.P_fact;
    
    def P_toSI(self,P_US):
        """
        convert psi to kPa
        """
        return P_US*self.P_fact;
    
    def e_toUS(self,e_SI):
        """
        convert specific energy term in kJ/kg to BTU/lbm
        """
        return e_SI/self.h_fact;
    
    def e_toSI(self,e_US):
        """
        convert specific energy term in BTU/lbm to kJ/kg
        
        """
        return e_US*self.h_fact;
    
    def s_toSI(self,s_US):
        """
        convert specific entropy term in BTU/lbm*R to kJ/kg*K
        """
        return s_US*self.s_fact;
    
    def s_toUS(self,s_SI):
        """
        convert specific entropy term in kJ/kg*K
        to BTU/lbm*R
        """
        return s_SI/self.s_fact;
    
    def rho_toUS(self,rho_SI):
        """
        convert kg/m^3 to lbm/ft^3
        """
        return rho_SI/self.rho_fact;
    
    def rho_toSI(self,rho_US):
        """
        convert lbm/ft^3 to kg/m^3
        """
        return rho_US*self.rho_fact;
    
    


class EasyProp(object):
    def __init__(self,fluidName='Water',UnitSystem='SI'):
        self.fluidName = fluidName;
        if UnitSystem=='SI': 
            self.ConvertUnits=False;
        else:
            self.ConvertUnits=True;
            
        self.converter = PropertyConverter();
    
    def hL_p(self,p):
        """
        get enthalpy of saturated liquid as a function of pressure.
        
        input:
        p - pressure.  SI units: kPa; USCS units: psia
        
        output:
        h - specific enthalpy.  SI units: kJ/kg; USCS units: BTU/lbm
        """
        
        if self.ConvertUnits==False:
            p*=1000.
        else:
            p = self.converter.P_toSI(p)*1000. 
        
        
        value = CP.PropsSI('H','P',p,'Q',0,self.fluidName)
        
        if self.ConvertUnits==False:
            value = value/1000.; #<- convert J/kg to kJ/kg
        else:
            value = self.converter.e_toUS(value/1000.)
        
        return value
    
    def hV_p(self,p):
        """
        get enthalpy of saturated vapor as a function of pressure.
        
        input:
        p - pressure.  SI units: kPa; USCS units: psia
        
        output:
        h - specific enthalpy.  SI units: kJ/kg; USCS units: BTU/lbm
        """
        
        if self.ConvertUnits==False:
            p*=1000. # from kPa to Pa
        else:
            p = self.converter.P_toSI(p)*1000. # from psi to kPa to Pa
        
        
        value = CP.PropsSI('H','P',p,'Q',1.0,self.fluidName)
        
        if self.ConvertUnits==False:
            value = value/1000.; #<- convert J/kg to kJ/kg
        else:
            value = self.converter.e_toUS(value/1000.)
        
        return value
    
    
    def sL_p(self,p):
        """
        get entropy of saturated liquid as a function of pressure.
        
        input:
        p - pressure.  SI units: kPa; USCS units: psia
        
        output:
        s - specific entropy.  SI units: kJ/kg*K; USCS units: BTU/lbm*R
        """
        
        if self.ConvertUnits==False:
            p*=1000. # from kPa to Pa
        else:
            p = self.converter.P_toSI(p)*1000. # from psi to kPa to Pa
        
        
        value = CP.PropsSI('S','P',p,'Q',0,self.fluidName)
        
        if self.ConvertUnits==False:
            value = value/1000.; #<- convert J/kg to kJ/kg
        else:
            value = self.converter.s_toUS(value/1000.)
        
        return value
        
    def sV_p(self,p):
        """
        get entropy of saturated vapor as a function of pressure.
        
        input:
        p - pressure.  SI units: kPa; USCS units: psia
        
        output:
        s - specific entropy.  SI units: kJ/kg*K; USCS units: BTU/lbm*R
        """
        
        if self.ConvertUnits==False:
            p*=1000. # from kPa to Pa
        else:
            p = self.converter.P_toSI(p)*1000. 
        
        
        value = CP.PropsSI('S','P',p,'Q',1.,self.fluidName)
        
        if self.ConvertUnits==False:
            value = value/1000.; #<- convert J/kg to kJ/kg
        else:
            value = self.converter.s_toUS(value/1000.)
        
        return value
    
    
    def h_Tx(self,T,x):
        """
        get enthalpy of saturated mixture as a function of temperature and quality
        """
        if self.ConvertUnits==False:
            T += 273.15; # C to K
        else:
            T = self.converter.F_toK(T)
            
        value = CP.PropsSI('H','T',T,'Q',x,self.fluidName)
        
        if self.ConvertUnits==False:
            value /= 1000.;
        else:
            value = self.converter.e_toUS(value/1000.)
        
        return value
    
    def h_px(self,P,x):
        """
        get enthalpy of a saturated mixture as a function of pressure and quality
        
        implement in terms of previously existing functions.
        """
        
        T = self.Tsat_p(P);
        h = self.h_Tx(T,x);
        
        return h
    
    def s_px(self,P,x):
        """
        get entropy of a saturated mixture as a function of pressure and quality
        
        implement in terms of previously existing functions.
        """
        T = self.Tsat_p(P);
        s = self.s_Tx(T,x);
        
        return s
    
    def u_Tx(self,T,x):
        """
        get internal energy of saturated mixture as a function of temperature and quality
        """
        if self.ConvertUnits==False:
            T += 273.15; # C to K
        else:
            T = self.converter.F_toK(T)
            
        value = CP.PropsSI('U','T',T,'Q',x,self.fluidName)
        
        if self.ConvertUnits==False:
            value /= 1000.;
        else:
            value = self.converter.e_toUS(value/1000.)
        
        return value
    
    def s_Tx(self,T,x):
        """
        get specific entropy of a saturated mixture as a function of temp and quality
        """
        if self.ConvertUnits==False:
            T += 273.15; # C to K
        else:
            T = self.converter.F_toK(T)
            
        value = CP.PropsSI('S','T',T,'Q',x,self.fluidName)
        
        if self.ConvertUnits==False:
            value /= 1000.;
        else:
            value = self.converter.s_toUS(value/1000.)
        
        return value
    
    def h_ps(self,p,s):
        """
        return enthalpy as a function of pressure and entropy
        s = kJ/kg*K or BTU/lbm
        p = kPa or psi
        """
        if self.ConvertUnits==False:
            p*=1000.  # from kPa to Pa
            s*=1000. # from kJ/kg*K to J/kg*K
        else:
            p = self.converter.P_toSI(p)*1000. # from psi to kPa to Pa
            s = self.converter.s_toSI(s)*1000. # from BTU/lbm to kJ/kg*K to J/kg*K
        
        value = CP.PropsSI('H','P',p,'S',s,self.fluidName)
        
        if self.ConvertUnits==False:
            value/=1000.; # from J/kg to kJ/kg
        else:
            value = self.converter.e_toUS(value/1000.)
            
        return value
    
    def h_pT(self,p,T):
        """
        return specific enthalpy as a function of pressure and temperature
        """
        
        if self.ConvertUnits==False:
            p*=1000. # from kPa to Pa
            T+=273.15 # from C to K
        else:
            p = self.converter.P_toSI(p)*1000 # from psi to kPa to Pa
            T = self.converter.F_toK(T)
            
        value = CP.PropsSI('H','P',p,'T',T,self.fluidName)
        
        if self.ConvertUnits==False:
            value/=1000.
        else:
            value = self.converter.e_toUS(value/1000.)
            
        return value
        
    def s_pT(self,p,T):
        """
        return specific entropy as a function of pressure and temperature
        """
        
        if self.ConvertUnits==False:
            p*=1000. # from kPa to Pa
            T+=273.15 # from C to K
        else:
            p = self.converter.P_toSI(p)*1000 # from psi to kPa to Pa
            T = self.converter.F_toK(T)
            
        value = CP.PropsSI('S','P',p,'T',T,self.fluidName)
        
        if self.ConvertUnits==False:
            value/=1000.
        else:
            value = self.converter.s_toUS(value/1000.)
            
        return value
    
    def T_ps(self,p,s):
        """
        return temperature as a function of pressure and entropy
        """
        if self.ConvertUnits==False:
            p*=1000. # from kPa to Pa
            s*=1000. # from kJ/kg*K to J/kg*K
        else:
            p = self.converter.P_toSI(p)*1000 # from psi to kPa to Pa
            s = self.converter.s_toSI(s)*1000 # from BTU/lbm*R to kJ/kg*k to J/kg*k
        
        value = CP.PropsSI('T','P',p,'S',s,self.fluidName)
        
        if self.ConvertUnits==False:
            value -= 273.15 # K to C
        else:
            value = self.converter.K_toF(value)
            
        return value
    
    def Tsat_p(self,p):
        """
        return saturation temperature as a function of pressure
        """
        if self.ConvertUnits==False:
            p*=1000. # from kPa to Pa
        else:
            p = self.converter.P_toSI(p)*1000 # from psi to kPa to Pa
            
        
        value = CP.PropsSI('T','P',p,'Q',0.5,self.fluidName)
        
        if self.ConvertUnits==False:
            value -= 273.15 # K to C
        else:
            value = self.converter.K_toF(value)
            
        return value
    
    
    def Psat_T(self,T):
        """
        return saturation pressure as a function of temperature
        """
        
        if self.ConvertUnits==False:
            T+=273.15 # convert C to K
        else:
            T = self.converter.F_toK(T)
            
        value = CP.PropsSI('P','T',T,'Q',0.5,self.fluidName);
        
        if self.ConvertUnits==False:
            value/=1000. # from Pa to kPa
        else:
            value = self.converter.P_toUS(value/1000.); # Pa to kPa to psi
            
        return value
    
    def T_ph(self,p,h):
        """
        return temperature as a function of pressure and enthalpy
        """
        if self.ConvertUnits==False:
            p*=1000. # from kPa to Pa
            h*=1000. # from kJ/kg to J/kg
        else:
            p = self.converter.P_toSI(p)*1000 # from psi to kPa to Pa
            h = self.converter.e_toSI(h)*1000 # from BTU/lbm to kJ/kg to J/kg
        
        value = CP.PropsSI('T','P',p,'H',h,self.fluidName)
        
        if self.ConvertUnits==False:
            value -= 273.15 # K to C
        else:
            value = self.converter.K_toF(value)
            
        return value
    
    def s_ph(self,p,h):
        """
        return entropy as a function of pressure and enthalpy
        """
        if self.ConvertUnits==False:
            p*=1000. # from kPa to Pa
            h*=1000. # from kJ/kg to J/kg
        else:
            p = self.converter.P_toSI(p)*1000 # from psi to kPa to Pa
            h = self.converter.e_toSI(h)*1000 # from BTU/lbm to kJ/kg to J/kg
        
        value = CP.PropsSI('S','P',p,'H',h,self.fluidName)
        
        if self.ConvertUnits==False:
            value /= 1000. # J/kg*K to kJ/kg*K
        else:
            value = self.converter.s_toUS(value/1000.) # J/kg*K to kJ/kg*K to BTU/lbm*R
            
        return value
    
    def x_ph(self,p,h):
        """
        return quality as a function of pressure and enthalpy
        """
        if self.ConvertUnits==False:
            p*=1000. # from kPa to Pa
            h*=1000. # from kJ/kg to J/kg
        else:
            p = self.converter.P_toSI(p)*1000 # from psi to kPa to Pa
            h = self.converter.e_toSI(h)*1000 # from BTU/lbm to kJ/kg to J/kg
        
        value = CP.PropsSI('Q','P',p,'H',h,self.fluidName)
        
        return value
    
    def v_ph(self,p,h):
        """
        return specific volume as a function of pressure and enthalpy
        """
        if self.ConvertUnits==False:
            p*=1000. # from kPa to Pa
            h*=1000. # from kJ/kg to J/kg
        else:
            p = self.converter.P_toSI(p)*1000 # from psi to kPa to Pa
            h = self.converter.e_toSI(h)*1000 # from BTU/lbm to kJ/kg to J/kg
        
        value = CP.PropsSI('D','P',p,'H',h,self.fluidName)
        
        if self.ConvertUnits==True:
            value = self.converter.rho_toUS(value)
            
        
        return (1./value)
    
    def v_Tx(self,T,x):
        """
        return specific volume as a function of temperature and quality
        """
        
        if self.ConvertUnits==False:
            T-=273.15 # C to K
        else:
            T = self.converter.F_toK(T)
        
        value = CP.PropsSI('D','T',T,'Q',x,self.fluidName)
        
        if self.ConvertUnits==True:
            value = self.converter.rho_toUS(value)
                    
        return (1./value)
    
    def v_pT(self,p,T):
        """
        return specific volume as a function of pressure and temperature
        """
        
        if self.ConvertUnits==False:
            p*=1000. # from kPa to Pa
            T+=273.15 # from C to K
        else:
            p = self.converter.P_toSI(p)*1000 # from psi to kPa to Pa
            T = self.converter.F_toK(T)
            
        value = CP.PropsSI('D','P',p,'T',T,self.fluidName)
        
        if self.ConvertUnits==True:
            value = self.converter.rho_toUS(value)
                
        return (1./value)
    
    def P_vT(self,v,T):
        """
        return pressure as a function of specific volume and Temperature
        """
        rho = 1./v # kg/m^3
        if self.ConvertUnits==False:
            T+=273.15 # from C to K
        else:
            T = self.converter.F_toK(T)
            rho = self.converter.rho_toSI(rho) # lbm/ft^3 to kg/m^3
            
        value = CP.PropsSI('P','D',rho,'T',T,self.fluidName)
        
        if self.ConvertUnits==False:
            value/=1000. # from Pa to kPa
        else:
            value /=1000.
            value = self.converter.P_toUS(value)
            
        return value
    
    def T_pv(self,p,v):
        """
        return temperature as a function of specific volume and pressure
        """
        rho = 1./v
        if self.ConvertUnits==False:
            p*=1000. # kPa to Pa
        else:
            p = self.converter.P_toSI(p)*1000. # from psi to kPa to Pa
            rho = self.converter.rho_toSI(rho) # from lbm/ft^3 to kg/m^3
            
        value = CP.PropsSI('T','P',p,'D',rho,self.fluidName)
        
        if self.ConvertUnits==False:
            value-=273.15 # from K to C
        else:
            value = self.converter.K_toF(value)
            
        return value
            
    def P_crit(self):
        """
        return critical pressure of the fluid
        """
        pC = CP.PropsSI('PCRIT',self.fluidName)
        
        if self.ConvertUnits==False:
            value = pC/1000. # convert Pa to kPa
        else:
            value = self.converter.P_toUS(pC/1000.)
            
        return value
    
    def T_crit(self):
        """
        return critical temperature of the fluid
        """
        
        cT = CP.PropsSI('TCRIT',self.fluidName)
        
        if self.ConvertUnits==False:
            value = cT-273.15; # convert K to C
        else:
            value = self.converter.K_toF(cT)
            
        return value
    
        
    def R_pT(self,p,T):
        """
        return Ideal Gas constant as a function of pressure and temperature
        
        """
        
        if self.ConvertUnits==False:
            p*=1000. # from kPa to Pa
            T+=273.15 # from C to K
        else:
            p = self.converter.P_toSI(p)*1000 # from psi to kPa to Pa
            T = self.converter.F_toK(T)
            
        R_bar = CP.PropsSI('GAS_CONSTANT','P',p,'T',T,self.fluidName) # J/mol-K
        M_bar = CP.PropsSI('M','P',p,'T',T,self.fluidName) # kg/mol
        
        value = R_bar/M_bar # J/kg-K
        value /=1000. #kJ/kg-K
        
        if self.ConvertUnits==False:
            pass
        else:
            value = self.converter.s_toUS(value) # note this is now BTU/lbm-R
            # not ft-lbf/lbm*R as may be desired for some ideal
            # gas applications.  
        
        
        return value
    
    def mu_pT(self,p,T):
        """
        return the viscosity (Pa-s, kg/m-s  or lbm/ft-s)
        """
        
        if self.ConvertUnits==False:
            p*=1000. # from kPa to Pa
            T+=273.15 # from C to K
        else:
            p = self.converter.P_toSI(p)*1000 # from psi to kPa to Pa
            T = self.converter.F_toK(T)
            
        value = CP.PropsSI('V','P',p,'T',T,self.fluidName) # kg/m-s
        
        if self.ConvertUnits==True:
            value = self.converter.mu_toUS(value)
            
        return value
    
    def Prandtl_pT(self,p,T):
        """
        return Prandtl number as a function of pressure and temperature
        """
        if self.ConvertUnits==False:
            p*=1000. # from kPa to Pa
            T+=273.15 # from C to K
        else:
            p = self.converter.P_toSI(p)*1000 # from psi to kPa to Pa
            T = self.converter.F_toK(T)
            
        value = CP.PropsSI('PRANDTL','P',p,'T',T,self.fluidName) # no units
        
        return value
       
class simpleFluid(EasyProp):
    """
    A new class inheriting EasyProp.  Partly to give a more logical class name, partly to
    make things less confusing when a simpleMixture class is added.
    
    """
    
    def __init__(self,fluidName,unitSystem):
        EasyProp.__init__(self,fluidName,unitSystem);
 
class fluidSpecies(object):
    """
    
    """
    def __init__(self,fluidName,unitSystem,weight):
        self.fluidName = fluidName
        self.unitSystem = unitSystem
        self.weight = weight
        
        self.fluid = simpleFluid(fluidName,unitSystem)
    
class simpleMixture(object):
    """
    A class to provide mass- or molar-weighted mixtures of fluids.  
    """
    def __init__(self,mixtureDict={'Water':1.0},weighting='w/o',units='SI'):
        """
        constructor:
        mixtureDict = dictionary containing the mixture fluid names and percent weighting
         -- initially only implement w/o; add a/o later.
         weighting = 'w/o' - weight percent; 'a/o' atom percent
         units = 'SI' | 'USCS'
        """
        self.mixtureDict = mixtureDict;
        self.weighting = weighting;
        self.units = 'SI'
        self.fluidDict = {}
        mixIndex = 0
        for species in mixtureDict.keys():
            self.fluidDict[mixIndex] = fluidSpecies(species,units,mixtureDict[species])
            mixIndex+=1
        
    def get_components(self):
        return self.fluidDict;
        
        