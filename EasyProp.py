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