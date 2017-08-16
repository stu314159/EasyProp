%RankineCycle2_analysis.m
%Use of EasyProp from MATLAB environment
%% Problem Description
%
%  Consider the Rankine cycle depicted in the figure below:
%
% 
% <<RankineCycle2_Schematic.png>>
% 
% 
%  The following fluid conditions are given for each numbered state point:
% 
% # state point 1: P = 1.5 psia, x = 0.0
% # state point 2: P = 164 psia, MCP efficiency = 84 percent
% # state point 3: P = 164 psia, x = 0.0
% # state point 4: P = 838 psia, MFP efficiency = 89 percent
% # state point 5: P = 838 psia, x = 1.0
% # state point 6: P = 164 psia, HP Turbine efficiency = 94%
% # state point 7: P = 164 psia, x = 1.0
% # state point 8: P = 164 psia, T = 490 F
% # state point 9: P = 1.5 psia, LP Turbine efficiency = 92%
% # state point 10: P = 838 psia, x = 0.0
% # state point 11: P = 164 psia, isenthalpic expansion in trap T1
% # state point 12: P = 164 psia, x = 0.0, liquid drain from moisture
% separator
% # state point 13: P = 2190 psia, T = 610 F
% # state point 14: P = 2170 psia, T = 535 F
% # state point 15: P = 14.7 psia, T = 91 F
% # state point 16: P = 14.7 psia, T = 114.9 F
%
% 
% 
%  Assume primary coolant flow rate to be 113.5 lbm/hr.  Flow fraction 
%  f1 and f2 represent a fraction of the total steam flow 
%  _at the respective extraction point_ and is constrained to be between 0
%  and 1.
%
%  Calculate the following:
% 
% # The enthalpy at all state points and the flow fraction f1 and f2.
% # The net specific work of all turbines and pumps (BTU/lbm)
% # The net specific heat supplied at the steam generator and condenser
% (BTU/lbm)
% # Cycle thermal efficiency
% # Flow rate of steam from the Steam Generator and the rate of cooling
% water flow through the condenser.
% 
% 


%% clear the environment
clear
clc
close 'all'

%% Put EasyProp on Python path
% put location of EasyProp.py module on the python search path
% here, I assume that EasyProp.py is up one level in the file system
if count(py.sys.path,'..') == 0  % <-- see if desired directory is on path
    insert(py.sys.path,int32(0),'..'); %<-- if not; add it.
end

%% Establish working fluid and unit system
fluid = 'Water';
unitSystem = 'USCS';
myFluid = py.EasyProp.EasyProp(fluid,unitSystem);