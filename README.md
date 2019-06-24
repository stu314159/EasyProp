# EasyProp
a simplified inferface to the CoolProp Python Wrapper 

This interface extends the CoolProp Python wrapper in that two unit systems are allowed:
1. SI (C,kJ/kg,kPa,etc..); or
2. USCS (psia, F, BTU/lbm, etc...)


CoolProp is also extended in that functionality is provided to model simple (non-interacting) mixtures
(using either a/o or w/o).

This package explicitly depends upon the CoolProp Python wrapper.  This dependency can be installed
with the following command: pip install CoolProp

EasyProp also presents a somewhat cleaner function calling syntax relative to CoolProp, however this admittedly
comes at the expense of incomplete functionality.  If there is something missing, however, it is quite 
straight-forward to add it.

Some example scripts are provided to illustrate how to use EasyProp for solving typical problems
in thermodynamic cycle analysis.  

Requires Python 3 (tested with 3.6.7), numpy, scipy, and coolprop libraries --- all installable via, for example, conda.
