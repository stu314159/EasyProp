# EasyProp
a simplified inferface to the CoolProp Python Wrapper 

This interface also provides for using SI (C,kJ/kg,kPa,etc..) or USCS (psia, F, BTU/lbm, etc...)
transparently.

This package depends on CoolProp Python wrapper and its dependencies.

Install by: pip install CoolProp

Of course, you can use CoolProp without this package.

The intent is to use this package with optional access through MATLAB.  The resulting
interface will provide fluid properties with a callin syntax similar to what you would find
in some other MATLAB fluid property libraries like XSteam.  Use of this
package will allow similar functionality but brings with it the broad range of fluid
properties provided by CoolProp.
