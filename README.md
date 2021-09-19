# FlexWingROM
 
The matlab codebase FlexWing-ROM models fluid structure interaction (FSI) of flexible wings and generates data-driven reduced order models for control of flexible wings.  

Details about the different methods can be found in the JOSS paper Fasel et al. (2021).


## How to run the code

[main.m](/main.m)  

## What you can do with this code

FlexWingROM is a matlab codebase that contains a fully parametrized wing model generator, a fluid-structure interaction solver, and three recent data-driven parametric reduced order modeling methods for flexible wings.

### Parametrized wing model generator and fluid structure interaction simulation of flexible wings

First, a flexible wing FSI model is generated that is coupling a finite element code with a 3D unsteady panel method. The wing design is fully parametrised and can be defined in [code/wingDesignAndSimParameters.m](/code/wingDesignAndSimParameters.m) (e.g. the airfoil shape, the planform, the material properties, structural design, ...).  
The open source FE-code YetAnotherFEcode and parts of the Apame 3D panel code and XFOIL are used:
* Finite element code: 
  * Shobhit Jain, Jacopo Marconi & Paolo Tiso (2020). YetAnotherFEcode. Zenodo. 
  * Paper: http://doi.org/10.5281/zenodo.4011281    
  * Code: https://github.com/jain-shobhit/YetAnotherFEcode
    
* 3D unsteady panel method
  * Based on J Katz, A Plotkin. Low-speed aerodynamics, Cambridge university press, 2001
  * Steady panel method matlab implementation: http://www.3dpanelmethod.com/ 
    
* XFOIL
  * M. Drela. XFOIL, https://web.mit.edu/drela/Public/web/xfoil/ (XFOIL MATLAB interface 2011 by Rafael Oliveira) 
   
   
### FSI test cases and comparisons

The FSI model of the flexible wing is then used to run some test cases 
* Modal vs. full FE-model FSI
* Unsteady FSI
* Comparison unsteady panel method vs. Theodorsen's function

### Data driven reduced order models

Finally, three data-driven (parametric) reduced order modeling approaches are implemented and compared:
* algebraic dynamic mode decomposition with control 
* input output dynamic mode decomposition
* balanced mode decomposition


## Code based on following publications

The code is mainly based on the following three publications:

* U. Fasel, P. Tiso, D. Keidel, G. Molinari, P. Ermanni. Reduced-order dynamic model of a morphing airborne wind energy aircraft.
  
  https://arc.aiaa.org/doi/abs/10.2514/1.J058019
     
* N. Fonzi, S. L. Brunton, U. Fasel. Data-driven nonlinear aeroelastic models of morphing wings for control
  
  https://royalsocietypublishing.org/doi/pdf/10.1098/rspa.2020.0079
     
* A. Iannelli, U. Fasel, R. S. Smith. The balanced mode decomposition algorithm for data-driven LPV low-order models of aeroservoelastic systems
  
  https://www.sciencedirect.com/science/article/pii/S127096382100331X
