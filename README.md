# FlexWing-ROM
 
The matlab codebase **FlexWing-ROM** models fluid structure interaction (FSI) of flexible wings and generates data-driven reduced order models that can be used for analysis and control.  

![FlexWing-ROM_OverviewFigure](/data/FlexWing-ROM_OverviewFigure.png)  

Details about the different methods can be found in the JOSS paper Fasel et al. (2021) -> link will follow.


## How to run the code

[MAIN.m](/MAIN.m)  

## What you can do with this code

FlexWingROM is a matlab codebase that contains a fully parametrized wing model generator, a fluid-structure interaction solver, and three recent data-driven parametric reduced order modeling methods for flexible wings.

### Parametrized wing model generator and fluid structure interaction solver

First, a flexible wing FSI model is generated that is coupling a finite element code with a 3D unsteady panel method. The wing design is fully parametrised and can be defined in [wingDesignAndSimParameters.m](/code/generateModel/wingDesignAndSimParameters.m) (e.g. the airfoil shape, the planform, the material properties, structural design, ...).  
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
* Modal vs. full FE-model FSI [runSteadyFSItestcases.m](/code/FSI/runSteadyFSItestcases.m)
* Unsteady FSI [runUnsteadyFSItestcases.m](/code/FSI/runUnsteadyFSItestcases.m)
* Comparison unsteady panel method vs. Theodorsen's function [runTheodorsenFSItestcases.m](/code/FSI/runTheodorsenFSItestcases.m)

### Data driven (parameter varying) reduced order models

The main features of the code are the three data-driven (parameter varying) reduced order modeling approaches ([MAIN_ROM.m](/code/ROM/MAIN_ROM.m)):
* algebraic dynamic mode decomposition with control 
* input output dynamic mode decomposition
* balanced mode decomposition

The FSI solver is used to generate synthetic impulse response data sets that are stored in data files. These are used to build and test the three different reduced order models. Each method can be used to first generate local linear models that are accurate around a single operating condition. Then, depending on each method, interpolation schemes are introduced to build parameter varying models that are accurate over multiple operating conditions. The accuracy of the method is evaluated and the different methods are compared.  
Most importantly, the model order reduction methods can easily be used with other fluid or fluid structure data sets (both numerical and experimental), by simply loading any external impulse response data. Therefore, the code should be widely applicable and useful for generating accurate and efficient reduced order models.


## This code based on the following publications

The code is mainly based on the following publications:
 
* N. Fonzi, S. L. Brunton, U. Fasel. Data-driven nonlinear aeroelastic models of morphing wings for control. Proceedings of the Royal Society A, 2020.
  
  https://royalsocietypublishing.org/doi/pdf/10.1098/rspa.2020.0079
     
* A. Iannelli, U. Fasel, R. S. Smith. The balanced mode decomposition algorithm for data-driven LPV low-order models of aeroservoelastic systems. Aerospace Science and Technology, 2021.
  
  https://www.sciencedirect.com/science/article/pii/S127096382100331X

* J. Annoni, P. Seiler. A method to construct reduced‐order parameter‐varying models. International Journal of Robust and Nonlinear Control, 2017.

  https://onlinelibrary.wiley.com/doi/am-pdf/10.1002/rnc.3586

* G. Molinari, A. F. Arrieta, and P. Ermanni. Aero-structural optimization of three-dimensional adaptive wings with embedded smart actuators. AIAA Journal, 2014
 
  https://arc.aiaa.org/doi/abs/10.2514/1.J052715

* U. Fasel, P. Tiso, D. Keidel, G. Molinari, P. Ermanni. Reduced-order dynamic model of a morphing airborne wind energy aircraft. AIAA Journal, 2019.
  
  https://arc.aiaa.org/doi/abs/10.2514/1.J058019
    
