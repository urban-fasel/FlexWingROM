# FlexWing-ROM
 
The matlab codebase **FlexWing-ROM** models fluid structure interaction (FSI) of flexible wings and generates data-driven reduced order models that can be used for analysis and control.  

The different methods are introduced in the arXiv paper Fasel et al. (2021) -> link will follow.

![FlexWing-ROM_OverviewFigure](/docs/FlexWing-ROM_OverviewFigure.png)  

## Table of contents
* [How to run the code](#how-to-run-the-code)
 	* [Main script](#main-script)
	* [Tutorials](#tutorials)
	* [Dependencies](#dependencies)
* [What you can do with this code](#what-you-can-do-with-this-code)
 	* [Parametrized wing model generator and fluid structure interaction solver](#parametrized-wing-model-generator-and-fluid-structure-interaction-solver)
 	* [FSI test cases and comparisons](#fsi-test-cases-and-comparisons)
 	* [Data driven (parameter varying) reduced order models](#data-driven-parameter-varying-reduced-order-models)
* [Community guidelines](#community-guidelines)
	* [Contribute to the software](#contribute-to-the-software)
	* [Report problems with the software or seek support](#report-problems-with-the-software-or-seek-support)
* [References](#references)


## How to run the code

### Main script
[MAIN.m](/MAIN.m)  
Generates a wing model (e.g. NACA0012 or NACA6418, morphing or non-morphing), runs different FSI test cases, and generates and compares reduced-order models.  


![FlexWing-ROM_CodeStructure](/docs/FlexWing-ROM_CodeStructure.png)  

&nbsp;  

### Tutorials

Two different flexible wings are used in the tutorials, shown in the figure below. On the left, a NACA6418 morphing wing (five compliant ribs on each side), and on the right, a NACA0012 non-morphing wing (no compliant ribs). 

![Flexible wings](/docs/FlexWing-ROM_WingExamples.png)

* [example1_NACA0012_FSI_modal_vs_displacement.m](/example1_NACA0012_FSI_modal_vs_displacement.m)   
	* Comparing a modal vs. displacement FE-model FSI (on a NACA0012).

* [example2_NACA6418_unsteadyFSI.m](/example2_NACA6418_unsteadyFSI.m)  
	* Runs an unsteady panel method FSI (on a NACA6418) and animates the wake and wing displacements. The animation shows sinusoidal pitching (left) and morphing (right)

![FlexWing-ROM_NACA6418_animation.gif](/docs/FlexWing-ROM_NACA6418_animation.gif) 


* [example3_NACA0012_unsteadyPM_vs_Theodorsen.m](/example3_NACA0012_unsteadyPM_vs_Theodorsen.m)   
	* Comparing an unsteady panel method (on a NACA0012) with Theodorsen's function.

* [example4_NACA6418_ROM.m](/example4_NACA6418_ROM.m)  
	* Comparing three data-driven reduced-order models using a NACA6412 morphing wing.


### Dependencies

The following matlab toolboxes are required to run the code:

* Control System Toolbox
* System Identification Toolbox
* Statistics and Machine Learning Toolbox
* Parallel Computing Toolbox

XFOIL is used to calculate the viscous drag and maximum lift coefficient of the airfoil. The XFOIL version only runs on Windows. In case the code is run on Linux, the airfoil coefficients are loaded from precalculated tables for NACA0012 and NACA6418.


## What you can do with this code

FlexWingROM is a matlab codebase that contains a fully parametrized wing model generator, a fluid-structure interaction solver, and three recent data-driven parametric reduced order modeling methods for flexible wings.

### Parametrized wing model generator and fluid structure interaction solver

The main code structure is shown in the figure above. First, a flexible wing FSI model is generated that is coupling a finite element code with a 3D unsteady panel method. The wing design is fully parametrised and can be defined in [wingDesignAndSimParameters.m](/code/generateModel/wingDesignAndSimParameters.m) (e.g. the airfoil shape (currently four digit NACAs implemented), the planform, the material properties, structural design, ...).  
The open source FE-code YetAnotherFEcode and parts of the Apame 3D panel code and XFOIL are used:
* Finite element code: 
  * Shobhit Jain, Jacopo Marconi & Paolo Tiso (2020). [YetAnotherFEcode](http://doi.org/10.5281/zenodo.4011281). Zenodo. 
    
* 3D unsteady panel method
  * Based on J Katz, A Plotkin. [Low-speed aerodynamics](https://www.cambridge.org/core/books/lowspeed-aerodynamics/077FAF851C4582F1B7593809752C44AE), Cambridge university press, 2001
  * Steady panel method matlab implementation: [Apame](http://www.3dpanelmethod.com/) 
    
* XFOIL
  * M. Drela. [XFOIL](https://web.mit.edu/drela/Public/web/xfoil/) 
  * [Matlab interface](https://www.mathworks.com/matlabcentral/fileexchange/30478-rafael-aero-xfoilinterface) 2011 by Rafael Oliveira
   
   
### FSI test cases and comparisons

The FSI model of the flexible wing is then used to run some test cases 
* Modal vs. full FE-model FSI [runSteadyFSItestcases.m](/code/FSI/runSteadyFSItestcases.m)
* Unsteady FSI [runUnsteadyFSItestcases.m](/code/FSI/runUnsteadyFSItestcases.m)
* Comparison unsteady panel method vs. Theodorsen's function [runTheodorsenFSItestcases.m](/code/FSI/runTheodorsenFSItestcases.m)

### Data driven (parameter varying) reduced order models

The main features of the code are three data-driven (parameter varying) reduced order modeling approaches ([ROM.m](/code/ROM/ROM.m)):
* algebraic dynamic mode decomposition with control 
* input output reduced-order model
* balanced mode decomposition

The FSI solver is used to generate synthetic impulse response data sets that are stored in data files. These are used to build and test the three different reduced order models. Each method can be used to first generate local linear models that are accurate around a single operating condition. Then, depending on each method, interpolation schemes are introduced to build parameter varying models that are accurate over multiple operating conditions. The accuracy of the method is evaluated and the different methods are compared.  
Most importantly, the model order reduction methods can be used with other fluid or fluid structure data sets (both numerical and experimental), by loading any external impulse response data. Therefore, the code should be widely applicable and useful for generating accurate and efficient reduced order models.


## Community guidelines

### Contribute to the software

If you would like to contribute an example or extend a method (e.g. new ROM variants or advanced functionalities), reach out to us by creating an issue!

### Report problems with the software or seek support

If you find a bug in the code or need help using FlexWingROM please create an issue.


## References

The code is mainly based on the following publications:
 
* N. Fonzi, S. L. Brunton, U. Fasel. [Data-driven nonlinear aeroelastic models of morphing wings for control](https://royalsocietypublishing.org/doi/pdf/10.1098/rspa.2020.0079). Proceedings of the Royal Society A, 2020. 
     
* A. Iannelli, U. Fasel, R. S. Smith. [The balanced mode decomposition algorithm for data-driven LPV low-order models of aeroservoelastic systems](https://www.sciencedirect.com/science/article/pii/S127096382100331X). Aerospace Science and Technology, 2021. 

* J. Annoni, P. Seiler. [A method to construct reduced‐order parameter‐varying models](https://onlinelibrary.wiley.com/doi/am-pdf/10.1002/rnc.3586). International Journal of Robust and Nonlinear Control, 2017. 

* G. Molinari, A. F. Arrieta, and P. Ermanni. [Aero-structural optimization of three-dimensional adaptive wings with embedded smart actuators](https://arc.aiaa.org/doi/abs/10.2514/1.J052715). AIAA Journal, 2014. 

* U. Fasel, P. Tiso, D. Keidel, G. Molinari, P. Ermanni. [Reduced-order dynamic model of a morphing airborne wind energy aircraft](https://arc.aiaa.org/doi/abs/10.2514/1.J058019). AIAA Journal, 2019. 
    
