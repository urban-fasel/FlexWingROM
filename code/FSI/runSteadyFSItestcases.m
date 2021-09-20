function results = runSteadyFSItestcases(paramFSI, simParam, wingModelStructure, wingDesign, wingModelAero)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% compare modal (K_MODES) vs. full (K_ASET) model, running FSI with steady panel method
%

paramFSI.aeroOnly = 0;      % run aerodynamics only (no wing deformation): false for FSI

% run steady FSI test case: comparison of full stiffness matrix with modal approach
simSteadyParam.V = 30;              % flight velocity
simSteadyParam.alpha = 6;           % angle of attack
simSteadyParam.forceACTinp.L = 1;   % actuator force unsymmetric actuation (roll actuation), [-1 1]
simSteadyParam.forceACTinp.fs = 0;  % actuator force symmetric actuation (cL control), [-1 1]

simSteadyParam.doPlot = paramFSI.plt; % plot deformed wing

% run FSI using modal stiffness matrix
FSI_steady_out = FSI_steady(paramFSI, simSteadyParam);

% run FSI using full stiffness matrix
FSI_steady_fullK_out = FSI_steady_fullK(paramFSI, simSteadyParam, simParam, wingModelStructure, wingDesign, wingModelAero);

% compare outputs
results.liftError_steady = abs(FSI_steady_out.cL - FSI_steady_fullK_out.cL)/FSI_steady_fullK_out.cL;
results.rollError_steady = abs(FSI_steady_out.cRoll - FSI_steady_fullK_out.cRoll)/FSI_steady_fullK_out.cRoll;
  
% save outputs
results.FSI_steady_out = FSI_steady_out;
results.FSI_steady_fullK_out = FSI_steady_fullK_out;


