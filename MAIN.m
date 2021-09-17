
%%%%%%%%%%
%
%% This code models fluid structure interaction of flexible wings and generates data-driven reduced order models for control of flexible wings
%
%   The code is mainly based on the following three publications:
%       - Reduced-order dynamic model of a morphing airborne wind energy aircraft
%            U. Fasel, P. Tiso, D. Keidel, G. Molinari, P. Ermanni.
%            https://arc.aiaa.org/doi/abs/10.2514/1.J058019
%       - Data-driven nonlinear aeroelastic models of morphing wings for control
%            N. Fonzi, S. L. Brunton, U. Fasel.
%            https://royalsocietypublishing.org/doi/pdf/10.1098/rspa.2020.0079
%       - The balanced mode decomposition algorithm for data-driven LPV low-order models of aeroservoelastic systems
%            A. Iannelli, U. Fasel, R. S. Smith
%            https://www.sciencedirect.com/science/article/pii/S127096382100331X
%
%
%   In this script, first, a flexible wing FSI model is generated that is coupling a finite element
%   code with a 3D unsteady panel method. The wing design is fully
%   parametrised and can be defined in the function wingDesignAndSimParameters 
%   (e.g. the airfoil shape, the planform, the material properties, structural design, ...).
%   The open source FE-code YetAnotherFEcode and parts of the Apame 3D panel code and XFOIL are used
%
%       - Finite element code: 
%           Shobhit Jain, Jacopo Marconi & Paolo Tiso (2020). YetAnotherFEcode. Zenodo. 
%           http://doi.org/10.5281/zenodo.4011281
%           https://github.com/jain-shobhit/YetAnotherFEcode
%       - 3D unsteady panel method
%           Based on J Katz, A Plotkin. Low-speed aerodynamics, Cambridge university press, 2001
%           Steady panel method matlab implementation: http://www.3dpanelmethod.com/ 
%       - XFOIL
%           M. Drela. XFOIL, https://web.mit.edu/drela/Public/web/xfoil/
%           XFOIL MATLAB interface 2011 by Rafael Oliveira 
%            
%
% 
%   The FSI model of the flexible wing is then used to run some test cases 
%        - Modal vs. full FE-model FSI
%        - Unsteady FSI
%        - Comparison unsteady panel method vs. Theodorsen's function
%
%   Finally, three data-driven (parametric) reduced order modeling approaches are implemented and compared:
%       - algebraic dynamic mode decomposition with control 
%       - input output dynamic mode decomposition
%       - balanced mode decomposition
%
%
%%%%%%%%%%

clear all
close all
clc

addpath(genpath('code'))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% WING DESIGN & SIMULATION PARAMETERS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the parameters and simulation results already stored in the folder are
% for a NACA6418 morphing wing and a NACA0012 non-morphing wing. 
% To reproduce these results, run either:
%   airfoil = 'NACA6418'; morphingWing = true; 
%   airfoil = 'NACA0012'; morphingWing = false;

generateWing = false;       % generate or load wing design parameters and FSI model
storeAllData = false;       % store all data: needed for runTestSteady
runTestSteady = false;      % run steady aerodynamics FSI test cases: need to also run generateWing and storeAllData, as the large system matrices are needed for comparison that are not stored 
runTestUnsteady = false;    % run unsteady aerodynamics FSI test cases
runTestTheodorsen = true;   % run Theodorsen comparison FSI test cases
runROM = false;             % generate reduced order models and compare different methods
 
airfoil = 'NACA0012'; % 'NACA6418';       % choose airfoil: coordinates are loaded from file in folder 'airfoils' 
morphingWing = false; % true;       % set true for wing design with compliant ribs -> morphing for roll and load control
plt = true;                % plot

% define the main wing design and simulation parameters
if generateWing
    [wingDesign, simParam] = wingDesignAndSimParameters(morphingWing, airfoil, plt);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% generate finite element and aerodynamic model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if generateWing
    % generate wing structure model 
    wingModelStructure = generateWingModelStructure(wingDesign,simParam);

    % generate wing aero model
    wingModelAero = generateWingModelAero(wingDesign,simParam,wingModelStructure);

    % plot finite element mesh, and vibration and morphing deformation modes
    if plt
        plotModel(wingDesign, simParam, wingModelAero, wingModelStructure)
    end
    
    % save parameters that are used in FSI
    paramFSI = saveParamFSI(wingDesign,simParam,wingModelStructure,wingModelAero,storeAllData);
else
%     % load parameters
%     parameters = load(['data/par_all_FSI_', airfoil, '.mat']);
%     paramFSI = parameters.paramFSI;
%     wingDesign = parameters.wingDesign;
%     simParam = parameters.simParam;
%     wingModelStructure = parameters.wingModelStructure;
%     wingModelAero = parameters.wingModelAero;
    
    parameters = load(['data/parsim_FSI_', airfoil, '.mat']);
    paramFSI = parameters.paramFSI;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% run FSI test cases:
% 
%   - comparing modal vs. full FE model with steady panel method
%   - run unsteady FSI test case
%   - compare unsteady panel method with Theodorsen's function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if runTest
%     FSI_testCases = runFSItestcases(paramFSI, simParam, wingModelStructure, wingDesign, wingModelAero);
% end
if runTestSteady
    SteadyTestCases = runSteadyFSItestcases(paramFSI, simParam, wingModelStructure, wingDesign, wingModelAero);
end
if runTestUnsteady
    UnsteadyTestCases = runUnsteadyFSItestcases(paramFSI);
end
if runTestTheodorsen
    TheodorsenTestCases = runTheodorsenFSItestcases(paramFSI);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% run (parametric) reduced order models:
% 
%   - generate reduced order models
%   - compare different methods: accuracy of model as a function of model rank
%
%   - currently, MAIN_ROM works with the morphing wing. However, the
%       ROMs can also be generated for: a) aero only, or b) flexible
%       non-morphing case. In these cases, the number of inputs need to be
%       changed (from 6 to 4, as the two morphing inputs are not used).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if runROM
    MAIN_ROM(paramFSI);
end

