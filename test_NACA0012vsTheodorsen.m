
%%%%%%%%%%
%
%% This test code runs a NACA0012 unsteady panel method and compares it with Theodorsens function
% 
%  The code loads all the data from previous runs.
%
%  To generate the data (run the unsteady panel method), go to runTheodorsenFSItestcases and set:
%       line 48:    initializeUnsteadySim = 1; % true will run initialization
%       line 63:    runSim = 1; % run simulations or load results 
%
%  Running the unsteady panel method will take a while. That's why we do reduced-order modeling!
%
%
%  Urban Fasel, 2021
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
 
airfoil = 'NACA0012';       % choose airfoil: coordinates are loaded from file in folder 'airfoils' 
morphingWing = false;       % set true for wing design with compliant ribs -> morphing for roll and load control
plt = false;                % plot

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

