
%%%%%%%%%%
%
%% This test code runs a NACA6412 morphing wing unsteady FSI case and compares three different data driven reduced order models
% 
%  The code loads all the data from previous runs.
%
%  To generate the data (run the FSI training and testing), go to MAIN_ROM and set:
%       line 4:    runTrainingPhase = 1;       
%       line 5:    runTestingFullModel = 1;  
%
%  Running the FSI will take a while. That's why we reduced the order of the model!
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

generateWing = false;       % generate or load wing design parameters and FSI model
generateNewData = false;    % run FSI to generate new data, or load data
storeAllData = false;       % store all data: needed for runTestSteady
runROM = true;              % generate reduced order models and compare different methods
 
airfoil = 'NACA_6418';      % choose airfoil: coordinates are generated with createNACA4.m 
morphingWing = true;        % set true for wing design with compliant ribs -> morphing for roll and load control
plt = true;                 % plot results of test cases

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

paramFSI.generateNewData = generateNewData;
paramFSI.plt = plt;

if runROM
    ROM(paramFSI);
end

