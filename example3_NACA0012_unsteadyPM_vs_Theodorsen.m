
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
%  Running the unsteady panel method will take a while. That's why we reduced the order of the model!
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

airfoil = 'NACA_0012';      % choose airfoil: coordinates are generated with createNACA4.m 
morphingWing = false;       % set true for wing design with compliant ribs -> morphing for roll and load control
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
%% run test case:
% 
%   - compare unsteady panel method with Theodorsen's function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paramFSI.generateNewData = generateNewData;
paramFSI.plt = plt;

TheodorsenTestCases = runTheodorsenFSItestcases(paramFSI);



