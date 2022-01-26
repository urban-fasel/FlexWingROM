
%%%%%%%%%%
%
%% This test code runs a flexible NACA6412 wing unsteady FSI and genereates an animation of the wing and wake displacements
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

generateWing = true;        % generate or load wing design parameters and FSI model
generateNewData = true;     % run FSI to generate new data, or load data
storeAllData = false;       % store all data: needed for runTestSteady

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
%% run test case:
% 
%   - run unsteady FSI test cases
%        1. Sinusoidal pitching
%        2. Sinusoidal unsymmetric morphing -> generating roll moment 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
% pitching
paramFSI.generateNewData = generateNewData;
paramFSI.initializeFSI = generateNewData;
paramFSI.plt = plt;
paramFSI.createAnimation = true;
paramFSI.animationName = 'pitching';
if paramFSI.createAnimation
    clearOldAnimationV(paramFSI)
end

simUnsteadyParam.V = 30;                % flight velocity
simUnsteadyParam.alpha = 6;             % angle of attack

simUnsteadyParam.time = 3.0;            % length of the simulation
simUnsteadyParam.dT = 0.006;            % timestep
simUnsteadyParam.iTest = round(simUnsteadyParam.time/simUnsteadyParam.dT); % total size of the simulated trajectory  
simUnsteadyParam.timeInit = 1.5;        % length of the initialization simulation

simUnsteadyParam.k = 0.2;               % Reduced frequency (angle of attack)
simUnsteadyParam.deltaAlpha = 2*pi/180; % Amplitude angle of attack
simUnsteadyParam.morphingRoll = 0;      % Amplitude morphing to generate roll moment

UnsteadyTestCasePitching = runUnsteadyFSItestcases(paramFSI,simUnsteadyParam);


%%%%%%%%%%%%%%%%
% morphing
paramFSI.generateNewData = generateNewData;
paramFSI.initializeFSI = false;         % load initialization from pitching
paramFSI.plt = plt;
paramFSI.createAnimation = true;
paramFSI.animationName = 'morphing';
if paramFSI.createAnimation
    clearOldAnimationV(paramFSI)
end

simUnsteadyParam.k = 0.1;               % Reduced frequency (morphing)
simUnsteadyParam.deltaAlpha = 0;        % Amplitude angle of attack
simUnsteadyParam.morphingRoll = 1;      % Amplitude morphing to generate roll moment

UnsteadyTestCaseRolling = runUnsteadyFSItestcases(paramFSI,simUnsteadyParam);



