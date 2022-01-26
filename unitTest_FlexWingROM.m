%  FlexWingROM test

clear all
close all
clc

tol = 1e-8;

addpath(genpath('code'))

% preconditions
assert(exist('MAIN','file') == 2,'Could not find MAIN script')

%% Test 1: generate NACA6418 wing design and simulaiton parameters
try 
    load data/unitTest1data
catch
    error('Could not load test 1 data')
end

try
    [wingDesign, simParam] = wingDesignAndSimParameters(1, 'NACA_6418', 1);
catch
    error('Could not run wingDesignAndSimParameters')
end

assert(abs(wingDesign.nRibsC - wingDesign_true.nRibsC) < tol)
assert(abs(simParam.addMM - simParam_true.addMM) < tol)

 
%% Test 2: generate finite element and aerodynamic model 
try 
    load data/unitTest1data
catch
    error('Could not load wingDesign and simParam data for test 2')
end

try 
    load data/unitTest2data
catch
    error('Could not load test 2 data')
end

% test generate wing structure model
try
    wingModelStructure = generateWingModelStructure(wingDesign,simParam);
catch
    error('Could not run generateWingModelStructure')
end

% test generate wing aero model
try
    wingModelAero = generateWingModelAero(wingDesign,simParam,wingModelStructure);
catch
    error('Could not run generateWingModelAero')
end

assert(sum(sum(abs(wingModelStructure.K_MODES - wingModel_true.K_MODES))) < tol)
assert(sum(sum(abs(wingModelStructure.C_MODES - wingModel_true.C_MODES))) < tol)
assert(sum(sum(abs(wingModelStructure.M_MODES - wingModel_true.M_MODES))) < tol)

assert(sum(sum(abs(wingModelAero.viscPre.cl - wingModel_true.viscPre.cl))) < tol)
assert(sum(sum(abs(wingModelAero.K_alpha_ind_ELLT - wingModel_true.K_alpha_ind_ELLT))) < tol)
assert(sum(abs(wingModelAero.fiACTMax - wingModel_true.fiACTMax)) < tol)
assert(sum(abs(wingModelAero.fiACTMin - wingModel_true.fiACTMin)) < tol)



%% To test the FSI code and reduced order modeling methods, follow the four examples

% % test the steady FSI
% example1_NACA0012_FSI_modal_vs_displacement

% % test the unsteady FSI
% example2_NACA6418_unsteadyFSI

% % compare the unsteady panel method with Theodorsen
% example3_NACA0012_unsteadyPM_vs_Theodorsen

% % test the ROM generation
% example4_NACA6418_ROM

