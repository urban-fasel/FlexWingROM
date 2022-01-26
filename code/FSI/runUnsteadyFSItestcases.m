function results = runUnsteadyFSItestcases(paramFSI,simUnsteadyParam)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% run FSI with unsteady panel method
%

simID = 'TEST_UNSTEADY';      % simulation name, used to save and load data
paramFSI.aeroOnly = 0;      % run aerodynamics only (no wing deformation): false for FSI

% FE-solver: Newark method parameters
gammaN = 1/2;
betaN = 1/4;
paramFSI.NM.gammaN = gammaN;
paramFSI.NM.betaN = betaN;
paramFSI.NM.b1 = simUnsteadyParam.dT;
paramFSI.NM.b2 = simUnsteadyParam.dT^2*(0.5-betaN);
paramFSI.NM.b3 = (1-gammaN)*simUnsteadyParam.dT;
paramFSI.NM.b4 = gammaN*simUnsteadyParam.dT;
paramFSI.NM.b5 = betaN*simUnsteadyParam.dT^2;


%% initialize FSI

% run or load initialization of unsteady simulation: initializes the wake. FSI is then restarted from the initialized wake
initializeUnsteadySim = paramFSI.initializeFSI; % true will run initialization

if initializeUnsteadySim
    paramFSI.inputCreate.restartIN = [];
    restartINsave = initializeUnsteadyFSI(simUnsteadyParam, paramFSI, simID);
    paramFSI.inputCreate.restartIN = restartINsave;
else
    training_restartIN = load(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,simID,filesep,'V',num2str(simUnsteadyParam.V),filesep,'FSI_restartIN.mat')); % load wing parameters
    paramFSI.inputCreate.restartIN = training_restartIN.restartINsave;
end


%% run FSI

% set inputs 
simInput = set_input_unsteadyFSI(paramFSI,simUnsteadyParam); % define FSI inputs

runSim = paramFSI.generateNewData; % run simulations or load results 

% run FSI
if runSim
    progressbarText = ['Run unsteady FSI V = ',num2str(simUnsteadyParam.V), ', alpha = ', num2str(simUnsteadyParam.alpha), '.'];
    simOutUnsteady = runFSI(simUnsteadyParam.V,simUnsteadyParam.dT,simUnsteadyParam.iTest,simInput,progressbarText);
    save(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,simID,filesep,...
        sprintf('simOutUnsteady_V%d_alpha%d.mat',simUnsteadyParam.V,simUnsteadyParam.alpha)), 'simOutUnsteady');
else
    simOut = load(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,simID,filesep,...
        sprintf('simOutUnsteady_V%d_alpha%d.mat',simUnsteadyParam.V,simUnsteadyParam.alpha)));
    simOutUnsteady = simOut.simOutUnsteady;
end

results.simOutUnsteady = simOutUnsteady;

% plot
if paramFSI.plt
    plotUnsteadyFSI(simOutUnsteady,simUnsteadyParam)
end

% animation
if paramFSI.createAnimation 
    animateUnsteadyFSI(simInput,simUnsteadyParam)
end

