function MAIN_ROM(paramFSI)

% Choose whether the data is generated running the FSI or data is loaded
runTrainingPhase = paramFSI.generateNewData;       % initialise wake, run training phase (impulse response in direction of all inputs), and store snapshots
runTestingFullModel = paramFSI.generateNewData;    % run testing phase with full FSI model, used for comparison with reduced order models

% Preliminary parameters definition
time = 3.0;     % length of the simulation: both training and testing phase
dT = 0.006;     % timestep
iTest = round(time/dT); % total size of the simulated trajectory  

% FE solver Newark method parameters
gammaN = 1/2;
betaN = 1/4;
paramFSI.NM.gammaN = gammaN;
paramFSI.NM.betaN = betaN;
paramFSI.NM.b1 = dT;
paramFSI.NM.b2 = dT^2*(0.5-betaN);
paramFSI.NM.b3 = (1-gammaN)*dT;
paramFSI.NM.b4 = gammaN*dT;
paramFSI.NM.b5 = betaN*dT^2;

paramFSI.dT = dT;
paramFSI.time = time;
paramFSI.iTest = iTest;


%% Simulation parameters:

% Velocity
% if only one velocity -> run non-parametric ROM
% if array -> run parametric ROM
% V0 = 30;
V0 = 30:2:40;
n_g = length(V0);

% wing base angle of attack
simInput.alpha = 6;
        
% Choose the type of input used in simulation
id_input = 1; 
% 1 --> sinusoid
% 2 --> chirp
% 3 --> PRBS
id_inputS = {'sinusoid', 'chirp', 'PRBS'}; 


%% training phase
% initialise wake and train model: random order of impulses 
% save snapshots and restart data in data folder
for k = 1:n_g
    if runTrainingPhase
        clearOldROMsV(paramFSI,V0(k))
        paramFSI.inputCreate.restartIN = [];
        simInput.VelocityROM = V0(k);
        restartINsave = ObtainSnapshots(simInput, paramFSI);
        paramFSI.inputCreate.restartIN = restartINsave;
    else
        if k == 1 % only load lowest velocity: that's the velocity at which the testing phase initiates
            training_restartIN = load(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,'ROM',...
                filesep,'V',num2str(V0(k)),filesep,'training_restartIN.mat')); % load wing parameters
            paramFSI.inputCreate.restartIN = training_restartIN.restartINsave;
        end
    end
end

%% testing phase
% set inputs 
simInput.VelocityROM = V0;
simInput = set_input(paramFSI,dT,iTest,time,id_input,simInput);

% run testing or load data
if runTestingFullModel
    progressbarText = ['Run testing phase with ',id_inputS{id_input}, ' input'];
    simOutFULL = runFSI(V0,dT,iTest,simInput,progressbarText);
    if n_g == 1
        save(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,'ROM',...
            filesep,sprintf('simOutFULL_input%d_V%d.mat',id_input,V0)), 'simOutFULL');
    else
        save(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,'ROM',...
            filesep,sprintf('simOutFULL_input%d_V%d-V%d.mat',id_input,min(V0),max(V0))), 'simOutFULL');
    end
else
    if n_g == 1
        simOutF = load(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,'ROM',...
            filesep,sprintf('simOutFULL_input%d_V%d.mat',id_input,V0)));
    else
        simOutF = load(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,'ROM',...
            filesep,sprintf('simOutFULL_input%d_V%d-V%d.mat',id_input,min(V0),max(V0))), 'simOutFULL');
    end
    simOutFULL = simOutF.simOutFULL;
end


%% Run the ROM generation + simulation + comparison with FSI solver routine for different algorithms

rankROM.rankMax = 50;               % maximum rank to be tested for aDMDc and aIODMD models
rankROM.rankDMDcTruncation1 = 40;   % truncation of first SVD in DMDc: does not affect the rank of the final system but has a large influence on accuracy of the model
rankROM.rankMaxBMD = 30;            % maximum rank for BMD model


%% aDMDc
type_string = 'aDMDc';
ROMsim(V0,rankROM,type_string,dT,iTest,simInput,simOutFULL);


%% aIODMD
type_string = 'aIODMD'; 
ROMsim(V0,rankROM,type_string,dT,iTest,simInput,simOutFULL);


%% BMD
type_string = 'BMD';
ROMsim(V0,rankROM,type_string,dT,iTest,simInput,simOutFULL);


%% Plot results

% plot model accuracy over model rank for different methods
Plot_uShape_comparison;

% ROM with minimum lift coefficient error
[minCLaDMDc,iminCLaDMDc]=min(aDMDc.lift_vec);
[minCLaIODMD,iminCLaIODMD]=min(aIODMD.lift_vec);
[minCLBMD,iminCLBMD]=min(BMD.lift_vec);

% ROM with minimum bending mode amplitude error
[minBMaDMDc,iminBMaDMDc]=min(aDMDc.bending_mode_vec);
[minBMaIODMD,iminBMaIODMD]=min(aIODMD.bending_mode_vec);
[minBMBMD,iminBMBMD]=min(BMD.bending_mode_vec);

% plot trajectories for one model: choose which rank is plotted inside Plot_trajectories_comparison
rankROMplot.aDMDc = iminCLaDMDc;
rankROMplot.aIODMD = iminCLaIODMD;
rankROMplot.BMD = iminCLBMD;
Plot_trajectories_comparison;

