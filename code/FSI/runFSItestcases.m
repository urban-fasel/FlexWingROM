function results = runFSItestcases(paramFSI, simParam, wingModelStructure, wingDesign, wingModelAero)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% compare modal (K_MODES) vs. full (K_ASET) model, running FSI with steady panel method
%

% run steady FSI test case: comparison of full stiffness matrix with modal approach
simSteadyParam.V = 30;              % flight velocity
simSteadyParam.alpha = 6;           % angle of attack
simSteadyParam.forceACTinp.L = 1;   % actuator force unsymmetric actuation (roll actuation), [-1 1]
simSteadyParam.forceACTinp.fs = 0;  % actuator force symmetric actuation (cL control), [-1 1]

simSteadyParam.doPlot = 1; % plot deformed wing

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% run FSI with unsteady panel method
%

simID = 'test_sim_us';      % simulation name, used to save and load data

paramFSI.aeroOnly = 0;      % run aerodynamics only (no wing deformation): true for comparison with Theodorsen

simUnsteadyParam.V = 30;    % flight velocity
simUnsteadyParam.alpha = 6; % angle of attack

simUnsteadyParam.time = 3.0; % length of the simulation
simUnsteadyParam.dT = 0.006; % timestep
simUnsteadyParam.iTest = round(simUnsteadyParam.time/simUnsteadyParam.dT); % total size of the simulated trajectory  

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
initializeUnsteadySim = 1; % true will run initialization

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
simUnsteadyParam.k = 0.2;               % Reduced frequency (angle of attack)
simUnsteadyParam.deltaAlpha = 2*pi/180; % Amplitude angle of attack
simInput = set_input_unsteadyFSI(paramFSI,simUnsteadyParam); % define FSI inputs

% run FSI
progressbarText = ['Run unsteady FSI V = ',num2str(simUnsteadyParam.V), ', alpha = ', num2str(simUnsteadyParam.alpha), '.'];
simOutUnsteady = runFSI(simUnsteadyParam.V,simUnsteadyParam.dT,simUnsteadyParam.iTest,simInput,progressbarText);
save(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,simID,filesep,sprintf('simOutUnsteady_V%d_alpha%d.mat',simUnsteadyParam.V,simUnsteadyParam.alpha)), 'simOutUnsteady');

results.simOutUnsteady = simOutUnsteady;

% plot
if plt
    plotUnsteadyFSI(simOutUnsteady,simUnsteadyParam)
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% comparison unsteady panel method with Theodorsen
% compare for symmetric airfoil, high aspect ratio, rigid wing at alpha = 0deg -> closest
% to 2D flat plat of Theodorsen
%

if paramFSI.wingParams.airfoil == 'NACA0012'
    
    simID = 'Theodorsen';      % simulation identification, used to save and load data

    paramFSI.aeroOnly = 1;      % run aerodynamics only (no wing deformation): true for comparison with Theodorsen

    simUnsteadyParam.V = 30;                % flight velocity
    simUnsteadyParam.alpha = 0;             % angle of attack
    simUnsteadyParam.deltaAlpha = 2*pi/180; % Amplitude angle of attack

    simUnsteadyParam.time = 3.0; % length of the simulation
    simUnsteadyParam.dT = 0.006; % timestep
    simUnsteadyParam.iTest = round(simUnsteadyParam.time/simUnsteadyParam.dT); % total size of the simulated trajectory  

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


    %% run steady for comparison
    simSteadyParam.V = simUnsteadyParam.V;              % flight velocity
    simSteadyParam.alpha = simUnsteadyParam.deltaAlpha;           % angle of attack
    simSteadyParam.forceACTinp.L = 0;   % actuator force unsymmetric actuation (roll actuation), [-1 1]
    simSteadyParam.forceACTinp.fs = 0;  % actuator force symmetric actuation (cL control), [-1 1]
    simSteadyParam.doPlot = 0; % plot deformed wing

    % run FSI using modal stiffness matrix
    FSI_steady_out_theo = FSI_steady(paramFSI, simSteadyParam);


    %% initialize unsteady FSI for Theodorsen

    % run or load initialization of unsteady simulation: initializes the wake. FSI is then restarted from the initialized wake
    initializeUnsteadySim = 1; % true will run initialization

    if initializeUnsteadySim
        paramFSI.inputCreate.restartIN = [];
        restartINsave = initializeUnsteadyFSI(simUnsteadyParam, paramFSI, simID);
        paramFSI.inputCreate.restartIN = restartINsave;
    else
        training_restartIN = load(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,simID,filesep,'V',num2str(simUnsteadyParam.V),filesep,'FSI_restartIN.mat')); % load wing parameters
        paramFSI.inputCreate.restartIN = training_restartIN.restartINsave;
    end



    %% comparison unsteady panel method with Theodorsen

    % run loop for different reduced frequencies
    kRange = [0.01:0.02:0.07 0.1:0.1:0.5];

    for ii = 1:length(kRange)
        % set inputs 
        simUnsteadyParam.k = kRange(ii);        % Reduced frequency (angle of attack)
    %     simUnsteadyParam.deltaAlpha = 2*pi/180; % Amplitude angle of attack
        simInput = set_input_unsteadyFSI(paramFSI,simUnsteadyParam);

        % run FSI
        progressbarText = ['Run unsteady PM V = ',num2str(simUnsteadyParam.V), ', alpha = ', num2str(simUnsteadyParam.alpha), ...
            ', k = ', num2str(simUnsteadyParam.k),'.'];
        simOutUnsteady = runFSI(simUnsteadyParam.V,simUnsteadyParam.dT,simUnsteadyParam.iTest,simInput,progressbarText);

        % save results
        simOutUnsteadyS{ii} = simOutUnsteady;

        if ~exist(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,'Theodorsen',filesep,'V',num2str(simUnsteadyParam.V)), 'dir')
           mkdir(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,'Theodorsen',filesep,'V',num2str(simUnsteadyParam.V)))
        end
        save(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,'Theodorsen',filesep,'V',num2str(simUnsteadyParam.V),filesep,...
            sprintf('simOutUnsteady_V%d_alpha%d_k0_%d.mat',simUnsteadyParam.V,simUnsteadyParam.alpha,100*simUnsteadyParam.k)), 'simOutUnsteady');
    end

    % figure
    % plot(simOutUnsteadyS{2}.cL); hold on
    % plot(simOutUnsteadyS{3}.cL); hold on
    % plot(simOutUnsteadyS{4}.cL); hold on
    % plot(simOutUnsteadyS{5}.cL); hold on
    % 
    % iii = 2;
    % k = kRange(iii);
    % AlphaShape = sin(2*pi*(k*simUnsteadyParam.V/paramFSI.wingParams.chord*2/2/pi)*simUnsteadyParam.dT*(0:simUnsteadyParam.iTest-1));       % Alpha signal: Reduced frequency k
    % 
    % xx = linspace(rx(j-1),rx(j+1),100);
    % yy = spline([rx(j-1),rx(j),rx(j+1)],[ry(j-1),ry(j),ry(j+1)],xx);
    % maxry = max(yy);
    % maxrx = xx(find(yy == maxry));
    % 
    % figure
    % plot(simOutUnsteadyS{iii}.cL/max(simOutUnsteadyS{iii}.cL),'b'); hold on
    % plot(AlphaShape/max(AlphaShape),'r'); hold on
    % 
    % % Theodorsen
    % 
    % C_theod= (besselh(1,2,reduced_frequency_perturbation))./...
    %     (besselh(1,2,reduced_frequency_perturbation)...
    %     +1i*besselh(0,2,reduced_frequency_perturbation));
    % first_derivative =  [zeros(1,iTestINIT),deltaVin*(2*pi*test.frequency)*(cos(2*pi*dT*(0:iTestT-1)*test.frequency))];
    % Lift_theo = (2*pi/2*1.225*0.2^2*span_of_wing*first_derivative + 2*pi*1.225*50*50*0.2*span_of_wing*...
    %     atan(deltaVin/50)*abs(C_theod).*sin([zeros(1,iTestINIT),2*pi*dT*(0:iTestT-1)*test.frequency+angle(C_theod)]))/(0.5*1.225*50*50*0.4*span_of_wing);
    % 
    % 
    % % W = logspace(-3,8,1000); % frequency range
    % % p = 0; % pitch point: 0 = around leading edge, 0.5 = around center, 1 = around trailing edge
    % % [magnitude,phase] = BodeTheoPitch(W,p); % pitch point around leading edge
end

