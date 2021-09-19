function results = runTheodorsenFSItestcases(paramFSI)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% comparison unsteady panel method with Theodorsen
% compare for symmetric airfoil, high aspect ratio, rigid wing at alpha = 0deg -> closest
% to 2D flat plat of Theodorsen
%

simID = 'Theodorsen';      % simulation identification, used to save and load data

paramFSI.aeroOnly = 1;      % run aerodynamics only (no wing deformation): true for comparison with Theodorsen

simUnsteadyParam.V = 50;                % flight velocity
simUnsteadyParam.alpha = 0;             % angle of attack
simUnsteadyParam.deltaAlpha = 1*pi/180; % Amplitude angle of attack

simUnsteadyParam.timeInit = 1.5;        % length of the simulation initialization
simUnsteadyParam.dT = 0.006;            % timestep
simUnsteadyParam.iTest = round(simUnsteadyParam.timeInit/simUnsteadyParam.dT); % total size of the simulated trajectory  

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


%% run quasi-steady FSI for comparison
simSteadyParam.V = simUnsteadyParam.V;                      % flight velocity
simSteadyParam.alpha = simUnsteadyParam.deltaAlpha*180/pi;	% angle of attack
simSteadyParam.forceACTinp.L = 0;   % actuator force unsymmetric actuation (roll actuation), [-1 1]
simSteadyParam.forceACTinp.fs = 0;  % actuator force symmetric actuation (cL control), [-1 1]
simSteadyParam.doPlot = 0;          % plot deformed wing

% run quasi steady aero
FSI_steady_out = FSI_steady(paramFSI, simSteadyParam);


%% initialize unsteady FSI for Theodorsen

% run or load initialization of unsteady simulation: initializes the wake. FSI is then restarted from the initialized wake
initializeUnsteadySim = 0; % true will run initialization

if initializeUnsteadySim
    paramFSI.inputCreate.restartIN = [];
    restartINsave = initializeUnsteadyFSI(simUnsteadyParam, paramFSI, simID);
    paramFSI.inputCreate.restartIN = restartINsave;
else
    training_restartIN = load(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,simID,filesep,'V',num2str(simUnsteadyParam.V),filesep,'FSI_restartIN.mat')); % load wing parameters
    paramFSI.inputCreate.restartIN = training_restartIN.restartINsave;
end



%% comparison unsteady panel method with Theodorsen

runSim = 0; % run simulations or load results 

% run loop for different reduced frequencies
kRange = [0.01:0.02:0.07 0.1:0.1:0.5];
simUnsteadyParam.time = 1.5; % length of the simulation 
simUnsteadyParam.iTest = round(simUnsteadyParam.time/simUnsteadyParam.dT); % total size of the simulated trajectory  

if runSim
    for ii = 1:length(kRange)
        % set inputs 
        simUnsteadyParam.k = kRange(ii);        % Reduced frequency (angle of attack)
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
else
    for ii = 1:length(kRange)
        simUnsteadyParam.k = kRange(ii);        % Reduced frequency (angle of attack)
        simOut = load(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,'Theodorsen',filesep,'V',num2str(simUnsteadyParam.V),filesep,...
            sprintf('simOutUnsteady_V%d_alpha%d_k0_%d.mat',simUnsteadyParam.V,simUnsteadyParam.alpha,100*simUnsteadyParam.k)));
        simOutUnsteadyS{ii} = simOut.simOutUnsteady;
    end
end


%% calculate amplitude and phase lag unsteady aero
startfrom = 2;

for ii = 1:length(kRange)
    
    frequency = kRange(ii)*simUnsteadyParam.V/paramFSI.wingParams.chord*2/2/pi;
    period = 1/frequency;
    alphaPhase = sin(2*pi*(frequency)*simUnsteadyParam.dT*(0:simUnsteadyParam.iTest-1)); % alpha input signal normalized
    alphaPhase = alphaPhase/max(alphaPhase);
    
    startCheck = 1;
    
    % quasi steady
    qsx = 0:simUnsteadyParam.dT:simUnsteadyParam.time-simUnsteadyParam.dT;
    qsy = alphaPhase*max(FSI_steady_out.cLroot);
    for j = length(qsx)-startCheck:-1:startfrom
        if qsy(j-1)<qsy(j) && qsy(j)>qsy(j+1)
            xx = linspace(qsx(j-2),qsx(j+2),100);
            yy = spline([qsx(j-2),qsx(j-1),qsx(j),qsx(j+1),qsx(j+2)],[qsy(j-2),qsy(j-1),qsy(j),qsy(j+1),qsy(j+2)],xx);
            maxqsy = max(yy); % max quasi steady lift coefficient
            maxqsx = xx(find(yy == maxqsy));
            break
        end
    end
    
    % unsteady
    rx = qsx;
    ry = simOutUnsteadyS{ii}.cLroot;
    for j = length(rx)-startCheck:-1:startfrom
        if ry(j-1)<ry(j) && ry(j)>ry(j+1)
            xx = linspace(rx(j-2),rx(j+2),100);
            yy = spline([rx(j-2),rx(j-1),rx(j),rx(j+1),rx(j+2)],[ry(j-2),ry(j-1),ry(j),ry(j+1),ry(j+2)],xx);
            maxry = max(yy); % max unsteady lift coefficient
            maxrx = xx(find(yy == maxry));
            break
        end
    end
    phase(ii) = -(maxrx-maxqsx)/period * 360;
    amplitude(ii) = abs(maxry/maxqsy);
end


%% Theodorsen function

for ii = 1:length(kRange)
    % Bessel function
    C_theod(ii) = (besselh(1,2,kRange(ii)))./...
        (besselh(1,2,kRange(ii))...
        +1i*besselh(0,2,kRange(ii)));
%     phaseTheodorsen(ii) = angle(C_theod(ii))*180/pi;
%     amplitudeTheodorsen(ii) = abs(C_theod(ii));
    
    rotAxis = 0.5; % pitch around half chord
    a = -1+rotAxis*2;
    ARG = (-pi/4)*a*(- kRange(ii)^2) + (pi/2)*(1i*kRange(ii)) + 2*pi*C_theod(ii) + pi*(.5-a)*(1i*kRange(ii)).*C_theod(ii);
    amplitudeTheodorsen(ii) = abs(ARG)/(2*pi);
    phaseTheodorsen(ii) = angle(ARG)*180/pi;
end

% plot
figure
plot(kRange,phaseTheodorsen,'b'); hold on
plot(kRange,phase,'r--')
xlabel('reduced frequency k')
ylabel('phase lag, deg')
legend({'Theodorsen', 'Unsteady panel method'})

figure
plot(kRange,amplitudeTheodorsen,'b'); hold on
plot(kRange,amplitude,'r--')
xlabel('reduced frequency k')
ylabel('amplitude')
legend({'Theodorsen', 'Unsteady panel method'})


results.kRange = kRange;
results.amplitudeTheodorsen = amplitudeTheodorsen;
results.phaseTheodorsen = phaseTheodorsen;
results.amplitude = amplitude;
results.phase = phase;

