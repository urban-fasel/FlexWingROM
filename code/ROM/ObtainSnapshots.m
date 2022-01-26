function restartINsave = ObtainSnapshots(simInput, paramFSI)

paramFSI.inputCreate.alphaWing = simInput.alpha; 
VelocityROM = simInput.VelocityROM; 

paramFSI.firstIt = 1;  % simulaiton is initialized, therefore no wake information loaded
paramFSI.aeroOnly = 0; % run fluid structure interaction

nModes = paramFSI.nModes;


%% Simulation parameters

timeINIT = 1.5;                 % Initialisation time
time = paramFSI.time;           % Simulation time 
deltaAct.DsimTime = paramFSI.dT;
dT = paramFSI.dT;

iTestT = round(time/dT);
iTestINIT = round(timeINIT/dT);

paramFSI.NM.xi = zeros(nModes,1);   % x: modal amplitudes
paramFSI.NM.xdi = zeros(nModes,1);  % x dot i
paramFSI.NM.xddi = zeros(nModes,1); % x dot dot i
paramFSI.NM.fi = zeros(nModes,1);   % external forces in modal coordinates: f i
paramFSI.NM.act = [0,0];

simOut.cL = zeros(iTestT,1); % lift coefficient
simOut.cD = zeros(iTestT,1); % drag coefficient
simOut.cRoll = zeros(iTestT,1); % roll coefficient
simOut.cYaw = zeros(iTestT,1); % yaw coefficient
simOut.cPitch = zeros(iTestT,1); % pitch coefficient


%% Initialise FSI

force_input = 0;               	% amplitude asymmetric morphing actuation in the inizialisation phase: control cRoll
inputMorphingSymmetric = 0; 	% amplitude symmetric morphing actuation in the inizialisation phase: controls cL
alphaIN0 = 0;                 	% aircraft angle of attack
alphaIN = alphaIN0;         	% angle of attack at the beginning of the simulation (initialisation phase)
Omega = [0,0,0];              	% amplitude rotation rates in the simulation
meanRotationRates = [0,0,0];

inputSnapshot0 = [force_input, inputMorphingSymmetric, meanRotationRates, alphaIN0]; % initial input that is subtracetd from input snapshots
    
rotOmega = zeros(1,iTestINIT);      % rotation rates trajectory
forceShape = zeros(1,iTestINIT);    % asymmetric morphing actuation
MorphingAct = zeros(1,iTestINIT);   % symmetric morphing actuation

saveSnapshots = false; % dont save snapshots in the initialization phase

progressbar('Run initialisation phase');

for i = 1:iTestINIT
    
    % input FSI
    inputCreate = paramFSI.inputCreate;
    input=[cos(alphaIN)*VelocityROM,0,sin(alphaIN)*VelocityROM,meanRotationRates+Omega*rotOmega(i),force_input*forceShape(i),-force_input*forceShape(i)];
    
    % run FSI
    FSI_OUT = FSI_unsteady_main(input,paramFSI,inputCreate,deltaAct,MorphingAct(i)*inputMorphingSymmetric,saveSnapshots,inputSnapshot0);
    
    % overwrite inputs
    inputCreate.restartIN = FSI_OUT.restartIN;
    paramFSI.inputCreate = inputCreate;
        
    paramFSI.NM = FSI_OUT.NM;
    paramFSI.NM.act = [force_input*forceShape(i),-force_input*forceShape(i)];

    progressbar(i/iTestINIT);
    
    paramFSI.firstIt = 0;
end

progressbar(1);

% save inputCreate.restartIN for restarting the testing phase without initialising the FSI
restartINsave = inputCreate.restartIN;
if ~exist(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,'ROM',filesep,'V',num2str(VelocityROM)), 'dir')
   mkdir(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,'ROM',filesep,'V',num2str(VelocityROM)))
end
save(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,'ROM',...
    filesep,'V',num2str(VelocityROM),filesep,'training_restartIN.mat'), 'restartINsave');


%% run training phase

force_input = 1;                % amplitude asymmetric morphing actuation impulse
inputMorphingSymmetric = 1;     % amplitude symmetric morphing actuation impulse
alphaIN = alphaIN0;             % angle of attack at the beginning of the simulation (initialisation phase)
Omega = pi/6;                   % rotation rates amplitude impulse
deltaAlpha = 8*pi/180;          % amplitude angle of attack impulse

% random order of impulses in the direction of the 6 inputs
% each impulse with amplitude +-1, therefore 12 impulses
rng(2,'twister');
NumberImpulses = 12;
inputSequence_1 = 1:12;
inputSequence_1 = inputSequence_1(randperm(length(inputSequence_1)));
inputSequence2(find(mod(1:iTestT,30)==0,NumberImpulses)) = inputSequence_1;

AlphaShape = zeros(1,iTestT);
AlphaShape(inputSequence2==1)=1;
AlphaShape(inputSequence2==9)=-1;
rotOmega1 = zeros(1,iTestT);
rotOmega1(inputSequence2==2)=1;
rotOmega1(inputSequence2==10)=-1;
rotOmega2 = zeros(1,iTestT);
rotOmega2(inputSequence2==3)=1;
rotOmega2(inputSequence2==11)=-1;
rotOmega3 = zeros(1,iTestT);
rotOmega3(inputSequence2==4)=1;
rotOmega3(inputSequence2==12)=-1;
forceShape1 = zeros(1,iTestT);
forceShape1(inputSequence2==5)=1;
forceShape1(inputSequence2==7)=-1;
forceShape2 = zeros(1,iTestT);
forceShape2(inputSequence2==6)=1;
forceShape2(inputSequence2==8)=-1;


saveSnapshots = true; % save the snapshots to generate the ROMs

progressbar('Run training phase and save snapshots');

for i = 1:iTestT
    
    % input FSI
    inputCreate = paramFSI.inputCreate;
    input=[cos(alphaIN+deltaAlpha*AlphaShape(i))*VelocityROM,0,sin(alphaIN+deltaAlpha*AlphaShape(i))*VelocityROM,...
        meanRotationRates+Omega*[rotOmega1(i),rotOmega2(i),rotOmega3(i)],force_input*forceShape1(i),-force_input*forceShape1(i)];
    
    % run FSI
    FSI_OUT = FSI_unsteady_main(input,paramFSI,inputCreate,deltaAct,inputMorphingSymmetric*forceShape2(i),saveSnapshots,inputSnapshot0);
    
    % overwrite inputs
    inputCreate.restartIN = FSI_OUT.restartIN;
    paramFSI.inputCreate = inputCreate;    
    
    paramFSI.NM = FSI_OUT.NM;
    paramFSI.NM.act = [force_input*forceShape1(i),-force_input*forceShape1(i)];
    
    % save outputs
    simOut.cL(i) = FSI_OUT.cL;
    simOut.cD(i) = FSI_OUT.cDi + FSI_OUT.cDvisc;
    simOut.cRoll(i) = FSI_OUT.cRoll;
    simOut.cYaw(i) = FSI_OUT.cYaw;
    simOut.cPitch(i) = FSI_OUT.cPitch;
    simOut.bendingModeAmplitude(i) = FSI_OUT.NM.xi(2); % output 610
    simOut.Snapshot(:,i) = FSI_OUT.Snapshot;
    
    progressbar(i/iTestT);
    
end

progressbar(1);


% save outputs 
save(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,'ROM',...
    filesep,'V',num2str(VelocityROM),filesep,'simOut.mat'), 'simOut');


end