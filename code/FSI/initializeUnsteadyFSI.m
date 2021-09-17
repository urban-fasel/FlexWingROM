function restartINsave = initializeUnsteadyFSI(simUnsteadyParam, paramFSI, simID)

% paramFSI.inputCreate.alphaWing = paramFSI.wingParams.alphaWing;
paramFSI.inputCreate.alphaWing = simUnsteadyParam.alpha;
paramFSI.firstIt = 1;
% paramFSI.aerodynamic_only = 0;

nModes = paramFSI.nModes;

V = simUnsteadyParam.V;


%% Simulation parameters

timeINIT = simUnsteadyParam.timeInit;%1.5;                 % Initialisation time
deltaAct.DsimTime = simUnsteadyParam.dT;
dT = simUnsteadyParam.dT;

iTestINIT = round(timeINIT/dT);

paramFSI.NM.xi = zeros(nModes,1);
paramFSI.NM.xdi = zeros(nModes,1); % x dot i
paramFSI.NM.xddi = zeros(nModes,1); % x dot dot i
paramFSI.NM.fi = zeros(nModes,1); % f i
paramFSI.NM.act = [0,0];

% Position to cut the wake
paramFSI.cutWakeAt = findCutPositionWake(paramFSI.wingParams,paramFSI.viscPre,dT);


%% Initialise FSI

force_input = 0;                                            % asymmetric morphing actuation in the inizialisation phase
inputMorphingSymmetric = 0;                                 % symmetric morphing actuation in the inizialisation phase
alphaIN0 = 0;                                               % aircraft angle of attack
alphaIN = alphaIN0;                                         % Angle of attack at the beginning of the simulation (initialisation phase)
Omega = [0,0,0];                                            % Rotation rates in the simulation
meanRotationRates = [0,0,0];

inputDMD0 = [force_input, inputMorphingSymmetric, meanRotationRates, alphaIN0];
    
MorphingAct = zeros(1,iTestINIT);
rotOmega = zeros(1,iTestINIT);
forceShape = zeros(1,iTestINIT);

saveSnapshots = false;

progressbar(['Run initialisation V = ', num2str(V)]);

for i = 1:iTestINIT
    
    inputCreate = paramFSI.inputCreate;
    input=[cos(alphaIN)*V,0,sin(alphaIN)*V,meanRotationRates+Omega*rotOmega(i),force_input*forceShape(i),-force_input*forceShape(i)];
    FSI_OUT = FSI_unsteady_main(input,paramFSI,inputCreate,deltaAct,MorphingAct(i)*inputMorphingSymmetric,saveSnapshots,inputDMD0);
    inputCreate.restartIN = FSI_OUT.restartIN;
    paramFSI.inputCreate = inputCreate;
        
    paramFSI.NM = FSI_OUT.NM;
    paramFSI.NM.act = [force_input*forceShape(i),-force_input*forceShape(i)];
    
    progressbar(i/iTestINIT);
    
    paramFSI.firstIt = 0;
end

progressbar(1);

% save inputCreate.restartIN for restarting the simulation without
% initialising the FSI
restartINsave = inputCreate.restartIN;
if ~exist(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,simID,filesep,'V',num2str(V)), 'dir')
   mkdir(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,simID,filesep,'V',num2str(V)))
end
save(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,simID,filesep,'V',num2str(V),filesep,'FSI_restartIN.mat'), 'restartINsave');
