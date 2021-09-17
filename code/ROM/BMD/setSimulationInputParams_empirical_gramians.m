%% Empirical Gramian choosing input
% input_string (basically e_i):
% 0: all input are zero
% 1: first input delta
% 2: second input delta
% 3: third input delta
% 4: fourth input delta
% 5: fifth input delta
% 6: sixth input delta
% c_m: Amplitude of input
% T_l: orthogonal matrix, but set it to I_n (maybe -I_n, but not implemented yet)
function simInput = setSimulationInputParams_empirical_gramians(paramFSI,VelocityROM,iTest,c_m,input)

% training_restartIN = load(strcat('temp',filesep,'V',num2str(VelocityROM),filesep,'training_restartIN.mat')); % load wing parameters
training_restartIN = load(strcat('data',filesep,paramFSI.wingParams.airfoil,filesep,'ROM',...
                filesep,'V',num2str(VelocityROM),filesep,'training_restartIN.mat')); % load wing parameters
paramFSI.inputCreate.restartIN = training_restartIN.restartINsave;


paramFSI.inputCreate.restartIN.Snapshot0 = paramFSI.inputCreate.restartIN.FSI_out.Snapshot;

paramFSI.inputCreate.alphaWing = paramFSI.wingParams.alphaWing;
paramFSI.firstIt = 0;
paramFSI.cutWakeAt = findCutPositionWake(paramFSI.wingParams,paramFSI.viscPre,paramFSI.dT);

simInput.paramFSI = paramFSI;

force_input = 0;                                            % asymmetric morphing actuation in the inizialisation phase
inputMorphingSymmetric = 0;                                 % symmetric morphing actuation in the inizialisation phase
alphaIN0 = -1*pi/180;                                       % zero aircraft angle of attack (true angle of attack = alphaIN0+paramFSI.wingParams.alphaWing = 5deg)
simInput.meanRotationRates = [0,0,0];
simInput.inputDMD0 = [force_input, inputMorphingSymmetric, simInput.meanRotationRates, alphaIN0]; % zero input aDMDc

simInput.inputMorphingSymmetric = 0;        % Morphing actuation input (symmetric -> load alleviation)
simInput.alphaIN = alphaIN0;                % Angle of attack at the beginning of the simulation (initialisation phase)
simInput.force_input = 0;                 % Morphing actuation input (unsymmetric -> roll)
simInput.Omega = [0,0,0];             % Rotation rates input
simInput.deltaAlpha = 0;             % Delta alpha input

simInput.rotOmega = zeros(1,iTest);         % Rotation rate signal: Reduced frequency 0.1
simInput.AlphaShape = zeros(1,iTest);       % Alpha signal: Reduced frequency 0.2
simInput.forceShape1 = zeros(1,iTest);     % Morphing actuation (unsymmetric -> roll): Reduced frequency 0.05
simInput.forceShape2 = zeros(1,iTest);

simInput.rotOmega(1) = 1;
simInput.AlphaShape(1) = 1;
simInput.forceShape1(1) = 1;
simInput.forceShape2(1) = 1;

switch input
    case 1
        simInput.force_input = c_m;
    case 2
        simInput.inputMorphingSymmetric = c_m; 
    case 3
        simInput.Omega = c_m*[1,0,0];
    case 4
        simInput.Omega = c_m*[0,1,0];
    case 5
        simInput.Omega = c_m*[0,0,1];
    case 6
        simInput.deltaAlpha = c_m; 
    otherwise % basically '0'
        % disp('set every input to 0');
end
                                                    
