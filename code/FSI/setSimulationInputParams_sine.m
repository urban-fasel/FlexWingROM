function simInput = setSimulationInputParams_sine(paramFSI,simInput)

paramFSI.inputCreate.alphaWing = simInput.alpha;
VelocityROM = min(simInput.VelocityROM);

paramFSI.firstIt = 0;

simInput.paramFSI = paramFSI;

force_input = 0;                                    % asymmetric morphing actuation in the inizialisation phase
inputMorphingSymmetric = 0;                         % symmetric morphing actuation in the inizialisation phase
alphaIN0 = 0;                                       % zero aircraft angle of attack (true angle of attack = alphaIN0+paramFSI.wingParams.alphaWing = 5deg)
simInput.meanRotationRates = [0,0,0];
simInput.inputSnapshot0 = [force_input, inputMorphingSymmetric, simInput.meanRotationRates, alphaIN0]; % zero input aDMDc

simInput.inputMorphingSymmetric = 1;        % Morphing actuation input (symmetric -> load alleviation)
simInput.alphaIN = alphaIN0;                % Angle of attack at the beginning of the simulation (initialisation phase)
simInput.force_input = 0.5;                 % Morphing actuation input (unsymmetric -> roll)
simInput.Omega = 1*[pi/10,0,0];             % Rotation rates input
simInput.deltaAlpha = 4*pi/180;             % Delta alpha input

simInput.rotOmega = sin(2*pi*(0.1*VelocityROM/paramFSI.wingParams.chord*2/2/pi)*paramFSI.dT*(0:paramFSI.iTest-1));          % Rotation rate signal: Reduced frequency 0.1
simInput.AlphaShape = sin(2*pi*(0.2*VelocityROM/paramFSI.wingParams.chord*2/2/pi)*paramFSI.dT*(0:paramFSI.iTest-1));        % Alpha signal: Reduced frequency 0.2
simInput.forceShape1 = sin(2*pi*(0.05*VelocityROM/paramFSI.wingParams.chord*2/2/pi)*paramFSI.dT*(0:paramFSI.iTest-1));      % Morphing actuation (unsymmetric -> roll): Reduced frequency 0.05
simInput.forceShape2 = zeros(1,paramFSI.iTest);                                                                             % Morphing actuation (symmetric -> load alleviation): Reduced frequency 0.05

