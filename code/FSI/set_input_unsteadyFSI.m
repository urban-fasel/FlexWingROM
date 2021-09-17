% Input signals used in the simulation of the ROM are defined

function simInput = set_input_unsteadyFSI(paramFSI,simUnsteadyParam)

paramFSI.inputCreate.alphaWing = simUnsteadyParam.alpha;
paramFSI.firstIt = 0;
paramFSI.cutWakeAt = findCutPositionWake(paramFSI.wingParams,paramFSI.viscPre,simUnsteadyParam.dT);

simInput.paramFSI = paramFSI;

force_input = 0;                                    % amplitude asymmetric morphing actuation in the inizialisation phase
inputMorphingSymmetric = 0;                         % amplitude symmetric morphing actuation in the inizialisation phase
alphaIN0 = 0;                                       % zero aircraft angle of attack 
simInput.meanRotationRates = [0,0,0];
simInput.inputSnapshot0 = [force_input, inputMorphingSymmetric, simInput.meanRotationRates, alphaIN0]; % zero input aDMDc

simInput.inputMorphingSymmetric = 0;                % Morphing actuation input (symmetric -> load alleviation)
simInput.alphaIN = alphaIN0;                        % Angle of attack at the beginning of the simulation (initialisation phase)
simInput.force_input = 0;                           % Morphing actuation input (unsymmetric -> roll)
simInput.Omega = [0,0,0];                           % Rotation rates input
simInput.deltaAlpha = simUnsteadyParam.deltaAlpha;  % Amplitude angle of attack


k = simUnsteadyParam.k; % reduced frequency

% Alpha signal: Reduced frequency k
simInput.AlphaShape = sin(2*pi*(k*simUnsteadyParam.V/paramFSI.wingParams.chord*2/2/pi)*simUnsteadyParam.dT*(0:simUnsteadyParam.iTest-1));       

simInput.rotOmega = zeros(1,simUnsteadyParam.iTest);    % Rotation rate signal
simInput.forceShape1 = zeros(1,simUnsteadyParam.iTest); % Morphing actuation (unsymmetric -> roll)
simInput.forceShape2 = zeros(1,simUnsteadyParam.iTest); % Morphing actuation (symmetric -> load alleviation)

