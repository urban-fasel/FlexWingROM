function simInput = setSimulationInputParams_AI_chirp(paramFSI,simInput,dT,iTest)

paramFSI.inputCreate.alphaWing = simInput.alpha;
VelocityROM = min(simInput.VelocityROM);
paramFSI.firstIt = 0;
paramFSI.cutWakeAt = findCutPositionWake(paramFSI.wingParams,paramFSI.viscPre,paramFSI.dT);

simInput.paramFSI = paramFSI;

force_input = 0;                                   % asymmetric morphing actuation in the inizialisation phase
inputMorphingSymmetric = 0;                        % symmetric morphing actuation in the inizialisation phase
alphaIN0 = 0;                                      % zero aircraft angle of attack (true angle of attack = alphaIN0+paramFSI.wingParams.alphaWing = 5deg)
simInput.meanRotationRates = [0,0,0];
simInput.inputSnapshot0 = [force_input, inputMorphingSymmetric, simInput.meanRotationRates, alphaIN0]; % zero input aDMDc

simInput.inputMorphingSymmetric = 0.5;      % Morphing actuation input (symmetric -> load alleviation)
simInput.alphaIN = alphaIN0;                % Angle of attack at the beginning of the simulation (initialisation phase)
simInput.force_input = 0.5;                 % Morphing actuation input (unsymmetric -> roll)
simInput.Omega = 1*[pi/10,pi/10,pi/10];     % Rotation rates input
simInput.deltaAlpha = 4*pi/180;             % Delta alpha input


f_ref=0.0025/0.006*VelocityROM/mean(paramFSI.wingParams.chord)*2/2/pi;

f0=.05*f_ref;   % Hz
f1=.4*f_ref;    % Hz
phase_chirp=0;  % to start from 0
t=dT*(0:iTest-1);
uin_chirp = mychirp_AI(t,f0,t(end),f1,phase_chirp);

uin_chirp2 = mychirp_AI(t,0.9*f0,t(end),1.1*f1,phase_chirp);

simInput.rotOmega = uin_chirp;          % Rotation rate signal
simInput.AlphaShape = uin_chirp2;       % Alpha signal
simInput.forceShape1 = uin_chirp;       % Morphing actuation (unsymmetric -> roll)
simInput.forceShape2 = uin_chirp2;    	% Morphing actuation (symmetric -> load alleviation)
simInput.t = t;                         % Morphing actuation (symmetric -> load alleviation)

end


function x=mychirp_AI(t,f0,t1,f1,phase_chirp)
%Y = mychirp(t,f0,t1,f1) generates samples of a linear swept-frequency
% signal at the time instances defined in timebase array t. The instantaneous
% frequency at time 0 is f0 Hertz. The instantaneous frequency f1
% is achieved at time t1.
% The argument 'phase' is optional. It defines the initial phase of the
% signal degined in radians. By default phase=0 radian

t0=t(1);
 T=t1-t0;
 k=(f1-f0)/T;
 x=cos(2*pi*(k/2*t+f0).*t+phase_chirp);
end


