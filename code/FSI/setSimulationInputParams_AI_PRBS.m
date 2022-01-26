function [simInput,iTest,time] = setSimulationInputParams_AI_PRBS(paramFSI,simInput,dT,iTest)

paramFSI.inputCreate.alphaWing = simInput.alpha;
paramFSI.firstIt = 0;

simInput.paramFSI = paramFSI;

force_input = 0;                            % asymmetric morphing actuation in the inizialisation phase
inputMorphingSymmetric = 0;                 % symmetric morphing actuation in the inizialisation phase
alphaIN0 = 0;                               % zero aircraft angle of attack
simInput.meanRotationRates = [0,0,0];
simInput.inputSnapshot0 = [force_input, inputMorphingSymmetric, simInput.meanRotationRates, alphaIN0]; % zero input aDMDc

simInput.inputMorphingSymmetric = 0.5;      % Morphing actuation input (symmetric -> load alleviation)
simInput.alphaIN = alphaIN0;                % Angle of attack at the beginning of the simulation (initialisation phase)
simInput.force_input = 0.5;                 % Morphing actuation input (unsymmetric -> roll)
simInput.Omega = 1*[pi/10,pi/10,pi/20];     % Rotation rates input
simInput.deltaAlpha = 4*pi/180;             % Delta alpha input


%% Input/Output DMD method for low-order model construction
[t_sim, T_sim, N_sim, u_temp] = my_PRBS(dT, 2^9-1, 1);
du = u_temp(:,2).';
uin_PRBS=du;
tin_PRBS=t_sim;
iTest=max(size(tin_PRBS));
time=tin_PRBS(end);

simInput.rotOmega = uin_PRBS(1:iTest);      % Rotation rate signal
simInput.AlphaShape = uin_PRBS(1:iTest);  	% Alpha signal
simInput.forceShape1 = uin_PRBS(1:iTest);   % Morphing actuation (unsymmetric -> roll)
simInput.forceShape2 = uin_PRBS(1:iTest);	% Morphing actuation (symmetric -> load alleviation)
simInput.t = tin_PRBS(1:iTest);          	


end




function [t_in, Tf, Nt, u] = my_PRBS(dt, period, num_period)
    t_in = zeros(1,period*num_period);
    for k = 1:period*num_period
        t_in(k) = dt*(k-1);
    end
    Tf = t_in(end);
    Nt = length(t_in);
    u_temp = idinput([period,1,num_period],'prbs');
    u = [t_in.' u_temp];
%     plot(autocorr(u(:,2))) % testing PRBS signal
%     plot(abs(fft(u(:,2))));
end



