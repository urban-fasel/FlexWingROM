function simOutEG = generateDataEmpiricalObservabilityGramians(VelocityROM, dT, simInput, InpAmp, maxIteration, activate_criteria, nStates)
%% run FSI to get data for empirical observability gramian

paramFSI = simInput.paramFSI;

%% Simulation parameters

deltaAct.DsimTime = dT;

paramFSI.NM.xi = paramFSI.inputCreate.restartIN.FSI_out.NM.xi;
paramFSI.NM.xdi = paramFSI.inputCreate.restartIN.FSI_out.NM.xdi; % x dot i
paramFSI.NM.xddi = paramFSI.inputCreate.restartIN.FSI_out.NM.xddi; % x dot dot i
paramFSI.NM.fi = paramFSI.inputCreate.restartIN.FSI_out.NM.fi; % f i
paramFSI.NM.act = [0,0];

% for convergence criteria, reach close to initial state
cc = paramFSI.inputCreate.restartIN.Snapshot0;

% % Position to cut the wake
% paramFSI.cutWakeAt = findCutPositionWake(paramFSI.wingParams,dT);


%% initialise state impulse (perturbation on state initial condition)
meanRotationRates = simInput.meanRotationRates;
inputDMD0 = simInput.inputDMD0;

inputMorphingSymmetric = 0;
alphaIN = simInput.alphaIN;
force_input = 0;
Omega = 0*simInput.Omega;
deltaAlpha = 0;

rotOmega = simInput.rotOmega;
AlphaShape = simInput.AlphaShape;
forceShape1 = simInput.forceShape1;
forceShape2 = simInput.forceShape2;

Velocity = min(VelocityROM);
% VelocitySWEEP = linspace(min(VelocityROM),max(VelocityROM),iTest);

saveSnapshots = false;

i = 1;
for n = 1:nStates
  
    paramFSI = simInput.paramFSI;
    paramFSI.NM.xi = paramFSI.inputCreate.restartIN.FSI_out.NM.xi;
    paramFSI.NM.xdi = paramFSI.inputCreate.restartIN.FSI_out.NM.xdi; % x dot i
    paramFSI.NM.xddi = paramFSI.inputCreate.restartIN.FSI_out.NM.xddi; % x dot dot i
    paramFSI.NM.fi = paramFSI.inputCreate.restartIN.FSI_out.NM.fi; % f i
    paramFSI.NM.act = [0,0];

%     NM = paramFSI.NM;
%     NM.NM0 = paramFSI.NM0;
    inputCreate = paramFSI.inputCreate;
    input=[cos(alphaIN+deltaAlpha*AlphaShape(i))*Velocity,0,sin(alphaIN+deltaAlpha*AlphaShape(i))*Velocity,...
        meanRotationRates+Omega*rotOmega(i),force_input*forceShape1(i),-force_input*forceShape1(i)];
%     FSI_OUT = FSI_main_observabilityGramian(input,paramFSI,inputCreate,NM,deltaAct,inputMorphingSymmetric*forceShape2(i),saveSnapshots,inputDMD0,n,InpAmp);   
    paramFSI.obsGramian.nState = n;
    paramFSI.obsGramian.InpAmp = InpAmp;
    paramFSI.obsGramian.run = true;
    FSI_OUT = FSI_unsteady_main(input,paramFSI,inputCreate,deltaAct,inputMorphingSymmetric*forceShape2(i),saveSnapshots,inputDMD0);   
    inputCreateINIT(n).restartIN = FSI_OUT.restartIN;
    
    simOutEG(n).cL(i) = FSI_OUT.cL;
    simOutEG(n).cD(i) = FSI_OUT.cDi + FSI_OUT.cDvisc;
    simOutEG(n).cRoll(i) = FSI_OUT.cRoll;
    simOutEG(n).cPitch(i) = FSI_OUT.cPitch;
    simOutEG(n).cYaw(i) = FSI_OUT.cYaw;
    simOutEG(n).bendingModeAmplitude(i) = FSI_OUT.NM.xi(2); % output 610

    simOutEG(n).Snapshot = zeros(618, maxIteration);
    simOutEG(n).Snapshot(:,i) = FSI_OUT.Snapshot;

end


%% run FSI for each initial state perturbation

eps = 0.5*10^-2; %1*10^-3; % convergence criteria
epsA = 1*10^-5; % convergence criteria on absolute value of perturbed state: some states are close to zero and not relevant
minItA = 25; % minimum number of iterations for very small values
minIt = 5; % minimum number of iterations

paramFSI.obsGramian.run = false;

for n = 1:nStates
    paramFSI.inputCreate.restartIN = inputCreateINIT(n).restartIN;
    paramFSI.NM.xi = paramFSI.inputCreate.restartIN.FSI_out.NM.xi;
    paramFSI.NM.xdi = paramFSI.inputCreate.restartIN.FSI_out.NM.xdi; % x dot i
    paramFSI.NM.xddi = paramFSI.inputCreate.restartIN.FSI_out.NM.xddi; % x dot dot i
    paramFSI.NM.fi = paramFSI.inputCreate.restartIN.FSI_out.NM.fi; % f i
    paramFSI.NM.act = [0,0];
    i = 2;
    notConverged = true;

    while notConverged 

%         NM = paramFSI.NM;
%         NM.NM0 = paramFSI.NM0;
        inputCreate = paramFSI.inputCreate;
        input=[cos(alphaIN+deltaAlpha*AlphaShape(i))*Velocity,0,sin(alphaIN+deltaAlpha*AlphaShape(i))*Velocity,...
            meanRotationRates+Omega*rotOmega(i),force_input*forceShape1(i),-force_input*forceShape1(i)];
        FSI_OUT = FSI_unsteady_main(input,paramFSI,inputCreate,deltaAct,inputMorphingSymmetric*forceShape2(i),saveSnapshots,inputDMD0);   
        inputCreate.restartIN = FSI_OUT.restartIN;
        paramFSI.inputCreate = inputCreate;    

%         NM.xi = FSI_OUT.NM.xi;
%         NM.xdi = FSI_OUT.NM.xdi;
%         NM.xddi = FSI_OUT.NM.xddi;
%         NM.fi = FSI_OUT.NM.fi;

        paramFSI.NM = FSI_OUT.NM;
        paramFSI.NM.act = [force_input*forceShape1(i),-force_input*forceShape1(i)];

        simOutEG(n).cL(i) = FSI_OUT.cL;
        simOutEG(n).cD(i) = FSI_OUT.cDi + FSI_OUT.cDvisc;
        simOutEG(n).cRoll(i) = FSI_OUT.cRoll;
        simOutEG(n).cPitch(i) = FSI_OUT.cPitch;
        simOutEG(n).cYaw(i) = FSI_OUT.cYaw;
        simOutEG(n).bendingModeAmplitude(i) = FSI_OUT.NM.xi(2); % output 610

        simOutEG(n).Snapshot(:,i) = FSI_OUT.Snapshot;
        
        if i > 1
            if(activate_criteria)
%                 abs((simOutEG(n).Snapshot(n,i)-simOutEG(n).Snapshot(n,i-1))./simOutEG(n).Snapshot(n,i))
%                 abs((simOutEG(n).Snapshot(n,i)-cc(n))./simOutEG(n).Snapshot(n,i))
                c1 = abs((simOutEG(n).Snapshot(n,i)-simOutEG(n).Snapshot(n,i-1))./simOutEG(n).Snapshot(n,i)) < eps; % relative error consecutive time step
                c2 = abs((simOutEG(n).Snapshot(n,i)-cc(n))./simOutEG(n).Snapshot(n,i)) < eps; % relative error initial state
                if (((c1 && c2) || (abs(simOutEG(n).Snapshot(n,i)) < epsA && i > minItA)) && i > minIt)
%                 c1 = abs((simOutEG(n).Snapshot(610,i)-simOutEG(n).Snapshot(610,i-1))./simOutEG(n).Snapshot(610,i)) < eps; % relative error consecutive time step
%                 c2 = abs((simOutEG(n).Snapshot(610,i)-cc(610))./simOutEG(n).Snapshot(610,i)) < eps; % relative error initial state
%                 if (c1 && c2 && i > minIt)
                    notConverged = 0;
                    simOutEG(n).stopped_after = i; % stopped after i iterations
                    disp(['State n = ',num2str(n),' converged after ',num2str(i),' iterations']);
                end
            end
            if i == maxIteration
                notConverged = 0;
                disp(['State n = ',num2str(n)]);
                warning('!!!! not converged !!!!!!')
                simOutEG(n).stopped_after = -1; % if maxIteration reached
            end
        end
        i = i+1;
    end
end

% [tt,iii]=max(abs((simOutEG(n).Snapshot(1:608,i)-simOutEG(n).Snapshot(1:608,i-1))./simOutEG(n).Snapshot(1:608,i)))
% 
% figure
% plot(simOutEG(n).Snapshot(iii,:))

% figure
% plot(abs((simOutEG(n).Snapshot(610,1:95)-simOutEG(n).Snapshot(610,2:96))./simOutEG(n).Snapshot(610,1:95)))




