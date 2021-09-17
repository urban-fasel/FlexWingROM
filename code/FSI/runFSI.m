% this function simulates the FSI solver (or "FULL") for given speeds and
% input signals

function simOut = runFSI(VelocityROM, dT, iTest, simInput, progressbarText)

paramFSI = simInput.paramFSI;

%% Simulation parameters

deltaAct.DsimTime = dT;

paramFSI.NM.xi = paramFSI.inputCreate.restartIN.FSI_out.NM.xi;
paramFSI.NM.xdi = paramFSI.inputCreate.restartIN.FSI_out.NM.xdi; 
paramFSI.NM.xddi = paramFSI.inputCreate.restartIN.FSI_out.NM.xddi; 
paramFSI.NM.fi = paramFSI.inputCreate.restartIN.FSI_out.NM.fi; 

simOut.cL = zeros(iTest,1);
simOut.cD = zeros(iTest,1);
simOut.cRoll = zeros(iTest,1);
simOut.cYaw = zeros(iTest,1);
simOut.cPitch = zeros(iTest,1);
simOut.bendingModeAmplitude = zeros(iTest,1); 
simOut.Snapshot = zeros(618, iTest);


%% run full simulation

meanRotationRates = simInput.meanRotationRates;
inputSnapshot0 = simInput.inputSnapshot0;

inputMorphingSymmetric = simInput.inputMorphingSymmetric;
alphaIN = simInput.alphaIN;
force_input = simInput.force_input;
Omega = simInput.Omega;
deltaAlpha = simInput.deltaAlpha;

rotOmega = simInput.rotOmega;
AlphaShape = simInput.AlphaShape;
forceShape1 = simInput.forceShape1;
forceShape2 = simInput.forceShape2;

VelocitySWEEP = linspace(min(VelocityROM),max(VelocityROM),iTest);

saveSnapshots = false;

progressbar(progressbarText);

for i = 1:iTest
    
    % input FSI
    inputCreate = paramFSI.inputCreate;
    input=[cos(alphaIN+deltaAlpha*AlphaShape(i))*VelocitySWEEP(i),0,sin(alphaIN+deltaAlpha*AlphaShape(i))*VelocitySWEEP(i),...
        meanRotationRates+Omega*rotOmega(i),force_input*forceShape1(i),-force_input*forceShape1(i)];
    
    % run FSI
    FSI_OUT = FSI_unsteady_main(input,paramFSI,inputCreate,deltaAct,inputMorphingSymmetric*forceShape2(i),saveSnapshots,inputSnapshot0);   
    
    % overwrite input
    inputCreate.restartIN = FSI_OUT.restartIN;
    paramFSI.inputCreate = inputCreate;   
    
    paramFSI.NM = FSI_OUT.NM;
    paramFSI.NM.act = [force_input*forceShape1(i),-force_input*forceShape1(i)];
    
    % save output
    simOut.cL(i) = FSI_OUT.cL;
    simOut.cLroot(i) = FSI_OUT.cLroot;
    simOut.cD(i) = FSI_OUT.cDi + FSI_OUT.cDvisc;
    simOut.cRoll(i) = FSI_OUT.cRoll;
    simOut.cPitch(i) = FSI_OUT.cPitch;
    simOut.cYaw(i) = FSI_OUT.cYaw;
    simOut.bendingModeAmplitude(i) = FSI_OUT.NM.xi(2); 
    simOut.Snapshot(:,i) = FSI_OUT.Snapshot;
    
    progressbar(i/iTest);
    
end

progressbar(1);


end