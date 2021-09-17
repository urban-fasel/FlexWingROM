%% loop to generate data for empirical controllability gramian

function empirical_gramian_controllability_loop(paramFSI,V0,dT,nIteration,nInputs,InpAmp)

c_m = InpAmp;

T_l = eye(6);

for i = 0:nInputs
    simInput = setSimulationInputParams_empirical_gramians(paramFSI,V0,nIteration,c_m,i);
    progressbarName = ['Run impulse response empirical gramian controllability, input ', num2str(i)];
    simOut = runFSI(V0, dT, nIteration, simInput, progressbarName); % -> rename testingPhase
%     save(['data/simOut',num2str(i),'_egc_V',num2str(V0),'.mat'],'simOut','c_m','T_l')
    save(strcat('data',filesep,paramFSI.wingParams.airfoil,filesep,'ROM',...
                filesep,'V',num2str(V0),filesep,'simOut',num2str(i),'_egc_V',num2str(V0),'.mat'),'simOut','c_m','T_l'); 
end 

