function empirical_gramian_observability_loop(paramFSI,V0,dT,nIteration,InpAmp)
%% loop to generate data for empirical observability gramian

c_m = pi/10;
input = 0;
simInput = setSimulationInputParams_empirical_gramians(paramFSI,V0,nIteration,c_m,input);

activate_criteria = true;
nStates = 618;
simOutEG = generateDataEmpiricalObservabilityGramians(V0, dT, simInput, InpAmp, nIteration, activate_criteria, nStates);

% save(['data/simOut_ego_V',num2str(V0),'.mat'],'simOutEG','InpAmp')
save(strcat('data',filesep,paramFSI.wingParams.airfoil,filesep,'ROM',...
            filesep,'V',num2str(V0),filesep,'simOut_ego_V',num2str(V0),'.mat'),'simOutEG','InpAmp'); 

    
