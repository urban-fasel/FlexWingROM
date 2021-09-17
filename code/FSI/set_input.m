% Input signals used in the simulation of the ROM are defined

function [simInput,iTest,time] = set_input(paramFSI,dT,iTest,time,id_input,simInput)

switch id_input
   
    case 1
        simInput = setSimulationInputParams_sine(paramFSI,simInput);

    case 2
        simInput = setSimulationInputParams_AI_chirp(paramFSI,simInput,dT,iTest);
        
    case 3
        simInput = setSimulationInputParams_AI_PRBS(paramFSI,simInput,dT,iTest);
 
end

end