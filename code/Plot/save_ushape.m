function save_ushape(simOutFULL,simOutROM_mat, type_string,V0,simInput)

%% save results for plotting

rank_X_0 = size(simOutROM_mat,2);
lift_vec = zeros(1,rank_X_0);
drag_vec = zeros(1,rank_X_0);
roll_vec = zeros(1,rank_X_0);
yaw_vec = zeros(1,rank_X_0);
pitch_vec = zeros(1,rank_X_0);
bending_mode_vec = zeros(1,rank_X_0);
bending_mode_states_vec = zeros(1, rank_X_0);
if(strcmp(type_string,'BMD') || strcmp(type_string,'aIODMD'))
    lift_states_vec = zeros(1, rank_X_0);
    drag_states_vec = zeros(1, rank_X_0);
    roll_states_vec = zeros(1, rank_X_0);
    yaw_states_vec = zeros(1, rank_X_0);
    pitch_states_vec = zeros(1, rank_X_0);
end
for r = 1:size(simOutROM_mat,2)
    simOutROM = simOutROM_mat{r};
    lift_vec(r) = norm(simOutFULL.cL-simOutROM.cL,2)/norm(simOutFULL.cL,2)*100;
    drag_vec(r) = norm(simOutFULL.cD-simOutROM.cD,2)/norm(simOutFULL.cD,2)*100;
    roll_vec(r) = norm(simOutFULL.cRoll-simOutROM.cRoll,2)/norm(simOutFULL.cRoll,2)*100;
    yaw_vec(r) = norm(simOutFULL.cYaw-simOutROM.cYaw,2)/norm(simOutFULL.cYaw,2)*100;
    pitch_vec(r) = norm(simOutFULL.cPitch-simOutROM.cPitch,2)/norm(simOutFULL.cPitch,2)*100;
    bending_mode_vec(r) = norm(simOutFULL.bendingModeAmplitude - simOutROM.bendingModeAmplitude,2)/norm(simOutFULL.bendingModeAmplitude,2)*100;
    bending_mode_states_vec(r) = norm(simOutFULL.bendingModeAmplitude - simOutROM.bendingModeAmplitude_states,2)/norm(simOutFULL.bendingModeAmplitude,2)*100;
    if(strcmp(type_string,'BMD') || strcmp(type_string,'aIODMD'))
        lift_states_vec(r) = norm(simOutFULL.cL-simOutROM.cL_states,2)/norm(simOutFULL.cL,2)*100;
        drag_states_vec(r) = norm(simOutFULL.cD-simOutROM.cD_states,2)/norm(simOutFULL.cD,2)*100;
        roll_states_vec(r) = norm(simOutFULL.cRoll-simOutROM.cRoll_states,2)/norm(simOutFULL.cRoll,2)*100;
        yaw_states_vec(r) = norm(simOutFULL.cYaw-simOutROM.cYaw_states,2)/norm(simOutFULL.cYaw,2)*100;
        pitch_states_vec(r) = norm(simOutFULL.cPitch-simOutROM.cPitch_states,2)/norm(simOutFULL.cPitch,2)*100;
    end
end


if ~exist(strcat(pwd,filesep,'data',filesep,simInput.paramFSI.wingParams.airfoil,filesep,'ROM',filesep,'V',num2str(min(V0)),filesep,'Plot'), 'dir')
   mkdir(strcat(pwd,filesep,'data',filesep,simInput.paramFSI.wingParams.airfoil,filesep,'ROM',filesep,'V',num2str(min(V0)),filesep,'Plot'))
end

if length(V0) == 1
    if(strcmp(type_string,'aIODMD')||strcmp(type_string,'BMD'))
        save(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(V0),'/Plot/u_shape_V',num2str(V0),'_',type_string,'.mat'],...
            'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','bending_mode_states_vec',...
            'lift_states_vec','drag_states_vec','roll_states_vec','yaw_states_vec','pitch_states_vec','rank_X_0');
    else
        save(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(V0),'/Plot/u_shape_V',num2str(V0),'_',type_string,'.mat'],...
            'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','bending_mode_states_vec','rank_X_0');
    end
else
    if(strcmp(type_string,'aIODMD')||strcmp(type_string,'BMD'))
        save(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(min(V0)),'/Plot/u_shape_V',...
            num2str(min(V0)),'-V',num2str(max(V0)),'_',type_string,'.mat'],...
            'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','bending_mode_states_vec',...
            'lift_states_vec','drag_states_vec','roll_states_vec','yaw_states_vec','pitch_states_vec','rank_X_0');
    else
        save(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(min(V0)),'/Plot/u_shape_V',...
            num2str(min(V0)),'-V',num2str(max(V0)),'_',type_string,'.mat'],...
            'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','bending_mode_states_vec','rank_X_0');
    end
end

end