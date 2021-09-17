function save_trajectory(simOutROM_mat, type_string,V0,simInput)

%% save results for plotting
rank_X_0 = size(simOutROM_mat,2);
size_traj_0 = size(simOutROM_mat{1}.cL,1);
lift_vec = zeros(size_traj_0,rank_X_0);
drag_vec = zeros(size_traj_0,rank_X_0);
roll_vec = zeros(size_traj_0,rank_X_0);
yaw_vec = zeros(size_traj_0,rank_X_0);
pitch_vec = zeros(size_traj_0,rank_X_0);
bending_mode_vec = zeros(size_traj_0,rank_X_0);
bending_mode_states_vec = zeros(size_traj_0, rank_X_0);
if(strcmp(type_string,'BMD') || strcmp(type_string,'aIODMD'))
    lift_states_vec = zeros(size_traj_0, rank_X_0);
    drag_states_vec = zeros(size_traj_0, rank_X_0);
    roll_states_vec = zeros(size_traj_0, rank_X_0);
    yaw_states_vec = zeros(size_traj_0, rank_X_0);
    pitch_states_vec = zeros(size_traj_0, rank_X_0);
end
for r = 1:size(simOutROM_mat,2)
    simOutROM = simOutROM_mat{r};
    lift_vec(:,r) = simOutROM.cL;
    drag_vec(:,r) = simOutROM.cD;
    roll_vec(:,r) = simOutROM.cRoll;
    yaw_vec(:,r) = simOutROM.cYaw;
    pitch_vec(:,r) = simOutROM.cPitch;
    bending_mode_vec(:,r) = simOutROM.bendingModeAmplitude;
    bending_mode_states_vec(:,r) = simOutROM.bendingModeAmplitude_states;
    if(strcmp(type_string,'BMD') || strcmp(type_string,'aIODMD'))
        lift_states_vec(:,r) = simOutROM.cL_states;
        drag_states_vec(:,r) = simOutROM.cD_states;
        roll_states_vec(:,r) = simOutROM.cRoll_states;
        yaw_states_vec(:,r) = simOutROM.cYaw_states;
        pitch_states_vec(:,r) = simOutROM.cPitch_states;
    end
end


if ~exist(strcat(pwd,filesep,'data',filesep,simInput.paramFSI.wingParams.airfoil,filesep,'ROM',filesep,'V',num2str(min(V0)),filesep,'Plot'), 'dir')
   mkdir(strcat(pwd,filesep,'data',filesep,simInput.paramFSI.wingParams.airfoil,filesep,'ROM',filesep,'V',num2str(min(V0)),filesep,'Plot'))
end

if length(V0) == 1
    if(strcmp(type_string,'aIODMD')||strcmp(type_string,'BMD'))
        save(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(V0),'/Plot/trajectory_fixV',num2str(V0),'_',type_string,'.mat'],...
            'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','bending_mode_states_vec',...
            'lift_states_vec','drag_states_vec','roll_states_vec','yaw_states_vec','pitch_states_vec','rank_X_0');
    else
        save(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(V0),'/Plot/trajectory_fixV',num2str(V0),'_',type_string,'.mat'],...
            'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','bending_mode_states_vec','rank_X_0');
    end
else
    if(strcmp(type_string,'aIODMD')||strcmp(type_string,'BMD'))
        save(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(min(V0)),'/Plot/trajectory_parametric_V',num2str(min(V0)),'-V',num2str(max(V0)),'_',type_string,'.mat'],...
            'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','bending_mode_states_vec',...
            'lift_states_vec','drag_states_vec','roll_states_vec','yaw_states_vec','pitch_states_vec','rank_X_0');
    else
        save(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(min(V0)),'/Plot/trajectory_parametric_V',num2str(min(V0)),'-V',num2str(max(V0)),'_',type_string,'.mat'],...
            'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','bending_mode_states_vec','rank_X_0');
    end
end


end