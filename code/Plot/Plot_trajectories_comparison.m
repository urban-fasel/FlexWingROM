%% Plotting comparison between 3 different ROMs

if length(V0) == 1
    aIODMD = load(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(V0),'/Plot/trajectory_fixV',num2str(V0),'_aIODMD.mat'],...
        'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','bending_mode_states_vec',...
        'lift_states_vec','drag_states_vec','roll_states_vec','yaw_states_vec','pitch_states_vec','rank_X_0');
    BMD = load(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(V0),'/Plot/trajectory_fixV',num2str(V0),'_BMD.mat'],...
        'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','bending_mode_states_vec',...
        'lift_states_vec','drag_states_vec','roll_states_vec','yaw_states_vec','pitch_states_vec','rank_X_0');
    aDMDc = load(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(V0),'/Plot/trajectory_fixV',num2str(V0),'_aDMDc.mat'],...
        'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','rank_X_0');
else
    aIODMD = load(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(min(V0)),'/Plot/trajectory_parametric_V',num2str(min(V0)),'-V',num2str(max(V0)),'_aIODMD.mat'],...
        'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','bending_mode_states_vec',...
        'lift_states_vec','drag_states_vec','roll_states_vec','yaw_states_vec','pitch_states_vec','rank_X_0');
    BMD = load(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(min(V0)),'/Plot/trajectory_parametric_V',num2str(min(V0)),'-V',num2str(max(V0)),'_BMD.mat'],...
        'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','bending_mode_states_vec',...
        'lift_states_vec','drag_states_vec','roll_states_vec','yaw_states_vec','pitch_states_vec','rank_X_0');
    aDMDc = load(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(min(V0)),'/Plot/trajectory_parametric_V',num2str(min(V0)),'-V',num2str(max(V0)),'_aDMDc.mat'],...
        'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','rank_X_0');
end


fontsize2 = 22; 
fontsize= 18;

timeP = dT:dT:time;
BMnorm = simOutFULL.bendingModeAmplitude(1); % normalize modal amplitudes with initial bending mode amplitude


%% aDMDc

% rank = rankROM.rankMax;
rank = rankROMplot.aDMDc;

figure

plot(timeP,aDMDc.bending_mode_vec(:,rank)/BMnorm,'b'); hold on
plot(timeP,simOutFULL.bendingModeAmplitude/BMnorm,'r--');
grid on;

set(gca,'Fontsize',fontsize,'TickLabelInterpreter', 'latex');
legend({'aDMDc','full'},'Interpreter','Latex','Fontsize',fontsize2);
xlabel('time','Interpreter','Latex','Fontsize',fontsize2);
ylabel('aDMDc bending mode','Interpreter','Latex','Fontsize',fontsize2);

title(sprintf('b-mode rank=%d aDMDc model',rank))


figure

plot(timeP,aDMDc.lift_vec(:,rank),'b'); hold on
plot(timeP,simOutFULL.cL,'r--');
grid on;

set(gca,'Fontsize',fontsize,'TickLabelInterpreter', 'latex');
legend({'aIODMDc','full'},'Interpreter','Latex','Fontsize',fontsize2);
xlabel('time','Interpreter','Latex','Fontsize',fontsize2);
ylabel('aDMDc cLift','Interpreter','Latex','Fontsize',fontsize2);

title(sprintf('cL rank=%d aDMDc model',rank))


%% aIODMD

rank = rankROMplot.aIODMD;

figure

plot(timeP,aIODMD.bending_mode_vec(:,rank)/BMnorm,'b'); hold on
plot(timeP,simOutFULL.bendingModeAmplitude/BMnorm,'r--');

grid on;

set(gca,'Fontsize',fontsize,'TickLabelInterpreter', 'latex');
legend({'aIODMDc','full'},'Interpreter','Latex','Fontsize',fontsize2);
xlabel('time','Interpreter','Latex','Fontsize',fontsize2);
ylabel('aIODMDc bending mode','Interpreter','Latex','Fontsize',fontsize2);

title(sprintf('b-mode rank=%d aIODMDc model',rank))


figure

plot(timeP,aIODMD.lift_vec(:,rank),'b'); hold on
plot(timeP,simOutFULL.cL,'r--');

grid on;

set(gca,'Fontsize',fontsize,'TickLabelInterpreter', 'latex');
legend({'aIODMDc','full'},'Interpreter','Latex','Fontsize',fontsize2);
xlabel('time','Interpreter','Latex','Fontsize',fontsize2);
ylabel('aIODMDc cLift','Interpreter','Latex','Fontsize',fontsize2);

title(sprintf('cL rank=%d aIODMDc model',rank))



%% BMD

% rank = rankROM.rankMaxBMD;
rank = rankROMplot.BMD;

figure

plot(timeP,BMD.bending_mode_vec(:,rank)/BMnorm,'b'); hold on
plot(timeP,simOutFULL.bendingModeAmplitude/BMnorm,'r--');

grid on;

set(gca,'Fontsize',fontsize,'TickLabelInterpreter', 'latex');
legend({'BMD','full'},'Interpreter','Latex','Fontsize',fontsize2);
xlabel('time','Interpreter','Latex','Fontsize',fontsize2);
ylabel('BMD bending mode','Interpreter','Latex','Fontsize',fontsize2);

title(sprintf('b-mode rank=%d BMD model',rank))


figure

plot(timeP,BMD.lift_vec(:,rank),'b'); hold on
plot(timeP,simOutFULL.cL,'r--');

grid on;

set(gca,'Fontsize',fontsize,'TickLabelInterpreter', 'latex');
legend({'BMD','full'},'Interpreter','Latex','Fontsize',fontsize2);
xlabel('time','Interpreter','Latex','Fontsize',fontsize2);
ylabel('BMD cLift','Interpreter','Latex','Fontsize',fontsize2);

title(sprintf('cL rank=%d BMD model',rank))

