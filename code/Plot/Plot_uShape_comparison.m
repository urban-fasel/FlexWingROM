%% Plotting comparison between 3 different ROMs

if length(V0) == 1
    aIODMD = load(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(V0),'/Plot/u_shape_V',num2str(V0),...
        '_aIODMD.mat'],'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','bending_mode_states_vec','rank_X_0');
    BMD = load(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(V0),'/Plot/u_shape_V',num2str(V0),...
        '_BMD.mat'],'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','bending_mode_states_vec',...
        'lift_states_vec','drag_states_vec','roll_states_vec','yaw_states_vec','pitch_states_vec','rank_X_0');
    aDMDc = load(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(V0),'/Plot/u_shape_V',num2str(V0),...
        '_aDMDc.mat'],'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','rank_X_0');
else
    aIODMD = load(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(min(V0)),'/Plot/u_shape_V',...
            num2str(min(V0)),'-V',num2str(max(V0)),...
            '_aIODMD.mat'],'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','bending_mode_states_vec','rank_X_0');
    BMD = load(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(min(V0)),'/Plot/u_shape_V',...
            num2str(min(V0)),'-V',num2str(max(V0)),...
            '_BMD.mat'],'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','bending_mode_states_vec',...
            'lift_states_vec','drag_states_vec','roll_states_vec','yaw_states_vec','pitch_states_vec','rank_X_0');
    aDMDc = load(['data/',simInput.paramFSI.wingParams.airfoil,'/ROM/V',num2str(min(V0)),'/Plot/u_shape_V',...
            num2str(min(V0)),'-V',num2str(max(V0)),...
            '_aDMDc.mat'],'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','rank_X_0');
end


fontsize2 = 22; 
fontsize= 18;



%% all outputs

nOutputs = 6;

figure

subplot(nOutputs,1,1)
plot(1:BMD.rank_X_0, BMD.lift_vec,'-o'); hold on
plot(1:aIODMD.rank_X_0, aIODMD.lift_vec,'-s'); hold on
plot(1:aDMDc.rank_X_0, aDMDc.lift_vec,'-d','Color',	[0, 0.5, 0]);
ylabel('cL','Interpreter','Latex','Fontsize',fontsize2);

grid on;
ylim([0 100]);
yticks(0:25:100)
set(gca, 'XTickLabel', [])
set(gca,'Fontsize',fontsize,'TickLabelInterpreter', 'latex');


subplot(nOutputs,1,2)
plot(1:BMD.rank_X_0, BMD.drag_vec,'-o'); hold on
plot(1:aIODMD.rank_X_0, aIODMD.drag_vec,'-s'); hold on
plot(1:aDMDc.rank_X_0, aDMDc.drag_vec,'-d','Color',	[0, 0.5, 0]);
ylabel('cD','Interpreter','Latex','Fontsize',fontsize2);

grid on;
ylim([0 100]);
yticks(0:25:100)
set(gca, 'XTickLabel', [])
set(gca,'Fontsize',fontsize,'TickLabelInterpreter', 'latex');


subplot(nOutputs,1,3)
plot(1:BMD.rank_X_0, BMD.roll_vec,'-o'); hold on
plot(1:aIODMD.rank_X_0, aIODMD.roll_vec,'-s'); hold on
plot(1:aDMDc.rank_X_0, aDMDc.roll_vec,'-d','Color',	[0, 0.5, 0]);
ylabel('cRoll','Interpreter','Latex','Fontsize',fontsize2);

grid on;
ylim([0 100]);
yticks(0:25:100)
set(gca, 'XTickLabel', [])
set(gca,'Fontsize',fontsize,'TickLabelInterpreter', 'latex');


subplot(nOutputs,1,4)
plot(1:BMD.rank_X_0, BMD.pitch_vec,'-o'); hold on
plot(1:aIODMD.rank_X_0, aIODMD.pitch_vec,'-s'); hold on
plot(1:aDMDc.rank_X_0, aDMDc.pitch_vec,'-d','Color',	[0, 0.5, 0]);
ylabel('cPitch','Interpreter','Latex','Fontsize',fontsize2);

grid on;
ylim([0 100]);
yticks(0:25:100)
set(gca, 'XTickLabel', [])
set(gca,'Fontsize',fontsize,'TickLabelInterpreter', 'latex');


subplot(nOutputs,1,5)
plot(1:BMD.rank_X_0, BMD.yaw_vec,'-o'); hold on
plot(1:aIODMD.rank_X_0, aIODMD.yaw_vec,'-s'); hold on
plot(1:aDMDc.rank_X_0, aDMDc.yaw_vec,'-d','Color',	[0, 0.5, 0]);
ylabel('cYaw','Interpreter','Latex','Fontsize',fontsize2);

grid on;
ylim([0 100]);
yticks(0:25:100)
set(gca, 'XTickLabel', [])
set(gca,'Fontsize',fontsize,'TickLabelInterpreter', 'latex');


subplot(nOutputs,1,6)
plot(1:BMD.rank_X_0, BMD.bending_mode_vec,'-o'); hold on
plot(1:aIODMD.rank_X_0, aIODMD.bending_mode_vec,'-s'); hold on
plot(1:aDMDc.rank_X_0, aDMDc.bending_mode_vec,'-d','Color',	[0, 0.5, 0]);
ylabel('b-mode','Interpreter','Latex','Fontsize',fontsize2);

grid on;
ylim([0 100]);
yticks(0:25:100)
% xlim([0 rankMax]);
set(gca,'Fontsize',fontsize,'TickLabelInterpreter', 'latex');

legend({'BMD','IOROM','aDMDc'},'Interpreter','Latex','Fontsize',fontsize2);
xlabel('Model rank','Interpreter','Latex','Fontsize',fontsize2);


