%% Plotting comparison between 3 different ROMs

aIODMD = load(['data/trajectory_fixV_V',num2str(V0),'_aIODMD.mat'],'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','bending_mode_states_vec',...
    'lift_states_vec','drag_states_vec','roll_states_vec','yaw_states_vec','pitch_states_vec','rank_X_0');
BMD = load(['data/trajectory_fixV_V',num2str(V0),'_BMD.mat'],'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','bending_mode_states_vec',...
    'lift_states_vec','drag_states_vec','roll_states_vec','yaw_states_vec','pitch_states_vec','rank_X_0');
aDMDc = load(['data/trajectory_fixV_V',num2str(V0),'_aDMDc.mat'],'lift_vec','drag_vec','roll_vec','yaw_vec','pitch_vec','bending_mode_vec','rank_X_0');


% font = 'Times New Roman';
% % figMA = figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
% figure('Name','Fig. 2','NumberTitle','off','units','normalized','DefaultTextFontName', font, 'DefaultAxesFontName', font,'outerposition',[0 0 1 1]);

figure

% x0=50;
% y0=10;
% width=800;
% height=600;%600;%1000;
% set(gcf,'units','points','position',[x0,y0,width,height])

fontsize2 = 22; 
fontsize= 18;

% plot(1:aIODMD_BALANCED_MIMO.rank_X_0, aIODMD_BALANCED_MIMO.bending_mode_vec,'-o'); hold on
% plot(1:aIODMD_MISO.rank_X_0, aIODMD_MISO.bending_mode_vec,'-s'); hold on
% plot(1:aDMDc.rank_X_0, aDMDc.bending_mode_vec,'-d','Color',	[0, 0.5, 0]);
plot(BMD.lift_vec(:,end),'b'); hold on
plot(aIODMD.lift_vec(:,end),'r--'); hold on
plot(aDMDc.lift_vec(:,end),'g');
plot(simOutFULL.cL,'k:');



grid on;

% ylim([0 100]);
% yticks(0:20:100)

% xlim([0 rankMax]);


set(gca,'Fontsize',fontsize,'TickLabelInterpreter', 'latex');

legend({'BMD','IOROM','aDMDc','full'},'Interpreter','Latex','Fontsize',fontsize2);

xlabel('time','Interpreter','Latex','Fontsize',fontsize2);

