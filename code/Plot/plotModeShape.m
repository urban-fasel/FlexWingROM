function plotModeShape(wingDesign, simParam, wingModelAero, q_ASET, modePlotScale, i)

xMesh_loop = wingModelAero.xMesh_loop;
yMesh_loop = wingModelAero.yMesh_loop;
zMesh_loop = wingModelAero.zMesh_loop;


%% Panel Method Grid

% TPS INTERPOLATION
TPSU_X = wingModelAero.TPSU_X*q_ASET;
TPSU_Y = wingModelAero.TPSU_Y*q_ASET;
TPSU_Z = wingModelAero.TPSU_Z*q_ASET;
TPSD_X = wingModelAero.TPSD_X*q_ASET;
TPSD_Y = wingModelAero.TPSD_Y*q_ASET;
TPSD_Z = wingModelAero.TPSD_Z*q_ASET;
    
size1U = simParam.n_seg_PM*2 + 1;
size2U = simParam.num_airfoil_nodes_panel;
xTPSsU = reshape(TPSU_X, size1U, size2U);
yTPSsU = reshape(TPSU_Y, size1U, size2U);
zTPSsU = reshape(TPSU_Z, size1U, size2U);
xTPSsD = reshape(TPSD_X, size1U, size2U);
yTPSsD = reshape(TPSD_Y, size1U, size2U);
zTPSsD = reshape(TPSD_Z, size1U, size2U);
    
yTPSU_tot = yMesh_loop(:,1:simParam.num_airfoil_nodes_panel) + rot90(yTPSsU,2); 
xTPSU_tot = xMesh_loop(:,1:simParam.num_airfoil_nodes_panel) + rot90(xTPSsU,2);
zTPSU_tot = zMesh_loop(:,1:simParam.num_airfoil_nodes_panel) + rot90(zTPSsU,2);

yTPSD_tot = yMesh_loop(:,simParam.num_airfoil_nodes_panel+1:end) + rot90(yTPSsD(:,2:end-1),2);  
xTPSD_tot = xMesh_loop(:,simParam.num_airfoil_nodes_panel+1:end) + rot90(xTPSsD(:,2:end-1),2);
zTPSD_tot = zMesh_loop(:,simParam.num_airfoil_nodes_panel+1:end) + rot90(zTPSsD(:,2:end-1),2);

xTPSs = [xTPSU_tot, xTPSD_tot];
yTPSs = [yTPSU_tot, yTPSD_tot];
zTPSs = [zTPSU_tot, zTPSD_tot];    

% Change wing reference system
aer_x = xTPSs;
aer_y = -zTPSs;
% aer_z = yTPSs;

%% plot modeshape

% only scale vibration modes for ploting, not morphing modes
if i > simParam.nmodes
    modePlotScale = 1;
end

yTPSU_totMODE = yMesh_loop(:,1:simParam.num_airfoil_nodes_panel) + rot90(yTPSsU,2)*modePlotScale; 
xTPSU_totMODE = xMesh_loop(:,1:simParam.num_airfoil_nodes_panel) + rot90(xTPSsU,2)*modePlotScale;
zTPSU_totMODE = zMesh_loop(:,1:simParam.num_airfoil_nodes_panel) + rot90(zTPSsU,2)*modePlotScale;

yTPSD_totMODE = yMesh_loop(:,simParam.num_airfoil_nodes_panel+1:end) + rot90(yTPSsD(:,2:end-1),2)*modePlotScale;  
xTPSD_totMODE = xMesh_loop(:,simParam.num_airfoil_nodes_panel+1:end) + rot90(xTPSsD(:,2:end-1),2)*modePlotScale;
zTPSD_totMODE = zMesh_loop(:,simParam.num_airfoil_nodes_panel+1:end) + rot90(zTPSsD(:,2:end-1),2)*modePlotScale;

xTPSsMODE = [xTPSU_totMODE, xTPSD_totMODE];
yTPSsMODE = [yTPSU_totMODE, yTPSD_totMODE];
zTPSsMODE = [zTPSU_totMODE, zTPSD_totMODE];  

aer_xMODE = xTPSsMODE;
aer_yMODE = -zTPSsMODE;
aer_zMODE = yTPSsMODE;

aer_z_UNDEF = [yMesh_loop(:,1:simParam.num_airfoil_nodes_panel), yMesh_loop(:,simParam.num_airfoil_nodes_panel+1:end)];

ScaleFactor = 1; % scale windows for save in better quality (semi-transparent matlab plots can't be saved as vectorgraphics)
dP = 0.001;
if i <= simParam.nmodes
    figM = figure('Position',[20 20 640*ScaleFactor 360*ScaleFactor],'Renderer','painters','name', sprintf('Vibration mode %d',i));
else
    figM = figure('Position',[20 20 640*ScaleFactor 360*ScaleFactor],'Renderer','painters','name', sprintf('Morphing actuation mode %d',i-simParam.nmodes));
end
hold on

xlabel('z');ylabel('x');zlabel('y');

axis equal

view([-20,12]);
ylim([-wingDesign.span/2-0.1,wingDesign.span/2+0.1])
zlim([-0.7,0.7])

surf(aer_xMODE(1:end-1,1:end),aer_yMODE(1:end-1,1:end),aer_zMODE(1:end-1,1:end),aer_zMODE(1:end-1,1:end)-aer_z_UNDEF(1:end-1,1:end),'EdgeAlpha',0.0,'FaceAlpha',0.8,'FaceColor','interp'); hold on

surf(aer_xMODE([end-1,end-1],:),[aer_yMODE(end-1,1:end);aer_yMODE(end-1,1:end)+dP],aer_zMODE([end-1,end-1],1:end),zeros(2,size(aer_y,2)),'FaceAlpha',0.0,'EdgeAlpha',1.0,'FaceColor','white') 
surf(aer_xMODE([1,1],:),[aer_yMODE(1,1:end);aer_yMODE(1,1:end)-dP],aer_zMODE([1,1],1:end),zeros(2,size(aer_y,2)),'FaceAlpha',0.0,'EdgeAlpha',1.0,'FaceColor','white') 

surf(aer_x([1:7,11:end-1],1:end),aer_y([1:7,11:end-1],1:end),aer_z_UNDEF([1:7,11:end-1],1:end),zeros(size(aer_y,1)-4,size(aer_y,2)),'FaceAlpha',1.0,'EdgeAlpha',0.0,'FaceColor',[0.6 0.6 0.6]) 

surf(aer_x([end-1,end-1],:),[aer_y(end-1,1:end);aer_y(end-1,1:end)+dP],aer_z_UNDEF([end-1,end-1],1:end),zeros(2,size(aer_y,2)),'FaceAlpha',0.0,'EdgeAlpha',1.0,'FaceColor','white') 
surf(aer_x([1,1],:),[aer_y(1,1:end);aer_y(1,1:end)-dP],aer_z_UNDEF([1,1],1:end),zeros(2,size(aer_y,2)),'FaceAlpha',0.0,'EdgeAlpha',1.0,'FaceColor','white') 

for ii = 1:16
    surf(aer_xMODE([ii,ii],:),[aer_yMODE(ii,1:end);aer_yMODE(ii,1:end)-dP],aer_zMODE([ii,ii],1:end),zeros(2,size(aer_y,2)),'FaceAlpha',0.0,'EdgeAlpha',1.0,'FaceColor','white') 
end

axis off

% % saveas(figM,sprintf('%s\\PlotModes\\mode%d_fig3D',pwd,modeNR),'epsc');
% saveas(figM,sprintf('%s\\PlotModes\\mode%d_fig3D',pwd,modeNR),'meta');
% 
% view([-90,90]);
% % saveas(figM,sprintf('%s\\PlotModes\\mode%d_figTop',pwd,modeNR),'epsc');
% saveas(figM,sprintf('%s\\PlotModes\\mode%d_figTop',pwd,modeNR),'meta');
% 
% view([-90,0]);
% % saveas(figM,sprintf('%s\\PlotModes\\mode%d_figFront',pwd,modeNR),'epsc');
% saveas(figM,sprintf('%s\\PlotModes\\mode%d_figFront',pwd,modeNR),'meta');
% 
% close all

end

