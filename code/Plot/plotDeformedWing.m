function plotDeformedWing(b,simSteadyParam,yMesh_loop,analysisParams,aer_x,aer_y,aer_z)

ScaleFactor = 1;%10; % scale windows to save in better quality (semi-transparent matlab plots can't be saved as vectorgraphics)
dP = 0.001;

figM = figure('Position',[20 20 640*ScaleFactor 360*ScaleFactor],'Renderer','painters','name', sprintf('FSI wing deformation, V = %dm/s, alpha = %d deg',simSteadyParam.V,simSteadyParam.alpha));

hold on

xlabel('z');ylabel('x');zlabel('y');

axis equal

view([-20,12]);
ylim([-b/2-0.1,b/2+0.1])
zl = max([0.7, 1.05*max(max(abs(aer_z)))]);
zlim([-zl,zl])

aer_z_UNDEF = [yMesh_loop(:,1:analysisParams.num_airfoil_nodes_panel), yMesh_loop(:,analysisParams.num_airfoil_nodes_panel+1:end)];

surf(aer_x(1:end-1,1:end),aer_y(1:end-1,1:end),aer_z(1:end-1,1:end),aer_z(1:end-1,1:end)-aer_z_UNDEF(1:end-1,1:end),'EdgeAlpha',0.0,'FaceAlpha',0.8,'FaceColor','interp'); hold on

surf(aer_x([end-1,end-1],:),[aer_y(end-1,1:end);aer_y(end-1,1:end)+dP],aer_z([end-1,end-1],1:end),zeros(2,size(aer_y,2)),'FaceAlpha',0.0,'EdgeAlpha',1.0,'FaceColor','white') 
surf(aer_x([1,1],:),[aer_y(1,1:end);aer_y(1,1:end)-dP],aer_z([1,1],1:end),zeros(2,size(aer_y,2)),'FaceAlpha',0.0,'EdgeAlpha',1.0,'FaceColor','white') 

surf(aer_x([1:7,11:end-1],1:end),aer_y([1:7,11:end-1],1:end),aer_z_UNDEF([1:7,11:end-1],1:end),zeros(size(aer_y,1)-4,size(aer_y,2)),'FaceAlpha',1.0,'EdgeAlpha',0.0,'FaceColor',[0.6 0.6 0.6]) 

surf(aer_x([end-1,end-1],:),[aer_y(end-1,1:end);aer_y(end-1,1:end)+dP],aer_z_UNDEF([end-1,end-1],1:end),zeros(2,size(aer_y,2)),'FaceAlpha',0.0,'EdgeAlpha',1.0,'FaceColor','white') 
surf(aer_x([1,1],:),[aer_y(1,1:end);aer_y(1,1:end)-dP],aer_z_UNDEF([1,1],1:end),zeros(2,size(aer_y,2)),'FaceAlpha',0.0,'EdgeAlpha',1.0,'FaceColor','white') 

for ii = 1:16
    surf(aer_x([ii,ii],:),[aer_y(ii,1:end);aer_y(ii,1:end)-dP],aer_z([ii,ii],1:end),zeros(2,size(aer_y,2)),'FaceAlpha',0.0,'EdgeAlpha',1.0,'FaceColor','white') 
end

axis off

c = colorbar('southoutside');
c.Label.String = 'deformation, m';
c.Position = [0.295104166666667 0.347969264544457 0.233020833333333 0.025477497255763];


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
