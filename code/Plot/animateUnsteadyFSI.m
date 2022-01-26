function animateUnsteadyFSI(simInput,simUnsteadyParam)

% load wing and wake data
cOut = load(strcat(pwd,filesep,'data',filesep,simInput.paramFSI.wingParams.airfoil,filesep,'Animation',filesep,simInput.paramFSI.animationName,filesep,'Animation.csv'));
cOut = cOut';

M = simInput.paramFSI.M1; % span
N = simInput.paramFSI.N1+1; % chord
nWake = simInput.paramFSI.cutWakeAt; % wake
nT = size(cOut,2); % time

x = cOut(1:(M*(N+nWake)),:);
y = cOut(M*(N+nWake)+1:2*M*(N+nWake),:);
z = cOut(2*M*(N+nWake)+1:3*M*(N+nWake),:);
cpW = cOut(3*M*(N+nWake)+1:end,:); % pressure coefficient

for i = 1:nT
    xPlot(:,:,i) = reshape(x(:,i),N+nWake,M);
    yPlot(:,:,i) = reshape(y(:,i),N+nWake,M);
    zPlot(:,:,i) = reshape(z(:,i),N+nWake,M);
    cpWPlot(:,:,i) = reshape(cpW(:,i),N,M);
end

clear cOut x y z cpW


%% plot figure of wing and wake after initialization

% check if initialization is stored in Animation.csv
if nT == (simUnsteadyParam.time+simUnsteadyParam.timeInit)/simUnsteadyParam.dT
    ts = simUnsteadyParam.timeInit/simUnsteadyParam.dT;
else
    ts = 1;
end

geomPlotWingX = xPlot(1:N,:,ts);
geomPlotWingY = yPlot(1:N,:,ts);
geomPlotWingZ = zPlot(1:N,:,ts);
geomPlotWingZ0 = zPlot(1:N,:,1);

geomPlotWakeX = xPlot(N+1:end,:,ts);
geomPlotWakeY = yPlot(N+1:end,:,ts);
geomPlotWakeZ = zPlot(N+1:end,:,ts);

% limits for plots
maxWingDef = max(max(max(abs(zPlot(1:N,:,:)-zPlot(1:N,:,1)))));
limitColorbar = [0 maxWingDef];
xl = [-1 max(max(max(abs(xPlot))))+1];
yl = [-max(max(max(abs(yPlot))))-1 max(max(max(abs(yPlot))))+1];
zl = [-max(max(max(abs(zPlot))))-1 max(max(max(abs(zPlot))))+1];

figWakef = figure;
x0=50;
y0=10;
width=600;
height=300;
set(gcf,'units','points','position',[x0,y0,width,height])

Cwing = geomPlotWingZ-geomPlotWingZ0; % plot wing deformation
% Cwing = cpWPlot; % plot wing pressure coefficient

gW = 0.6;
nPw = nWake;
nP=2;
gS = 0.2;%1;
eA = 1;
lW = 0.5;%%1.2;
fA = 0;

mesh(geomPlotWingX,geomPlotWingY,geomPlotWingZ,Cwing,'FaceColor','interp','EdgeColor','interp','EdgeAlpha',gW,'FaceAlpha',gW); hold on
mesh(geomPlotWakeX(1:nP:nPw,:),geomPlotWakeY(1:nP:nPw,:),geomPlotWakeZ(1:nP:nPw,:),...
        'EdgeColor',gS*[1 1 1],'EdgeAlpha',eA,'FaceColor',gS*[1 1 1],'FaceAlpha',fA,'LineWidth',lW)
    
axis equal
grid off
axis off

% colormap(autumn)
% colormap(jet)
colormap(turbo)
caxis(limitColorbar)

xlim(xl)
ylim(yl)
zlim(zl)

view([-60 10]);
% view([-50 5]); 
set(gca, 'Projection','perspective')

% set(gcf,'color','k');
% set(gca,'color','none')
set(gca,'color','w')

% save figure
folderName = strcat(pwd,filesep,'data',filesep,simInput.paramFSI.wingParams.airfoil,filesep,'Animation',filesep,simInput.paramFSI.animationName,filesep);
figureName = [simInput.paramFSI.wingParams.airfoil '_animation'];
% saveas(figWakef,[folderName figureName],'epsc')
saveas(figWakef,[folderName figureName],'meta')
saveas(figWakef,[folderName figureName],'fig')



%% generate animation

baseFileName = [simInput.paramFSI.wingParams.airfoil '_animation'];
fullFileName = fullfile(folderName, baseFileName);

nS = 1;
numberOfFrames = nT;
figure;

% Set up the movie structure.
allTheFrames = cell(numberOfFrames,1);
vidHeight = 1500;
vidWidth = 3000;
allTheFrames(:) = {zeros(vidHeight, vidWidth, 3, 'uint8')};
allTheColorMaps = cell(numberOfFrames,1);
allTheColorMaps(:) = {zeros(256, 3)};
myMovie = struct('cdata', allTheFrames, 'colormap', allTheColorMaps);
set(gcf, 'renderer', 'zbuffer','units','points','position',[x0,y0,width,height])

% Create the movie.
% After this loop starts, BE SURE NOT TO RESIZE THE WINDOW AS IT'S SHOWING THE FRAMES, or else you won't be able to save it.
for frameIndex = 1:nS:nT
	
    geomPlotWingX = xPlot(1:N,:,frameIndex);
    geomPlotWingY = yPlot(1:N,:,frameIndex);
    geomPlotWingZ = zPlot(1:N,:,frameIndex);
    geomPlotWingZ0 = zPlot(1:N,:,1);

    geomPlotWakeX = xPlot(N+1:end,:,frameIndex);
    geomPlotWakeY = yPlot(N+1:end,:,frameIndex);
    geomPlotWakeZ = zPlot(N+1:end,:,frameIndex);
    Cwing = geomPlotWingZ-geomPlotWingZ0; % plot wing deformation
    
	cla reset; % reset plot

    mesh(geomPlotWingX,geomPlotWingY,geomPlotWingZ,Cwing,'FaceColor','interp','EdgeColor','interp','EdgeAlpha',gW,'FaceAlpha',gW); hold on
    mesh(geomPlotWakeX(1:nP:nPw,:),geomPlotWakeY(1:nP:nPw,:),geomPlotWakeZ(1:nP:nPw,:),...
            'EdgeColor',gS*[1 1 1],'EdgeAlpha',eA,'FaceColor',gS*[1 1 1],'FaceAlpha',fA,'LineWidth',lW)

    axis equal
    grid off
    axis off

    colormap(turbo)

    caxis(limitColorbar)

    % view([-50 5]); 
    view([-60 10]); 
    
    xlim(xl)
    ylim(yl)
    zlim(zl)
    
    set(gca, 'Projection','perspective')

    % set(gcf,'color','k');
    % set(gca,'color','none')
    set(gca,'color','w')
    
	drawnow;
	thisFrame = getframe(gca);
	myMovie(frameIndex) = thisFrame;
    
    frameGif = getframe(gcf);
    im = frame2im(frameGif);
    [imind,cm] = rgb2ind(im,256);
    delayTgif = 0.05;
    % Write to the GIF File
    if frameIndex == 1
        imwrite(imind,cm,[fullFileName '.gif'],'gif', 'Loopcount',inf,'DelayTime',delayTgif);
    else
        imwrite(imind,cm,[fullFileName '.gif'],'gif','WriteMode','append','DelayTime',delayTgif);
    end
end


% save video
profile = 'MPEG-4';
% profile = 'Uncompressed AVI';

writerObj = VideoWriter(fullFileName, profile);
writerObj.FrameRate = round(nT/(simUnsteadyParam.time+simUnsteadyParam.timeInit)/nS); % video length is equal to simulation time
open(writerObj);
for frameNumber = 1:nS:nT
   writeVideo(writerObj, myMovie(frameNumber));
end
close(writerObj);


