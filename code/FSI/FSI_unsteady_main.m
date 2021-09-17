function FSI_OUT = FSI_unsteady_main(input,paramFSI,inputCreate,deltaAct,inputMorphing,createROM,inputSnapshot0)

% store snapshots
if (~exist('createROM', 'var'))
    createROM = false;
end

% check if this FSI run generates data for the calculation of the empirical
% observaibility gramian
if ~isfield(paramFSI,'obsGramian')
    paramFSI.obsGramian.run = false;
end

% check if aero only or if FSI is run
if ~isfield(paramFSI,'aeroOnly')
    paramFSI.aeroOnly = false;
end

analysisParams = paramFSI.analysisParams;

% stiffness, mass, and damping matrix
K_MODES = paramFSI.K_MODES;
M_MODES = paramFSI.M_MODES;
C_MODES = paramFSI.C_MODES;

analysisParams.V = norm(input(1:3)); % m/s
analysisParams.alpha = atan2(input(3),input(1))+(inputCreate.alphaWing)*pi/180;
inputSnapshot0(6) = inputSnapshot0(6)+(inputCreate.alphaWing)*pi/180;
analysisParams.cutWakeAt = paramFSI.cutWakeAt;

analysisParams.Omega = input(4:6);
% % unsteady aero: delta alpha induced pitch rate
% if isfield(inputCreate.restartIN,'alpha')
%     analysisParams.Omega(2) = analysisParams.Omega(2) + (analysisParams.alpha - inputCreate.restartIN.alpha)/deltaAct.DsimTime;
% end
    
analysisParams.vinf = input(1:3);
analysisParams.deltaT = deltaAct.DsimTime;

forceACTinp.L = input(7);
forceACTinp.fs = inputMorphing;

analysisParams.cM_pos = paramFSI.cog; % center of gravity, used to calculate pitch moment

%% Take care of restarting or not from previous time step

if ~isempty(inputCreate.restartIN)
    analysisParams.oldCoordinatesL = inputCreate.restartIN.oldCoordinatesL;
    analysisParams.cP_conv_old_up = inputCreate.restartIN.cP_conv_old_up;
    analysisParams.cP_conv_old_dn = inputCreate.restartIN.cP_conv_old_dn;
    analysisParams.fiT = inputCreate.restartIN.fiT;
    analysisParams.fiACT = inputCreate.restartIN.fiACT;
    analysisParams.restart = true;
else
    analysisParams.restart = false;
    analysisParams.fiT = zeros(size(K_MODES,1),1);
    analysisParams.fiACT = zeros(size(K_MODES,1),1);
end

% initialize wake in case it is the first iteration
if paramFSI.firstIt == 1
    inputCreate.restartIN.Wake_Information.X = [];
    inputCreate.restartIN.Wake_Information.Y = [];
    inputCreate.restartIN.Wake_Information.Z = [];
    inputCreate.restartIN.Wake_Information.Strength = [];
    inputCreate.restartIN.Wake_Information.Gamma_Old = zeros(analysisParams.num_airfoil_nodes_panel*2-2,analysisParams.n_seg_PM*2);
    for i = 1:paramFSI.cutWakeAt
        for j = 1:analysisParams.n_seg_PM*2+1
            inputCreate.restartIN.Wake_Information.X(i,j) = paramFSI.wingParams.chord + analysisParams.deltaT*analysisParams.V*(i+1);
        end
    end
    for i = 1:paramFSI.cutWakeAt
        inputCreate.restartIN.Wake_Information.Y(i,:) = linspace(paramFSI.wingParams.b/2,-paramFSI.wingParams.b/2,(analysisParams.n_seg_PM*2+1));
    end
    inputCreate.restartIN.Wake_Information.Z = zeros(paramFSI.cutWakeAt,analysisParams.n_seg_PM*2+1);
    inputCreate.restartIN.Wake_Information.Strength = zeros((paramFSI.cutWakeAt)*(analysisParams.n_seg_PM*2),1);
    
end



%% Perform the FSI
FSI_out  = FSI_unsteady(analysisParams, K_MODES, forceACTinp, paramFSI, M_MODES, C_MODES,inputCreate.restartIN.Wake_Information,createROM,inputSnapshot0);


%% Take care of the outputs
restartOUT.oldCoordinatesL = FSI_out.oldCoordinatesL;
restartOUT.cP_conv_old_up = FSI_out.cP_conv_old_up;
restartOUT.cP_conv_old_dn = FSI_out.cP_conv_old_dn;
restartOUT.FSI_out = FSI_out;

restartOUT.fiT = FSI_out.fiT;
restartOUT.fiACT = FSI_out.fiACT;
restartOUT.Wake_Information = FSI_out.Wake_Information;
restartOUT.alpha = analysisParams.alpha;

FSI_OUT = FSI_out;
FSI_OUT.restartIN = restartOUT;

end