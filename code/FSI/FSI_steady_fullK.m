function wing_def = FSI_steady_fullK(paramFSI, simSteadyParam, simParam, wingModelStructure, wingDesign, wingModelAero)

FSI_curIT = 0;

forceACTinp = simSteadyParam.forceACTinp;

analysisParams = paramFSI.analysisParams;

b = paramFSI.wingParams.b;
meanChord = paramFSI.wingParams.chord;

y_segLimit = paramFSI.y_segLimit;
y_segMidpoint = paramFSI.y_segMidpoint;

FE_Grid = wingModelStructure.FE_Grid;

xMesh_loop = paramFSI.xMesh_loop;
yMesh_loop = paramFSI.yMesh_loop;
zMesh_loop = paramFSI.zMesh_loop;

cLdistr_yPos = paramFSI.cLdistr_yPos;

pInterU = wingModelAero.wingModelAeroFullK.pInterU;
pInterD = wingModelAero.wingModelAeroFullK.pInterD;

GtpsSU = wingModelAero.wingModelAeroFullK.GtpsSU;
GtpsSD = wingModelAero.wingModelAeroFullK.GtpsSD;

IDWu = wingModelAero.wingModelAeroFullK.IDWu;
IDWd = wingModelAero.wingModelAeroFullK.IDWd;

V = simSteadyParam.V; % m/s
alpha = simSteadyParam.alpha*pi/180;
vinf = [V*cos(alpha); 0; V*sin(alpha)];
restart = 0;

K_ASET = wingModelStructure.K_ASET;
nDOF = length(wingModelStructure.DOF_ASET) + length(wingModelStructure.DOF_RSET);
nDOF_ASET = length(wingModelStructure.DOF_ASET);
forces_ASET_L = zeros(nDOF_ASET, 1);

analysisParams.fiT = zeros(size(paramFSI.K_MODES,1),1);

while true
	FSI_curIT = FSI_curIT + 1;
	
    %% apply actuation
    if FSI_curIT == 1
        forces_ACT = zeros(nDOF/6, 6);
        % left wing
        ActFront_Pos1 = zeros(2*wingDesign.nRibs,3);
        ActRear_Pos1 = zeros(2*wingDesign.nRibs,3);
        nV = zeros(2*wingDesign.nRibs,3);

        force_ACT_t_L = -simParam.actForce*forceACTinp.L + simParam.actForce*forceACTinp.fs;
        force_ACT_t_R = simParam.actForce*forceACTinp.L + simParam.actForce*forceACTinp.fs;

        for iAct = 2*wingDesign.nRibs:-2:2*wingDesign.nRibs-2*wingDesign.nRibsC+2
            ActFront_ID(iAct) = find(FE_Grid.X(:) == wingModelStructure.actuation.actuationSparAP(iAct,1) & FE_Grid.Y(:) == wingModelStructure.actuation.actuationSparAP(iAct,2) & FE_Grid.Z(:) == wingModelStructure.actuation.actuationSparAP(iAct,3));
            ActRear_ID(iAct) = find(FE_Grid.X(:) == wingModelStructure.actuation.actuationsdvAP(iAct,1) & FE_Grid.Y(:) == wingModelStructure.actuation.actuationsdvAP(iAct,2) & FE_Grid.Z(:) == wingModelStructure.actuation.actuationsdvAP(iAct,3));
            ActFront_Pos1(iAct,:) = [FE_Grid.X(ActFront_ID(iAct)), FE_Grid.Y(ActFront_ID(iAct)), FE_Grid.Z(ActFront_ID(iAct))];
            ActRear_Pos1(iAct,:) = [FE_Grid.X(ActRear_ID(iAct)), FE_Grid.Y(ActRear_ID(iAct)), FE_Grid.Z(ActRear_ID(iAct))];
            nV(iAct,:) = (ActFront_Pos1(iAct,:)-ActRear_Pos1(iAct,:))./norm(ActFront_Pos1(iAct,:)-ActRear_Pos1(iAct,:));
            forces_ACT(ActFront_ID(iAct),1:3) = -nV(iAct,:)*force_ACT_t_L;
            forces_ACT(ActRear_ID(iAct),1:3) = nV(iAct,:)*force_ACT_t_L;
        end
        % right wing
        ActFront_Pos1R = zeros(2*wingDesign.nRibs,3);
        ActRear_Pos1R = zeros(2*wingDesign.nRibs,3);
        nVR = zeros(2*wingDesign.nRibs,3);
        for iAct = 2*wingDesign.nRibs:-2:2*wingDesign.nRibs-2*wingDesign.nRibsC+2
            ActFront_IDR(iAct) = find(FE_Grid.X(:) == wingModelStructure.actuation.actuationSparAP(iAct,1) & FE_Grid.Y(:) == wingModelStructure.actuation.actuationSparAP(iAct,2) & FE_Grid.Z(:) == -wingModelStructure.actuation.actuationSparAP(iAct,3));
            ActRear_IDR(iAct) = find(FE_Grid.X(:) == wingModelStructure.actuation.actuationsdvAP(iAct,1) & FE_Grid.Y(:) == wingModelStructure.actuation.actuationsdvAP(iAct,2) & FE_Grid.Z(:) == -wingModelStructure.actuation.actuationsdvAP(iAct,3));
            ActFront_Pos1R(iAct,:) = [FE_Grid.X(ActFront_IDR(iAct)), FE_Grid.Y(ActFront_IDR(iAct)), FE_Grid.Z(ActFront_IDR(iAct))];
            ActRear_Pos1R(iAct,:) = [FE_Grid.X(ActRear_IDR(iAct)), FE_Grid.Y(ActRear_IDR(iAct)), FE_Grid.Z(ActRear_IDR(iAct))];
            nVR(iAct,:) = (ActFront_Pos1R(iAct,:)-ActRear_Pos1R(iAct,:))./norm(ActFront_Pos1R(iAct,:)-ActRear_Pos1R(iAct,:));
            forces_ACT(ActFront_IDR(iAct),1:3) = -nVR(iAct,:)*force_ACT_t_R;
            forces_ACT(ActRear_IDR(iAct),1:3) = nVR(iAct,:)*force_ACT_t_R;
        end
%         ACTparam.ActFront_ID = ActFront_ID;
%         ACTparam.ActRear_ID = ActRear_ID;
%         ACTparam.ActFront_IDR = ActFront_IDR;
%         ACTparam.ActRear_IDR = ActRear_IDR;
% 
%         ACTparam.nV = nV;
%         ACTparam.nVR = nVR;
% 
%         wing_def.ACTparam = ACTparam;

        forces_ACT_transp = forces_ACT';
        forces_ACT_L = forces_ACT_transp(:);
        forces_ACT_L(wingModelStructure.DOF_RSET) = [];
    end

    %% calculate displacement
    forces_ASET_L = forces_ASET_L + forces_ACT_L;

    x_ASET = K_ASET\forces_ASET_L;
    newCoordinatesL = x_ASET;

    x_GLOB = zeros(nDOF, 1);
    x_GLOB(wingModelStructure.DOF_ASET) = x_ASET;
    x_GLOB = reshape(x_GLOB', 6, nDOF/6)';

    
    if restart && NASPAN_curIT == 1
%         iteration(NASPAN_curIT).iterations_error = norm(newCoordinatesLq - NM.xi)/sqrt(length(newCoordinatesLq));
        newCoordinatesLqConv = newCoordinatesL(abs(newCoordinatesL)>10^-5); 
        oldCoordinatesLqConv = NM.xi(abs(newCoordinatesL)>10^-5); 
        iteration(NASPAN_curIT).iterations_error = max(abs((newCoordinatesLqConv - oldCoordinatesLqConv)./newCoordinatesLqConv));
    end
    
    if restart==0 && FSI_curIT == 1
        iteration(FSI_curIT).iterations_error = inf;                    % We just started the simulation
    end
    
    if FSI_curIT > 1
%         iteration(NASPAN_curIT).iterations_error = norm(newCoordinatesLq - oldCoordinatesLq)/sqrt(length(newCoordinatesLq));
        % don't consider values too close to zero for convergence
        newCoordinatesLqConv = newCoordinatesL(abs(newCoordinatesL)>10^-5); 
        oldCoordinatesLqConv = oldCoordinatesL(abs(newCoordinatesL)>10^-5); 
        iteration(FSI_curIT).iterations_error = max(abs((newCoordinatesLqConv - oldCoordinatesLqConv)./newCoordinatesLqConv));
%         iteration(NASPAN_curIT).iterations_errorRel
    end
    oldCoordinatesL = newCoordinatesL;
    

    %% Run Panel Method
    
    % TPS INTERPOLATION
    skinNodesXU_def_tps = x_GLOB(pInterU,1); %up
    skinNodesYU_def_tps = x_GLOB(pInterU,2);
    skinNodesZU_def_tps = x_GLOB(pInterU,3);
    skinNodesXD_def_tps = x_GLOB(pInterD,1); %lower
    skinNodesYD_def_tps = x_GLOB(pInterD,2);
    skinNodesZD_def_tps = x_GLOB(pInterD,3);
    
    xTPSU = GtpsSU*skinNodesXU_def_tps;
    yTPSU = GtpsSU*skinNodesYU_def_tps;
    zTPSU = GtpsSU*skinNodesZU_def_tps;
    xTPSD = GtpsSD*skinNodesXD_def_tps;
    yTPSD = GtpsSD*skinNodesYD_def_tps;
    zTPSD = GtpsSD*skinNodesZD_def_tps;

    size1U = analysisParams.n_seg_PM*2 + 1;  
    size2U = analysisParams.num_airfoil_nodes_panel;
    
    xTPSsU = reshape(xTPSU, size1U, size2U);
    yTPSsU = reshape(yTPSU, size1U, size2U);
    zTPSsU = reshape(zTPSU, size1U, size2U);
    xTPSsD = reshape(xTPSD, size1U, size2U);
    yTPSsD = reshape(yTPSD, size1U, size2U);
    zTPSsD = reshape(zTPSD, size1U, size2U);

%     yTPSU_tot = yMesh_loop(:,1:analysisParams.num_airfoil_nodes_panel) + rot90(yTPSsU,2); 
%     xTPSU_tot = xMesh_loop(:,1:analysisParams.num_airfoil_nodes_panel) + rot90(xTPSsU,2);
%     zTPSU_tot = zMesh_loop(:,1:analysisParams.num_airfoil_nodes_panel) + rot90(zTPSsU,2);
%     
%     yTPSD_tot = yMesh_loop(:,analysisParams.num_airfoil_nodes_panel+1:end) + rot90(yTPSsD(:,2:end-1),2);  
%     xTPSD_tot = xMesh_loop(:,analysisParams.num_airfoil_nodes_panel+1:end) + rot90(xTPSsD(:,2:end-1),2);
%     zTPSD_tot = zMesh_loop(:,analysisParams.num_airfoil_nodes_panel+1:end) + rot90(zTPSsD(:,2:end-1),2);
    yTPSU_tot = yMesh_loop(:,1:analysisParams.num_airfoil_nodes_panel) + flip(flip(yTPSsU,1),2);
    xTPSU_tot = xMesh_loop(:,1:analysisParams.num_airfoil_nodes_panel) + flip(flip(xTPSsU,1),2);
    zTPSU_tot = zMesh_loop(:,1:analysisParams.num_airfoil_nodes_panel) + flip(flip(zTPSsU,1),2);
    
    yTPSD_tot = yMesh_loop(:,analysisParams.num_airfoil_nodes_panel+1:end) + flip(flip(yTPSsD(:,2:end-1),1),2);
    xTPSD_tot = xMesh_loop(:,analysisParams.num_airfoil_nodes_panel+1:end) + flip(flip(xTPSsD(:,2:end-1),1),2);
    zTPSD_tot = zMesh_loop(:,analysisParams.num_airfoil_nodes_panel+1:end) + flip(flip(zTPSsD(:,2:end-1),1),2);
     
    xTPSs = [xTPSU_tot, xTPSD_tot];
    yTPSs = [yTPSU_tot, yTPSD_tot];
    zTPSs = [zTPSU_tot, zTPSD_tot];    

    % Change wing reference system
    aer_x = xTPSs;
    aer_y = -zTPSs;
    aer_z = yTPSs;
    
    [M1,N1] = size(aer_x);
    x = [];
    y = [];
    z = [];
    x(1:N1,1:M1) = aer_x';
    y(1:N1,1:M1) = aer_y';
    z(1:N1,1:M1) = aer_z';

    % close trailing edge
    N1 = N1+1;
    x = [x;x(1,:)];
    y = [y;y(1,:)];
    z = [z;z(1,:)];

    % Flip for the correct orientation
    x = flip(x);
    y = flip(y);
    z = flip(z);
    
    M = M1-1;
    N = N1-1;     
    N2 = N+2;
    farpoint = 1000;
    error = 1e-6;               % minimal distance (error)
    farfieldfaktor = 5;             % "far field" distance in longer panel diagonal

    for j = 1:M1
        x(N2,j) = x(N1,j)+farpoint;
        y(N2,j) = y(N1,j);
        z(N2,j) = z(N1,j);
    end
       
    Ax = (x(2:N1+1,2:M+1)-x(1:N1,1:M));
    Ay = (y(2:N1+1,2:M+1)-y(1:N1,1:M)); 
    Az = (z(2:N1+1,2:M+1)-z(1:N1,1:M));
    Bx = (x(2:N1+1,1:M)-x(1:N1,2:M+1));
    By = (y(2:N1+1,1:M)-y(1:N1,2:M+1));
    Bz = (z(2:N1+1,1:M)-z(1:N1,2:M+1));

    FF = farfieldfaktor*max(realsqrt(Ax.^2+Ay.^2+Az.^2),realsqrt(Bx.^2+By.^2+Bz.^2));
    
    C1 = Ay.*Bz-Az.*By;
    C2 = Az.*Bx-Ax.*Bz;
    C3 = Ax.*By-Ay.*Bx;
    modc = sqrt(C1.^2 + C2.^2 + C3.^2);
    S = modc/2;
    modc(modc==0) = 10^-20;
%     n1 = -C1./modc;
%     n2 = -C2./modc;
%     n3 = -C3./modc;
    n1 = C1./modc;
    n2 = C2./modc;
    n3 = C3./modc;
    % calculation of colocation points
    cx = (x(1:N1,1:M) + x(1:N1,2:M+1) + x(2:N1+1,1:M) + x(2:N1+1,2:M+1))/4;
    cy = (y(1:N1,1:M) + y(1:N1,2:M+1) + y(2:N1+1,1:M) + y(2:N1+1,2:M+1))/4;
    cz = (z(1:N1,1:M) + z(1:N1,2:M+1) + z(2:N1+1,1:M) + z(2:N1+1,2:M+1))/4;
    % calculation of u(longitudinal), p(transversal) and o(perpendicular) vectors
    ux = ((x(2:N1+1,1:M)+x(2:N1+1,2:M+1))-(x(1:N1,1:M)+x(1:N1,2:M+1)))/2;
    uy = ((y(2:N1+1,1:M)+y(2:N1+1,2:M+1))-(y(1:N1,1:M)+y(1:N1,2:M+1)))/2;
    uz = ((z(2:N1+1,1:M)+z(2:N1+1,2:M+1))-(z(1:N1,1:M)+z(1:N1,2:M+1)))/2;
    uu = sqrt(ux.^2 + uy.^2 + uz.^2);
    u1 = ux./uu;
    u2 = uy./uu;
    u3 = uz./uu;
%     px = ((x(1:N1,2:M+1)+x(2:N1+1,2:M+1))-(x(1:N1,1:M)+x(2:N1+1,1:M)))/2;
%     py = ((y(1:N1,2:M+1)+y(2:N1+1,2:M+1))-(y(1:N1,1:M)+y(2:N1+1,1:M)))/2;
%     pz = ((z(1:N1,2:M+1)+z(2:N1+1,2:M+1))-(z(1:N1,1:M)+z(2:N1+1,1:M)))/2;
    px = -((x(1:N1,2:M+1)+x(2:N1+1,2:M+1))-(x(1:N1,1:M)+x(2:N1+1,1:M)))/2;
    py = -((y(1:N1,2:M+1)+y(2:N1+1,2:M+1))-(y(1:N1,1:M)+y(2:N1+1,1:M)))/2;
    pz = -((z(1:N1,2:M+1)+z(2:N1+1,2:M+1))-(z(1:N1,1:M)+z(2:N1+1,1:M)))/2;
    pp = sqrt(px.^2 + py.^2 + pz.^2);
    p1 = px./pp;
    p2 = py./pp;
    p3 = pz./pp;
    o1 = n2(1:N1,1:M).*u3(1:N1,1:M) - n3(1:N1,1:M).*u2(1:N1,1:M);
    o2 = u1(1:N1,1:M).*n3(1:N1,1:M) - n1(1:N1,1:M).*u3(1:N1,1:M);
    o3 = n1(1:N1,1:M).*u2(1:N1,1:M) - n2(1:N1,1:M).*u1(1:N1,1:M);

    % Clockwise in order to have a positive doublet in the normal direction
    x1 = (x(1:end-1,2:end)-cx).*u1+(y(1:end-1,2:end)-cy).*u2+(z(1:end-1,2:end)-cz).*u3;
    y1 = (x(1:end-1,2:end)-cx).*o1+(y(1:end-1,2:end)-cy).*o2+(z(1:end-1,2:end)-cz).*o3;
    x2 = (x(1:end-1,1:end-1)-cx).*u1+(y(1:end-1,1:end-1)-cy).*u2+(z(1:end-1,1:end-1)-cz).*u3;
    y2 = (x(1:end-1,1:end-1)-cx).*o1+(y(1:end-1,1:end-1)-cy).*o2+(z(1:end-1,1:end-1)-cz).*o3;
    x3 = (x(2:end,1:end-1)-cx).*u1+(y(2:end,1:end-1)-cy).*u2+(z(2:end,1:end-1)-cz).*u3;
    y3 = (x(2:end,1:end-1)-cx).*o1+(y(2:end,1:end-1)-cy).*o2+(z(2:end,1:end-1)-cz).*o3;
    x4 = (x(2:end,2:end)-cx).*u1+(y(2:end,2:end)-cy).*u2+(z(2:end,2:end)-cz).*u3;
    y4 = (x(2:end,2:end)-cx).*o1+(y(2:end,2:end)-cy).*o2+(z(2:end,2:end)-cz).*o3;
    
    d1 = realsqrt((x2-x1).^2 + (y2-y1).^2);
    d2 = realsqrt((x3-x2).^2 + (y3-y2).^2);
    d3 = realsqrt((x4-x3).^2 + (y4-y3).^2);
    d4 = realsqrt((x1-x4).^2 + (y1-y4).^2);
    
    % setting source strength
    sigma = zeros(M*N,1);
    K = 0;
    for i = 1:N
        for j = 1:M
            K = K+1;
%             sigma(K) = n1(i,j)*vinf(1) + n2(i,j)*vinf(2) + n3(i,j)*vinf(3);
            sigma(K) = -(n1(i,j)*vinf(1) + n2(i,j)*vinf(2) + n3(i,j)*vinf(3));
        end
    end

    aeroMesh.aer_x = aer_x;
    aeroMesh.aer_y = aer_y;
    aeroMesh.aer_z = aer_z;
%     [aPM,bPM] = apame_panel_method(aeroMesh);  % the function computes the AIC
%     error = 0.000001;               % minimal distance (error)
    [aPM,bPM]=influence_self_lift(FF(1:N1,1:M),error,S(1:N1,1:M),M,N,N1,M*N,u1(1:N1,1:M),u2(1:N1,1:M),u3(1:N1,1:M),o1(1:N1,1:M),o2(1:N1,1:M),o3(1:N1,1:M),n1(1:N1,1:M),n2(1:N1,1:M),n3(1:N1,1:M),cx(1:N1,1:M),cy(1:N1,1:M),cz(1:N1,1:M),x1(1:N1,1:M),y1(1:N1,1:M),x2(1:N1,1:M),y2(1:N1,1:M),x3(1:N1,1:M),y3(1:N1,1:M),x4(1:N1,1:M),y4(1:N1,1:M),d1(1:N1,1:M),d2(1:N1,1:M),d3(1:N1,1:M),d4(1:N1,1:M));

    %%%%%%%%%%%%% SOLVING SYSTEM OF EQUATIONS %%%%%%%%%%%%%
%     RHS = bPM*sigma;
    RHS = -bPM*sigma;
    gamma = aPM\RHS;
    wing_def.gammaAPAME = gamma;
    
    % returning dipole values to i,j matrix
   	gamF = zeros(N,M);
    K = 0;
    for k = 1:N
        for l = 1:M
            K = K+1;
            gamF(k,l) = gamma(K);
        end
    end
    
    % calculation of speeds, pressures and forces
    q = 0.5*1.225*V^2;
    % calculation of induced speeds N-direction
    dx = zeros(N,M);
    dy = zeros(N,M);
    dz = zeros(N,M);
    quF = zeros(N,M);
    % boundary panels - first order interpolation
    dx(1,:) = cx(2,:)-cx(1,:);
    dy(1,:) = cy(2,:)-cy(1,:);
    dz(1,:) = cz(2,:)-cz(1,:);
    quF(1,:) = gamF(2,:)-gamF(1,:);
    dx(N,:) = cx(N,:)-cx(N-1,:);
    dy(N,:) = cy(N,:)-cy(N-1,:);
    dz(N,:) = cz(N,:)-cz(N-1,:);
    quF(N,:) = gamF(N,:)-gamF(N-1,:);
    % longitudinal
    dx(2:N-1,:) = cx(3:N,:)-cx(1:N-2,:);
    dy(2:N-1,:) = cy(3:N,:)-cy(1:N-2,:);
    dz(2:N-1,:) = cz(3:N,:)-cz(1:N-2,:);
    quF(2:N-1,:) = gamF(3:N,:)-gamF(1:N-2,:);
 
    dru = sqrt(dx.^2+dy.^2+dz.^2);
    quF = quF./dru;

    % calculation of induced speeds M-direction
    dx = zeros(N,M);
    dy = zeros(N,M);
    dz = zeros(N,M);
    qpF = zeros(N,M);
    % boundary panels - first order interpolation
    dx(:,1) = cx(1:N,2)-cx(1:N,1);
    dy(:,1) = cy(1:N,2)-cy(1:N,1);
    dz(:,1) = cz(1:N,2)-cz(1:N,1);
    qpF(:,1) = gamF(:,2)-gamF(:,1);
    dx(:,M) = cx(1:N,M)-cx(1:N,M-1);
    dy(:,M) = cy(1:N,M)-cy(1:N,M-1);
    dz(:,M) = cz(1:N,M)-cz(1:N,M-1);
    qpF(:,M) = gamF(:,M)-gamF(:,M-1);
    % transversal
    dx(:,2:M-1) = cx(1:N,3:M)-cx(1:N,1:M-2);
    dy(:,2:M-1) = cy(1:N,3:M)-cy(1:N,1:M-2);
    dz(:,2:M-1) = cz(1:N,3:M)-cz(1:N,1:M-2);
    qpF(:,2:M-1) = gamF(:,3:M)-gamF(:,1:M-2);
    
    drp = sqrt(dx.^2+dy.^2+dz.^2);
%     qpF = qpF./drp;
    qpF = -qpF./drp;
    
    qoF = (p1(1:N,:).*o1(1:N,:) + p2(1:N,:).*o2(1:N,:) + p3(1:N,:).*o3(1:N,:)).*qpF;

    gu = u1(1:N,:)*vinf(1) + u2(1:N,:)*vinf(2) + u3(1:N,:)*vinf(3);
    go = o1(1:N,:)*vinf(1) + o2(1:N,:)*vinf(2) + o3(1:N,:)*vinf(3);

    vxF = (-quF + gu).*u1(1:N,:) + (-qoF + go).*o1(1:N,:);
    vyF = (-quF + gu).*u2(1:N,:) + (-qoF + go).*o2(1:N,:);
    vzF = (-quF + gu).*u3(1:N,:) + (-qoF + go).*o3(1:N,:);
    vF = sqrt(vxF.^2+vyF.^2+vzF.^2);
 
    cp = 1-vF.^2/V^2;
    dX = -cp.*S(1:N,:).*n1(1:N,:)*q;
    dY = -cp.*S(1:N,:).*n2(1:N,:)*q;  %in the wing frame of ref
    dZ = -cp.*S(1:N,:).*n3(1:N,:)*q;
    
    FX = sum(sum(dX));
    FY = sum(sum(dY));
    FZ = sum(sum(dZ));
    FL = sum(sum(-dY.*cz(1:N,:)+dZ.*cy(1:N,:)));  % roll
    FM = sum(sum( dX.*cz(1:N,:)-dZ.*(cx(1:N,:)+paramFSI.cog(1)))); % pitch
    FN = sum(sum(-dX.*cy(1:N,:)+dY.*cx(1:N,:)));  % yaw

    dXU = dX(1:size(dX,1)/2+1,:);
    dYU = dY(1:size(dX,1)/2+1,:);
    dZU = dZ(1:size(dX,1)/2+1,:);
    dXD = dX(size(dX,1)/2:end,:);
    dYD = dY(size(dX,1)/2:end,:);
    dZD = dZ(size(dX,1)/2:end,:);
    
    % inverse distance weighting
    FAeroUIDW = [dXU(:),dZU(:),dYU(:)];
    FAeroDIDW = [dXD(:),dZD(:),dYD(:)];

    FStructUIDW = IDWu*FAeroUIDW;
    FStructDIDW = IDWd*FAeroDIDW;
    
    forcesIDWi = zeros(6,nDOF/6);
    forcesIDWi(1:3,pInterU) = FStructUIDW';
    forcesIDWi(1:3,pInterD) = FStructDIDW';
    forcesIDW = forcesIDWi(:);
    
    forcesIDW_tot = forcesIDW;
    forcesIDW_tot(wingModelStructure.DOF_RSET) = [];  % R_set is the constrained nodes
    
    forces_ASET_L = forcesIDW_tot;
    
    Sw = meanChord*b;
    FZA = FZ*cos(alpha)+FX*sin(alpha);
%     FZA = FZ*cos(alpha)-FX*sin(alpha);

    wing_def.alpha = alpha(end);
    wing_def.cL = FZA/q/Sw;
    wing_def.cD0 = 0;
    wing_def.cRoll = FL/q/Sw/b;
    wing_def.cPitch = FM/q/Sw/meanChord;
    wing_def.cYaw = FN/q/Sw/b;

	%% Convergency check on the pressure distribution   
    cP_conv_new_up = dZU;
    cP_conv_new_dn = dZD;
	
    if restart==0 && FSI_curIT == 1
%         iteration(FSI_curIT).iterations_error_cP = inf;                 % We just started the simulation
        iteration(FSI_curIT).iterations_error_cPrel = inf;                 % We just started the simulation
    else
%         iteration(FSI_curIT).iterations_error_cP = norm((cP_conv_new_up(2:end-1,:) - flipud(cP_conv_new_dn(2:end-1,:)) - (cP_conv_old_up(2:end-1,:) - flipud(cP_conv_old_dn(2:end-1,:)))),2); % do not compare cP on leading and trailing edge...
        convUp = max(max(abs((cP_conv_new_up(2:end-1,:) - cP_conv_old_up(2:end-1,:))./cP_conv_new_up(2:end-1,:)))); % do not compare cP on leading and trailing edge...
        convDn = max(max(abs((cP_conv_new_dn(2:end-1,:) - cP_conv_old_dn(2:end-1,:))./cP_conv_new_dn(2:end-1,:)))); % do not compare cP on leading and trailing edge...
        iteration(FSI_curIT).iterations_error_cPrel = max([convUp,convDn]); % do not compare cP on leading and trailing edge...
    end
    
    if convergencyCheck(iteration(FSI_curIT), analysisParams)
        wing_def.oldCoordinatesL = newCoordinatesL;
        wing_def.cP_conv_old_up = cP_conv_new_up;
        wing_def.cP_conv_old_dn = cP_conv_new_dn;
        break;
    end
    
    if FSI_curIT > 100
        wing_def.oldCoordinatesL = newCoordinatesL;
        wing_def.cP_conv_old_up = cP_conv_new_up;
        wing_def.cP_conv_old_dn = cP_conv_new_dn;
        warning('!!!!!! FSI not converged !!!!!!!!!!')
        break
    end
    
	cP_conv_old_up = cP_conv_new_up;
	cP_conv_old_dn = cP_conv_new_dn;
    
    % convergence in first loop for aerodynamics only simulation
    if paramFSI.aeroOnly
        wing_def.oldCoordinatesL = newCoordinatesL;
        wing_def.cP_conv_old_up = cP_conv_new_up;
        wing_def.cP_conv_old_dn = cP_conv_new_dn;
        break
    end

end

%% recalculate induced drag with ellt
CX = sum(dX,1)./sum(S(1:end-1,:),1)/q*2;
CZ = sum(dZ,1)./sum(S(1:end-1,:),1)/q*2;

cLdistr_val = CZ*cos(alpha)-CX*sin(alpha);
cLdistr_val = cLdistr_val';

% if simSteadyParam.doPlot
%     figure('name', sprintf('Spanwise lift distribution from panel method analysis, it. %d', FSI_curIT))
%     plot(cLdistr_yPos, cLdistr_val);
%     xlabel('y')
%     ylabel('c_L')
%     title('Spanwise lift distribution from panel method analysis')
% end

cLdistr_val_interp = interp1(cLdistr_yPos, cLdistr_val, y_segMidpoint, 'linear', 'extrap');		
K_alpha_ind_ELLT = paramFSI.K_alpha_ind_ELLT;
[eLLT_cDi,wing_def.cYaw,wing_def.FN] = eLLT_PM_cDi(y_segLimit, meanChord, V, cLdistr_val_interp, K_alpha_ind_ELLT,q,Sw,b);

% Overwrite panel method induced drag with the one calculated by means of the lifting line
wing_def.cDi = eLLT_cDi;

end

