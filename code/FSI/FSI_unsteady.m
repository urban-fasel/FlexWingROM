function wing_def = FSI_unsteady(analysisParams, K_MODES, forceACTinp, paramFSI, M_MODES, C_MODES, Wake_Information , createROM, inputSnapshot0)

cutWakeAt = analysisParams.cutWakeAt;

FSI_curIT = 0; % current iteration of FSI loop
b = paramFSI.wingParams.b;
y_segLimit = paramFSI.y_segLimit;
y_segMidpoint = paramFSI.y_segMidpoint;

meanChord = paramFSI.wingParams.chord;
xMesh_loop = paramFSI.xMesh_loop;
yMesh_loop = paramFSI.yMesh_loop;
zMesh_loop = paramFSI.zMesh_loop;
cLdistr_yPos = paramFSI.cLdistr_yPos;

V = analysisParams.V; % m/s
alpha = analysisParams.alpha;
vinf = [V*cos(alpha); 0; V*sin(alpha)];
mu = 1.81*10^-5; % viscosity air kg/(m*s)

restart = analysisParams.restart;

if restart
    cP_conv_old_up = analysisParams.cP_conv_old_up;
    cP_conv_old_dn = analysisParams.cP_conv_old_dn;
end

NM = paramFSI.NM;
b1 = NM.b1;
b2 = NM.b2;
b3 = NM.b3;
b4 = NM.b4;
b5 = NM.b5;

Wake_Information.X(end,:) = [];
Wake_Information.Y(end,:) = [];
Wake_Information.Z(end,:) = [];
Wake_Information.Strength(end-(paramFSI.M1-1)+1:end)=[];

error = 1e-6;           % minimal distance (error)
farfieldfaktor = 5;     % "far field" distance in longer panel diagonal


while true
    FSI_curIT = FSI_curIT + 1;
    
    %% Calculate displacements as sum of actuation + deformation
    
    % actuation force in modal coordinates, roll control
    if forceACTinp.L >= 0
        fiACT = paramFSI.fiACTMax*forceACTinp.L;
    else
        fiACT = -paramFSI.fiACTMin*forceACTinp.L;
    end
    
    % actuation force in modal coordinates, lift control
    if forceACTinp.fs >= 0
        fiACT_LA = paramFSI.fiACT_LAMax*forceACTinp.fs;
    else
        fiACT_LA = -paramFSI.fiACT_LAMin*forceACTinp.fs;
    end
    
    % sum of roll and lift control actuation
    fiACT = fiACT + fiACT_LA;
    
    % sum of actuation and aerodynamic forces ikn modal coordinates
    if FSI_curIT == 1
        fi1 = analysisParams.fiT + fiACT;
    else
        fi1 = fiT + fiACT;
    end
    
    % Newmark beta method: book rixen p.523, book hughes p.491
    dtilde = (NM.xi + b1*NM.xdi + b2*NM.xddi);
    vtilde = (NM.xdi + b3*NM.xddi);
    xddi1 = (M_MODES + b4*C_MODES + b5*K_MODES)\(fi1 -  C_MODES*vtilde - K_MODES*dtilde); % modal accelerations
    xi1 = dtilde + b5*xddi1; % modal amplitudes
    xdi1 = vtilde + b4*xddi1; % modal velocities
    
    % initial conditions in case this FSI run generates data for the calculation of the empirical observaibility gramian
    if paramFSI.obsGramian.run
        nStatesAero = 608;
        if paramFSI.obsGramian.nState > nStatesAero && FSI_curIT == 1
            xi1(paramFSI.obsGramian.nState-nStatesAero) = xi1(paramFSI.obsGramian.nState-nStatesAero)*paramFSI.obsGramian.InpAmp(1); % 1.2
            xi1_0 = xi1(paramFSI.obsGramian.nState-nStatesAero);
        elseif paramFSI.obsGramian.nState > nStatesAero && FSI_curIT > 1
            xi1(paramFSI.obsGramian.nState-nStatesAero) = xi1_0;
        end
    end
    
    % modal amplitudes saved for convergence check
    newCoordinatesLq = xi1;
    if restart && FSI_curIT == 1
        newCoordinatesLqConv = newCoordinatesLq(abs(newCoordinatesLq)>10^-5); 
        oldCoordinatesLqConv = NM.xi(abs(newCoordinatesLq)>10^-5); 
        iteration(FSI_curIT).iterations_error = max(abs((newCoordinatesLqConv - oldCoordinatesLqConv)./newCoordinatesLqConv));
    end
    if restart==0 && FSI_curIT == 1
        iteration(FSI_curIT).iterations_error = inf;                    % We just started the simulation
    end
    if FSI_curIT > 1
        % don't consider values too close to zero for convergence
        newCoordinatesLqConv = newCoordinatesLq(abs(newCoordinatesLq)>10^-5); 
        oldCoordinatesLqConv = oldCoordinatesLq(abs(newCoordinatesLq)>10^-5); 
        iteration(FSI_curIT).iterations_error = max(abs((newCoordinatesLqConv - oldCoordinatesLqConv)./newCoordinatesLqConv));
    end
    oldCoordinatesLq = newCoordinatesLq;
    
    
    %% Run Unsteady Panel Method
    
    % thin plate spline (TPS) INTERPOLATION
    size1U = analysisParams.n_seg_PM*2 + 1;
    size2U = analysisParams.num_airfoil_nodes_panel;
    
    TPSU_X = paramFSI.TPSU_X*xi1;
    TPSU_Y = paramFSI.TPSU_Y*xi1;
    TPSU_Z = paramFSI.TPSU_Z*xi1;
    TPSD_X = paramFSI.TPSD_X*xi1;
    TPSD_Y = paramFSI.TPSD_Y*xi1;
    TPSD_Z = paramFSI.TPSD_Z*xi1;
    
    xTPSsU = reshape(TPSU_X, size1U, size2U);
    yTPSsU = reshape(TPSU_Y, size1U, size2U);
    zTPSsU = reshape(TPSU_Z, size1U, size2U);
    xTPSsD = reshape(TPSD_X, size1U, size2U);
    yTPSsD = reshape(TPSD_Y, size1U, size2U);
    zTPSsD = reshape(TPSD_Z, size1U, size2U);
    
    yTPSU_tot = yMesh_loop(:,1:analysisParams.num_airfoil_nodes_panel) + flip(flip(yTPSsU,1),2);
    xTPSU_tot = xMesh_loop(:,1:analysisParams.num_airfoil_nodes_panel) + flip(flip(xTPSsU,1),2);
    zTPSU_tot = zMesh_loop(:,1:analysisParams.num_airfoil_nodes_panel) + flip(flip(zTPSsU,1),2);
    
    yTPSD_tot = yMesh_loop(:,analysisParams.num_airfoil_nodes_panel+1:end) + flip(flip(yTPSsD(:,2:end-1),1),2);
    xTPSD_tot = xMesh_loop(:,analysisParams.num_airfoil_nodes_panel+1:end) + flip(flip(xTPSsD(:,2:end-1),1),2);
    zTPSD_tot = zMesh_loop(:,analysisParams.num_airfoil_nodes_panel+1:end) + flip(flip(zTPSsD(:,2:end-1),1),2);
    
    xTPSs = [xTPSU_tot, xTPSD_tot];
    yTPSs = [yTPSU_tot, yTPSD_tot];
    zTPSs = [zTPSU_tot, zTPSD_tot];
    
    
    % TPS interpolation VELOCITY unsteady aero
    TPSUv2 = paramFSI.TPSU_Y*xdi1;
    TPSDv2 = paramFSI.TPSD_Y*xdi1;
    vTPSsU = reshape(TPSUv2, size1U, size2U);
    vTPSsD = reshape(TPSDv2, size1U, size2U);
    
    vTPSU_tot = flip(flip(vTPSsU,1),2);
    vTPSD_tot = flip(flip(vTPSsD(:,2:end-1),1),2);
    vTPSs = [vTPSU_tot, vTPSD_tot];
    
    aer_z_V = vTPSs;
    
    % Change wing reference system
    aer_x = xTPSs;
    aer_y = -zTPSs;
    aer_z = yTPSs;
    
    M1 = paramFSI.M1;
    N1a = paramFSI.N1;
    x = zeros(N1a+1,M1);
    y = zeros(N1a+1,M1);
    z = zeros(N1a+1,M1);
    x(1:N1a,1:M1) = aer_x';
    y(1:N1a,1:M1) = aer_y';
    z(1:N1a,1:M1) = aer_z';
    
    % close trailing edge
    N1 = N1a+1;
    x(end,:) = x(1,:);
    y(end,:) = y(1,:);
    z(end,:) = z(1,:);
    
    % Flip for the correct orientation
    x = flip(x);
    y = flip(y);
    z = flip(z);
    
    M = M1-1;
    N = N1-1;
    N2 = N+2;
    newWakeDistanceFactor = 0.25; % The new wake will extend up to 25% of the distance covered by the trailing edge in the time step. Refer to Katz for that
    
    x(N2,1:M1) = x(N1,1:M1)+ newWakeDistanceFactor*(vinf(1) - analysisParams.Omega(2)*z(N1,1:M1) + analysisParams.Omega(3)*y(N1,1:M1))*analysisParams.deltaT;
    y(N2,1:M1) = y(N1,1:M1);%+ newWakeDistanceFactor*(vinf(2) - analysisParams.Omega(3)*x(N1,1:M1) + analysisParams.Omega(1)*z(N1,1:M1))*analysisParams.deltaT;
    z(N2,1:M1) = z(N1,1:M1);%+ newWakeDistanceFactor*(vinf(3) - analysisParams.Omega(1)*y(N1,1:M1) + analysisParams.Omega(2)*x(N1,1:M1))*analysisParams.deltaT;
    
    
    %% Here we add all th other x,y,z's due to the continuos wake shed
    
    x = [x;Wake_Information.X];
    y = [y;Wake_Information.Y];
    z = [z;Wake_Information.Z];
    
    % We also compute the number of wake panels
    Nwakes = cutWakeAt-1;
    
    % Compute the aero grid
    Ax = (x(2:end,2:end)-x(1:(end-1),1:(end-1)));
    Ay = (y(2:end,2:end)-y(1:(end-1),1:(end-1)));
    Az = (z(2:end,2:end)-z(1:(end-1),1:(end-1)));
    Bx = (x(2:end,1:(end-1))-x(1:(end-1),2:end));
    By = (y(2:end,1:(end-1))-y(1:(end-1),2:end));
    Bz = (z(2:end,1:(end-1))-z(1:(end-1),2:end));
    
    FF = farfieldfaktor*max(realsqrt(Ax.^2+Ay.^2+Az.^2),realsqrt(Bx.^2+By.^2+Bz.^2));
    
    C1 = Ay.*Bz-Az.*By;
    C2 = Az.*Bx-Ax.*Bz;
    C3 = Ax.*By-Ay.*Bx;
    
    modc = sqrt(C1.^2 + C2.^2 + C3.^2);
    S = modc/2;
    modc(modc==0) = 10^-20;
    n1 = C1./modc;
    n2 = C2./modc;
    n3 = C3./modc;
    
    % calculation of colocation points
    cx = (x(1:(end-1),1:(end-1)) + x(1:(end-1),2:end) + x(2:end,1:(end-1)) + x(2:end,2:end))/4;
    cy = (y(1:(end-1),1:(end-1)) + y(1:(end-1),2:end) + y(2:end,1:(end-1)) + y(2:end,2:end))/4;
    cz = (z(1:(end-1),1:(end-1)) + z(1:(end-1),2:end) + z(2:end,1:(end-1)) + z(2:end,2:end))/4;
    
    % calculation of u(longitudinal), p(transversal) and o(perpendicular) vectors
    ux = ((x(2:end,1:(end-1))+x(2:end,2:end))-(x(1:(end-1),1:(end-1))+x(1:(end-1),2:end)))/2;
    uy = ((y(2:end,1:(end-1))+y(2:end,2:end))-(y(1:(end-1),1:(end-1))+y(1:(end-1),2:end)))/2;
    uz = ((z(2:end,1:(end-1))+z(2:end,2:end))-(z(1:(end-1),1:(end-1))+z(1:(end-1),2:end)))/2;
    uu = sqrt(ux.^2 + uy.^2 + uz.^2);
    u1 = ux./uu;
    u2 = uy./uu;
    u3 = uz./uu;
    px = -((x(1:(end-1),2:end)+x(2:end,2:end))-(x(1:(end-1),1:(end-1))+x(2:end,1:(end-1))))/2;
    py = -((y(1:(end-1),2:end)+y(2:end,2:end))-(y(1:(end-1),1:(end-1))+y(2:end,1:(end-1))))/2;
    pz = -((z(1:(end-1),2:end)+z(2:end,2:end))-(z(1:(end-1),1:(end-1))+z(2:end,1:(end-1))))/2;
    pp = sqrt(px.^2 + py.^2 + pz.^2);
    p1 = px./pp;
    p2 = py./pp;
    p3 = pz./pp;
    o1 = n2.*u3 - n3.*u2;
    o2 = u1.*n3 - n1.*u3;
    o3 = n1.*u2 - n2.*u1;
    
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

    % calculation of rotational velocity for unsteady aero: Omega x r on
    % the surface of the body
    OxR1 = analysisParams.Omega(2)*cz(1:N,1:M) - analysisParams.Omega(3)*cy(1:N,1:M);
    OxR2 = analysisParams.Omega(3)*cx(1:N,1:M) - analysisParams.Omega(1)*cz(1:N,1:M);
    OxR3 = analysisParams.Omega(1)*cy(1:N,1:M) - analysisParams.Omega(2)*cx(1:N,1:M);
    % caclulation of structural velocity for unsteady aero
    vUSi = zeros(N1,M1);
    vUSi(1:N1a,1:M1) = aer_z_V';
    vUSi = flip(vUSi);
    vUS = (vUSi(1:(end-1),1:(end-1)) + vUSi(1:(end-1),2:end) + vUSi(2:end,1:(end-1)) + vUSi(2:end,2:end))/4;
    
    if abs(max(max(vUS(:,1:M/2)))) > abs(min(min(vUS(:,1:M/2))))
        wing_def.vUSmaxL = max(max(vUS(:,1:M/2)));
    else
        wing_def.vUSmaxL = min(min(vUS(:,1:M/2)));
    end
    if abs(max(max(vUS(:,M/2+1:end)))) > abs(min(min(vUS(:,M/2+1:end))))
        wing_def.vUSmaxR = max(max(vUS(:,M/2+1:end)));
    else
        wing_def.vUSmaxR = min(min(vUS(:,M/2+1:end)));
    end
    
    vinf1 = vinf(1) - OxR1;
    vinf2 = vinf(2) - OxR2;
    vinf3 = vinf(3) - OxR3 - vUS;
    
    % setting source strength
    sigmaq = -(n1(1:N,1:M).*vinf1 + n2(1:N,1:M).*vinf2 + n3(1:N,1:M).*vinf3); % Omega x r
    sigma = reshape(sigmaq',N*M,1);
    
    
    [AM,BM]=influence_self_lift(FF,error,S,M,N,N1,M*N,u1,u2,u3,o1,o2,o3,n1,n2,n3,cx,cy,cz,x1,y1,x2,y2,x3,y3,x4,y4,d1,d2,d3,d4);
    
    % We now need to compute also the influence coefficients for the other wake panels
    if Nwakes~=0
        [CM]=influence_wake(FF,error,S,M,N,Nwakes,u1,u2,u3,o1,o2,o3,n1,n2,n3,cx,cy,cz,x1,y1,x2,y2,x3,y3,x4,y4,d1,d2,d3,d4);
        RHS = -BM*sigma-CM*Wake_Information.Strength;
    else
        RHS = -BM*sigma;
    end
    gamma = AM\RHS;
    
    % returning dipole values to i,j matrix
    gamma_bounded = zeros(N,M);
    K = 0;
    for k = 1:N
        for l = 1:M
            K = K+1;
            gamma_bounded(k,l) = gamma(K);
        end
    end
    
    gamma_FirstWake = (gamma_bounded(end,:)-gamma_bounded(1,:));
    gamma_FirstWake = gamma_FirstWake(:);
    
    % calculation of speeds, pressures and forces
    q = 0.5*analysisParams.rho*V^2;
    
    % calculation of induced speeds N-direction
    dx = zeros(N,M);
    dy = zeros(N,M);
    dz = zeros(N,M);
    qu = zeros(N,M);
    
    % boundary panels - first order interpolation
    dx(1,:) = cx(2,:)-cx(1,:);
    dy(1,:) = cy(2,:)-cy(1,:);
    dz(1,:) = cz(2,:)-cz(1,:);
    qu(1,:) = gamma_bounded(2,:)-gamma_bounded(1,:);
    dx(N,:) = cx(N,:)-cx(N-1,:);
    dy(N,:) = cy(N,:)-cy(N-1,:);
    dz(N,:) = cz(N,:)-cz(N-1,:);
    qu(N,:) = gamma_bounded(N,:)-gamma_bounded(N-1,:);
    
    % longitudinal
    dx(2:N-1,:) = cx(3:N,:)-cx(1:N-2,:);
    dy(2:N-1,:) = cy(3:N,:)-cy(1:N-2,:);
    dz(2:N-1,:) = cz(3:N,:)-cz(1:N-2,:);
    qu(2:N-1,:) = gamma_bounded(3:N,:)-gamma_bounded(1:N-2,:);
    
    dru = sqrt(dx.^2+dy.^2+dz.^2);
    qu = qu./dru;
    
    % calculation of induced speeds M-direction
    dx = zeros(N,M);
    dy = zeros(N,M);
    dz = zeros(N,M);
    qp = zeros(N,M);
    
    % boundary panels - first order interpolation
    dx(:,1) = cx(1:N,2)-cx(1:N,1);
    dy(:,1) = cy(1:N,2)-cy(1:N,1);
    dz(:,1) = cz(1:N,2)-cz(1:N,1);
    qp(:,1) = gamma_bounded(:,2)-gamma_bounded(:,1);
    dx(:,M) = cx(1:N,M)-cx(1:N,M-1);
    dy(:,M) = cy(1:N,M)-cy(1:N,M-1);
    dz(:,M) = cz(1:N,M)-cz(1:N,M-1);
    qp(:,M) = gamma_bounded(:,M)-gamma_bounded(:,M-1);
    
    % transversal
    dx(:,2:M-1) = cx(1:N,3:M)-cx(1:N,1:M-2);
    dy(:,2:M-1) = cy(1:N,3:M)-cy(1:N,1:M-2);
    dz(:,2:M-1) = cz(1:N,3:M)-cz(1:N,1:M-2);
    qp(:,2:M-1) = gamma_bounded(:,3:M)-gamma_bounded(:,1:M-2);
    
    drp = sqrt(dx.^2+dy.^2+dz.^2);
    qp = -qp./drp;
    
    qo = (p1(1:N,:).*o1(1:N,:) + p2(1:N,:).*o2(1:N,:) + p3(1:N,:).*o3(1:N,:)).*qp;
    
    gu = u1(1:N,:).*vinf1(1:N,:) + u2(1:N,:).*vinf2(1:N,:) + u3(1:N,:).*vinf3(1:N,:); % Omega x r
    go = o1(1:N,:).*vinf1(1:N,:) + o2(1:N,:).*vinf2(1:N,:) + o3(1:N,:).*vinf3(1:N,:); % Omega x r
    
    vx = (-qu + gu).*u1(1:N,:) + (-qo + go).*o1(1:N,:);
    vy = (-qu + gu).*u2(1:N,:) + (-qo + go).*o2(1:N,:);
    vz = (-qu + gu).*u3(1:N,:) + (-qo + go).*o3(1:N,:);
    v = sqrt(vx.^2+vy.^2+vz.^2);
    
    % We use simple backward differences for the time derivative of the
    % pressure coefficient cp
    cp = 1 - v.^2/V^2 + 2/V^2*(gamma_bounded-Wake_Information.Gamma_Old)./analysisParams.deltaT;      
    
    
    % initial conditions in case this FSI run generates data for the calculation of the empirical observaibility gramian
    if paramFSI.obsGramian.run
        if paramFSI.obsGramian.nState <= nStatesAero && FSI_curIT == 1
            cpRS = reshape(cp',N*M,1);
            cpRS(paramFSI.obsGramian.nState) = cpRS(paramFSI.obsGramian.nState)*paramFSI.obsGramian.InpAmp(2); % 1.05
            cpRS_0 = cpRS(paramFSI.obsGramian.nState);
            cpNew = reshape(cpRS,M,N)';
            cp = cpNew;
        elseif paramFSI.obsGramian.nState <= nStatesAero && FSI_curIT > 1
            cpRS = reshape(cp',N*M,1);
            cpRS(paramFSI.obsGramian.nState) = cpRS_0;
            cpNew = reshape(cpRS,M,N)';
            cp = cpNew;
        end
    end
    
    dX = -cp.*S(1:N,:).*n1(1:N,:)*q;
    dY = -cp.*S(1:N,:).*n2(1:N,:)*q;
    dZ = -cp.*S(1:N,:).*n3(1:N,:)*q;
    
    FX = sum(sum(dX));
    FZ = sum(sum(dZ));
    FL = sum(sum(-dY.*cz(1:N,:)+dZ.*cy(1:N,:)));
    
    FM = sum(sum( dX.*cz(1:N,:)-dZ.*(cx(1:N,:)+analysisParams.cM_pos(1))));
    
    FN = sum(sum(-dX.*cy(1:N,:)+dY.*cx(1:N,:)));
    
    dXU = dX(size(dX,1)/2:end,:);
    dYU = dY(size(dX,1)/2:end,:);
    dZU = dZ(size(dX,1)/2:end,:);
    dXD = dX(1:size(dX,1)/2+1,:);
    dYD = dY(1:size(dX,1)/2+1,:);
    dZD = dZ(1:size(dX,1)/2+1,:);
   
    
    % inverse distance weighting interpolation
    FAeroUIDW = [dXU(:),dZU(:),dYU(:)];
    FAeroDIDW = [dXD(:),dZD(:),dYD(:)];
    
    FStructUIDW_X = paramFSI.IDWuX*FAeroUIDW(:,1);
    FStructUIDW_Y = paramFSI.IDWuY*FAeroUIDW(:,2);
    FStructUIDW_Z = paramFSI.IDWuZ*FAeroUIDW(:,3);
    
    FStructDIDW_X = paramFSI.IDWdX*FAeroDIDW(:,1);
    FStructDIDW_Y = paramFSI.IDWdY*FAeroDIDW(:,2);
    FStructDIDW_Z = paramFSI.IDWdZ*FAeroDIDW(:,3);
    
    fiT = FStructUIDW_X + FStructUIDW_Y + FStructUIDW_Z + FStructDIDW_X + FStructDIDW_Y + FStructDIDW_Z;

    % calculate aerodynamic coefficients
    Sw = meanChord*b;
    FZA = FZ*cos(alpha)-FX*sin(alpha);
    wing_def.cL = FZA/q/Sw;
    wing_def.cRoll = FL/q/Sw/b;
    wing_def.cPitch = FM/q/Sw/meanChord;
    wing_def.cYaw = FN/q/Sw/b;

    
    %% Convergency check on the pressure distribution
    cP_conv_new_up = dZU;
    cP_conv_new_dn = dZD;
    
    if restart==0 && FSI_curIT ==1
        iteration(FSI_curIT).iterations_error_cP = inf;                 % We just started the simulation
        iteration(FSI_curIT).iterations_error_cPrel = inf;
    else
        iteration(FSI_curIT).iterations_error_cP = norm((cP_conv_new_up(2:end-1,:) - flipud(cP_conv_new_dn(2:end-1,:)) - (cP_conv_old_up(2:end-1,:) - flipud(cP_conv_old_dn(2:end-1,:)))),2); % do not compare cP on leading and trailing edge...
        convUp = max(max(abs((cP_conv_new_up(2:end-1,:) - cP_conv_old_up(2:end-1,:))./cP_conv_new_up(2:end-1,:)))); % do not compare cP on leading and trailing edge...
        convDn = max(max(abs((cP_conv_new_dn(2:end-1,:) - cP_conv_old_dn(2:end-1,:))./cP_conv_new_dn(2:end-1,:)))); % do not compare cP on leading and trailing edge...
        iteration(FSI_curIT).iterations_error_cPrel = max([convUp,convDn]); % do not compare cP on leading and trailing edge...
    end
    if convergencyCheck(iteration(FSI_curIT), analysisParams)
        wing_def.status = 3;
        wing_def.oldCoordinatesL = newCoordinatesLq;
        wing_def.cP_conv_old_up = cP_conv_new_up;
        wing_def.cP_conv_old_dn = cP_conv_new_dn;
        break
    end
    
    if FSI_curIT > 100
        wing_def.oldCoordinatesL = newCoordinatesLq;
        wing_def.cP_conv_old_up = cP_conv_new_up;
        wing_def.cP_conv_old_dn = cP_conv_new_dn;
        break
    end
    
    cP_conv_old_up = cP_conv_new_up;
    cP_conv_old_dn = cP_conv_new_dn;
    
    % convergence in first loop for aerodynamics only simulation
    if paramFSI.aeroOnly
        wing_def.oldCoordinatesL = newCoordinatesLq;
        wing_def.cP_conv_old_up = cP_conv_new_up;
        wing_def.cP_conv_old_dn = cP_conv_new_dn;
        break
    end
    
end


%% recalculate induced drag with ellt
CX = sum(dX,1)./sum(S(1:N,:),1)/q*2;
CZ = sum(dZ,1)./sum(S(1:N,:),1)/q*2;

cLdistr_val = CZ*cos(alpha)-CX*sin(alpha);
cLdistr_val = cLdistr_val';

cLdistr_val_interp = interp1(cLdistr_yPos, cLdistr_val, y_segMidpoint, 'linear', 'extrap');

K_alpha_ind_ELLT = paramFSI.K_alpha_ind_ELLT;
[eLLT_cDi,wing_def.cYaw,wing_def.FN] = eLLT_PM_cDi(y_segLimit, meanChord, V, cLdistr_val_interp, K_alpha_ind_ELLT,q,Sw,b);

% Overwrite panel method induced drag with the one calculated by means of the lifting line
wing_def.cDi = eLLT_cDi;
wing_def.FXA = eLLT_cDi*q*Sw;
wing_def.cLroot = cLdistr_val(M/2);


%% calculate viscous drag with XFOIL precalulated cDvisc
Re = paramFSI.viscPre.chordMean*analysisParams.rho*V/mu;

try
    if Re > max(max(paramFSI.viscPre.ReINT))
        wing_def.cDvisc = interp1(paramFSI.viscPre.alpha(2,:),paramFSI.viscPre.cd(2,:),alpha*180/pi);
    elseif Re < max(max(paramFSI.viscPre.ReINT))
        wing_def.cDvisc = interp1(paramFSI.viscPre.alpha(1,:),paramFSI.viscPre.cd(1,:),alpha*180/pi);
    else
        wing_def.cDvisc = interp2(paramFSI.viscPre.alpha,paramFSI.viscPre.ReINT,paramFSI.viscPre.cd,alpha*180/pi,Re);
    end
catch
    wing_def.cDvisc = max(max(paramFSI.viscPre.cd));
end
if isnan(wing_def.cDvisc)
    wing_def.cDvisc = max(max(paramFSI.viscPre.cd));
end


%% output modal amplitudes, velocities, and accelerations
% NM-alpha
wing_def.NM = NM;
wing_def.NM.xi = xi1;
wing_def.NM.xdi = xdi1; % x dot i
wing_def.NM.xddi = xddi1; % x dot dot i

wing_def.NM.fi = fi1; % f i
if ~paramFSI.aeroOnly
    wing_def.fiT = fiT;
else
    wing_def.fiT = fiT*0;
end
wing_def.fiACT = fiACT;

wing_def.Snapshot = [reshape(cp',N*M,1); xi1];

wing_def.FSI_curIT = FSI_curIT;


%% Treat wake rollup and the necessary output for the next step

wing_def.Wake_Information.Strength = [gamma_FirstWake ; Wake_Information.Strength];

wing_def.Wake_Information.Gamma_Old = gamma_bounded;



%%% For the moment, debugging reasons, just a flat wake -> to be updated
% wing_def.Wake_Information.X = x(N2:end,:)+vinf(1)*analysisParams.deltaT;
% wing_def.Wake_Information.Y = y(N2:end,:);%+vinf(2)*analysisParams.deltaT;
% wing_def.Wake_Information.Z = z(N2:end,:);%+vinf(3)*analysisParams.deltaT;


%%% Here we compute the real rollup

%We use the RK2 method because of stability reasons (Bhagart and Leishman)

x_to_move = zeros((M+1)*(Nwakes+1),1);
y_to_move = zeros((M+1)*(Nwakes+1),1);
z_to_move = zeros((M+1)*(Nwakes+1),1);
K = 0;
for k = N+2:N+2+Nwakes
    for l = 1:M+1
        K = K+1;
        x_to_move(K) = x(k,l);
        y_to_move(K) = y(k,l);
        z_to_move(K) = z(k,l);
    end
end

%Predict
[us,vs,ws,ud,vd,wd]=induced_velocity(FF,error,S,M,N,Nwakes,u1,u2,u3,o1,o2,o3,n1,n2,n3,cx,cy,cz,x1,y1,x2,y2,x3,y3,x4,y4,d1,d2,d3,d4,x,y,z,analysisParams.deltaT,'predict');

OxR1w = analysisParams.Omega(2)*z(N+2:N+2+Nwakes,1:M+1) - analysisParams.Omega(3)*y(N+2:N+2+Nwakes,1:M+1);
OxR2w = analysisParams.Omega(3)*x(N+2:N+2+Nwakes,1:M+1) - analysisParams.Omega(1)*z(N+2:N+2+Nwakes,1:M+1);
OxR3w = analysisParams.Omega(1)*y(N+2:N+2+Nwakes,1:M+1) - analysisParams.Omega(2)*x(N+2:N+2+Nwakes,1:M+1);
OxR1w = reshape(OxR1w',(Nwakes+1)*(M+1),1);
OxR2w = reshape(OxR2w',(Nwakes+1)*(M+1),1);
OxR3w = reshape(OxR3w',(Nwakes+1)*(M+1),1);

x_moved = x_to_move + (us*sigma + ud*[gamma ; wing_def.Wake_Information.Strength] + vinf(1) - OxR1w)* analysisParams.deltaT;
y_moved = y_to_move + (vs*sigma + vd*[gamma ; wing_def.Wake_Information.Strength] + vinf(2) - OxR2w)* analysisParams.deltaT;
z_moved = z_to_move + (ws*sigma + wd*[gamma ; wing_def.Wake_Information.Strength] + vinf(3) - OxR3w)* analysisParams.deltaT;

X_new = zeros(Nwakes+2,M+1);
Y_new = X_new;
Z_new = Y_new;

K = 0;
for k = N+2:N+2+Nwakes
    i = k-N;
    for l = 1:M+1
        K = K+1;
        X_new(i,l) = x_moved(K);
        Y_new(i,l) = y_moved(K);
        Z_new(i,l) = z_moved(K);
    end
end
X_new(1,:) = x(N1,:);
Y_new(1,:) = y(N1,:);         % We need all the wake panel edges, so we recover the trailing edge
Z_new(1,:) = z(N1,:);
clear x_moved y_moved z_moved

Ax_new = (X_new(2:end,2:end)-X_new(1:(end-1),1:(end-1)));
Ay_new = (Y_new(2:end,2:end)-Y_new(1:(end-1),1:(end-1)));
Az_new = (Z_new(2:end,2:end)-Z_new(1:(end-1),1:(end-1)));
Bx_new = (X_new(2:end,1:(end-1))-X_new(1:(end-1),2:end));
By_new = (Y_new(2:end,1:(end-1))-Y_new(1:(end-1),2:end));
Bz_new = (Z_new(2:end,1:(end-1))-Z_new(1:(end-1),2:end));

FF_new = farfieldfaktor*max(realsqrt(Ax_new.^2+Ay_new.^2+Az_new.^2),realsqrt(Bx_new.^2+By_new.^2+Bz_new.^2));

C1_new = Ay_new.*Bz_new-Az_new.*By_new;
C2_new = Az_new.*Bx_new-Ax_new.*Bz_new;
C3_new = Ax_new.*By_new-Ay_new.*Bx_new;

modc_new = sqrt(C1_new.^2 + C2_new.^2 + C3_new.^2);
S_new = modc_new/2;
modc_new(modc_new==0) = 10^-20;
n1_new = C1_new./modc_new;
n2_new = C2_new./modc_new;
n3_new = C3_new./modc_new;

% calculation of colocation points
cx_new = (X_new(1:(end-1),1:(end-1)) + X_new(1:(end-1),2:end) + X_new(2:end,1:(end-1)) + X_new(2:end,2:end))/4;
cy_new = (Y_new(1:(end-1),1:(end-1)) + Y_new(1:(end-1),2:end) + Y_new(2:end,1:(end-1)) + Y_new(2:end,2:end))/4;
cz_new = (Z_new(1:(end-1),1:(end-1)) + Z_new(1:(end-1),2:end) + Z_new(2:end,1:(end-1)) + Z_new(2:end,2:end))/4;
% calculation of u(longitudinal), p(transversal) and o(perpendicular) vectors
ux_new = ((X_new(2:end,1:(end-1))+X_new(2:end,2:end))-(X_new(1:(end-1),1:(end-1))+X_new(1:(end-1),2:end)))/2;
uy_new = ((Y_new(2:end,1:(end-1))+Y_new(2:end,2:end))-(Y_new(1:(end-1),1:(end-1))+Y_new(1:(end-1),2:end)))/2;
uz_new = ((Z_new(2:end,1:(end-1))+Z_new(2:end,2:end))-(Z_new(1:(end-1),1:(end-1))+Z_new(1:(end-1),2:end)))/2;
uu_new = sqrt(ux_new.^2 + uy_new.^2 + uz_new.^2);
u1_new = ux_new./uu_new;
u2_new = uy_new./uu_new;
u3_new = uz_new./uu_new;
px_new = -((X_new(1:(end-1),2:end)+X_new(2:end,2:end))-(X_new(1:(end-1),1:(end-1))+X_new(2:end,1:(end-1))))/2;
py_new = -((Y_new(1:(end-1),2:end)+Y_new(2:end,2:end))-(Y_new(1:(end-1),1:(end-1))+Y_new(2:end,1:(end-1))))/2;
pz_new = -((Z_new(1:(end-1),2:end)+Z_new(2:end,2:end))-(Z_new(1:(end-1),1:(end-1))+Z_new(2:end,1:(end-1))))/2;
pp_new = sqrt(px_new.^2 + py_new.^2 + pz_new.^2);
p1_new = px_new./pp_new;
p2_new = py_new./pp_new;
p3_new = pz_new./pp_new;
o1_new = n2_new.*u3_new - n3_new.*u2_new;
o2_new = u1_new.*n3_new - n1_new.*u3_new;
o3_new = n1_new.*u2_new - n2_new.*u1_new;


% Clockwise in order to have a positive doublet in the normal direction
x1_new = (X_new(1:end-1,2:end)-cx_new).*u1_new+(Y_new(1:end-1,2:end)-cy_new).*u2_new+(Z_new(1:end-1,2:end)-cz_new).*u3_new;
y1_new = (X_new(1:end-1,2:end)-cx_new).*o1_new+(Y_new(1:end-1,2:end)-cy_new).*o2_new+(Z_new(1:end-1,2:end)-cz_new).*o3_new;
x2_new = (X_new(1:end-1,1:end-1)-cx_new).*u1_new+(Y_new(1:end-1,1:end-1)-cy_new).*u2_new+(Z_new(1:end-1,1:end-1)-cz_new).*u3_new;
y2_new = (X_new(1:end-1,1:end-1)-cx_new).*o1_new+(Y_new(1:end-1,1:end-1)-cy_new).*o2_new+(Z_new(1:end-1,1:end-1)-cz_new).*o3_new;
x3_new = (X_new(2:end,1:end-1)-cx_new).*u1_new+(Y_new(2:end,1:end-1)-cy_new).*u2_new+(Z_new(2:end,1:end-1)-cz_new).*u3_new;
y3_new = (X_new(2:end,1:end-1)-cx_new).*o1_new+(Y_new(2:end,1:end-1)-cy_new).*o2_new+(Z_new(2:end,1:end-1)-cz_new).*o3_new;
x4_new = (X_new(2:end,2:end)-cx_new).*u1_new+(Y_new(2:end,2:end)-cy_new).*u2_new+(Z_new(2:end,2:end)-cz_new).*u3_new;
y4_new = (X_new(2:end,2:end)-cx_new).*o1_new+(Y_new(2:end,2:end)-cy_new).*o2_new+(Z_new(2:end,2:end)-cz_new).*o3_new;

d1_new = realsqrt((x2_new-x1_new).^2 + (y2_new-y1_new).^2);
d2_new = realsqrt((x3_new-x2_new).^2 + (y3_new-y2_new).^2);
d3_new = realsqrt((x4_new-x3_new).^2 + (y4_new-y3_new).^2);
d4_new = realsqrt((x1_new-x4_new).^2 + (y1_new-y4_new).^2);

% Recompute
[us_new,vs_new,ws_new,ud_new,vd_new,wd_new]=induced_velocity([FF(1:N,:);FF_new],error,[S(1:N,:);S_new],M,N,Nwakes,[u1(1:N,:);u1_new],[u2(1:N,:);u2_new]...
    ,[u3(1:N,:);u3_new],[o1(1:N,:);o1_new],[o2(1:N,:);o2_new],[o3(1:N,:);o3_new],[n1(1:N,:);n1_new],...
    [n2(1:N,:);n2_new],[n3(1:N,:);n3_new],[cx(1:N,:);cx_new],[cy(1:N,:);cy_new],[cz(1:N,:);cz_new],[x1(1:N,:);x1_new]...
    ,[y1(1:N,:);y1_new],[x2(1:N,:);x2_new],[y2(1:N,:);y2_new],[x3(1:N,:);x3_new],[y3(1:N,:);y3_new],[x4(1:N,:);x4_new]...
    ,[y4(1:N,:);y4_new],[d1(1:N,:);d1_new],[d2(1:N,:);d2_new],[d3(1:N,:);d3_new],[d4(1:N,:);d4_new],...
    [x(1:N,:);X_new],[y(1:N,:);Y_new],[z(1:N,:);Z_new],analysisParams.deltaT,'compute');

% Correct
U_corr = 0.5*((us*sigma + ud*[gamma ; wing_def.Wake_Information.Strength] + vinf(1) - OxR1w)+...
    (us_new*sigma + ud_new*[gamma ; wing_def.Wake_Information.Strength] + vinf(1) - OxR1w));
V_corr = 0.5*((vs*sigma + vd*[gamma ; wing_def.Wake_Information.Strength] + vinf(2) - OxR2w)+...
    (vs_new*sigma + vd_new*[gamma ; wing_def.Wake_Information.Strength] + vinf(2) - OxR2w));
W_corr = 0.5*((ws*sigma + wd*[gamma ; wing_def.Wake_Information.Strength] + vinf(3) - OxR3w)+...
    (ws_new*sigma + wd_new*[gamma ; wing_def.Wake_Information.Strength] + vinf(3) - OxR3w));


x_moved = x_to_move + U_corr*analysisParams.deltaT;
y_moved = y_to_move + V_corr*analysisParams.deltaT;
z_moved = z_to_move + W_corr*analysisParams.deltaT;


K = 0;
for k = N+2:N+2+Nwakes
    i = k-N-1;
    for l = 1:M+1
        K = K+1;
        wing_def.Wake_Information.X(i,l) = x_moved(K);
        wing_def.Wake_Information.Y(i,l) = y_moved(K);
        wing_def.Wake_Information.Z(i,l) = z_moved(K);
    end
end


if 0
    % We plot here the shedding procedure, for duoble checking. Comment
    % this about for increased performances
    figure(101)
    surf(cx(1:N,1:M),cy(1:N,1:M),cz(1:N,1:M),cp(1:N,1:M)) % cp
    hold on
    % axis equal
    axis([-0.5 0.5 -2.5 2.5 -0.5 0.5])
    %mesh(x(N1:end,:),y(N1:end,:),z(N1:end,:))
    hold off
    figure(102)
    plot(x(1:N,M/2+1),(cp(:,M/2)+cp(:,M/2))*0.5)
    
    figure(103)
    mesh([x(1:N+1,:);wing_def.Wake_Information.X],[y(1:N+1,:);wing_def.Wake_Information.Y],[z(1:N+1,:);wing_def.Wake_Information.Z])
    axis equal
end


%% Here we save the snapshot of the state in order to use DMD
if createROM
    Snapshot = [reshape(cp',N*M,1); xi1];  % The states are: the pressure coefficients and the modal amplitudes
    Snapshot_input = [ forceACTinp.L ; forceACTinp.fs ; analysisParams.Omega(1) ; analysisParams.Omega(2) ; analysisParams.Omega(3) ; alpha ] - inputSnapshot0';

    fid = fopen(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,'ROM',filesep,'V',num2str(V),filesep,'Snapshot',filesep,'Snapshot.csv'), 'a');
    if fid == -1
        mkdir(strcat(pwd,filesep,'data'));
        mkdir(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,'ROM',filesep,'V',num2str(V),filesep,'Snapshot'));
        fid = fopen(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,...
            'ROM',filesep,'V',num2str(V),filesep,'Snapshot',filesep,'Snapshot.csv'), 'a');
    end
    fprintf(fid, strcat('%e',repmat(', %e',1,length(Snapshot)-1),'\n'),Snapshot');
    fclose(fid);

    fid = fopen(strcat(pwd,filesep,'data',filesep,paramFSI.wingParams.airfoil,filesep,'ROM',filesep,...
        'V',num2str(V),filesep,'Snapshot',filesep,'Snapshot_input.csv'), 'a');
    fprintf(fid, strcat('%e',repmat(', %e',1,length(Snapshot_input)-1),'\n'),Snapshot_input');
    fclose(fid);
end

end