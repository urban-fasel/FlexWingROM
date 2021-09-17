function out = ObtainOutputROM(paramFSI,Input_ROMs,VEL,XXX)

      
    %% Recompute the physiscal quantitites starting from the states obtained in XXX

    AlphaShape = Input_ROMs(end);
        
    xi = XXX(end-(2*paramFSI.nModes-1):end-paramFSI.nModes);
  
    xMesh_loop = paramFSI.xMesh_loop;
    yMesh_loop = paramFSI.yMesh_loop;
    zMesh_loop = paramFSI.zMesh_loop;
    
    size1U = paramFSI.analysisParams.n_seg_PM*2 + 1;
    size2U = paramFSI.analysisParams.num_airfoil_nodes_panel;
    
    V = VEL;
    
    % TPS INTERPOLATION
    TPSU_X = paramFSI.TPSU_X*xi;
    TPSU_Y = paramFSI.TPSU_Y*xi;
    TPSU_Z = paramFSI.TPSU_Z*xi;
    TPSD_X = paramFSI.TPSD_X*xi;
    TPSD_Y = paramFSI.TPSD_Y*xi;
    TPSD_Z = paramFSI.TPSD_Z*xi;
    
    
    
    xTPSsU = reshape(TPSU_X, size1U, size2U);
    yTPSsU = reshape(TPSU_Y, size1U, size2U);
    zTPSsU = reshape(TPSU_Z, size1U, size2U);
    xTPSsD = reshape(TPSD_X, size1U, size2U);
    yTPSsD = reshape(TPSD_Y, size1U, size2U);
    zTPSsD = reshape(TPSD_Z, size1U, size2U);
    
    
    yTPSU_tot = yMesh_loop(:,1:paramFSI.analysisParams.num_airfoil_nodes_panel) + flip(flip(yTPSsU,1),2);
    xTPSU_tot = xMesh_loop(:,1:paramFSI.analysisParams.num_airfoil_nodes_panel) + flip(flip(xTPSsU,1),2);
    zTPSU_tot = zMesh_loop(:,1:paramFSI.analysisParams.num_airfoil_nodes_panel) + flip(flip(zTPSsU,1),2);
    
    yTPSD_tot = yMesh_loop(:,paramFSI.analysisParams.num_airfoil_nodes_panel+1:end) + flip(flip(yTPSsD(:,2:end-1),1),2);
    xTPSD_tot = xMesh_loop(:,paramFSI.analysisParams.num_airfoil_nodes_panel+1:end) + flip(flip(xTPSsD(:,2:end-1),1),2);
    zTPSD_tot = zMesh_loop(:,paramFSI.analysisParams.num_airfoil_nodes_panel+1:end) + flip(flip(zTPSsD(:,2:end-1),1),2);
    
    xTPSs = [xTPSU_tot, xTPSD_tot];
    yTPSs = [yTPSU_tot, yTPSD_tot];
    zTPSs = [zTPSU_tot, zTPSD_tot];

    
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
    
    Ax = (x(2:end,2:end)-x(1:(end-1),1:(end-1)));
    Ay = (y(2:end,2:end)-y(1:(end-1),1:(end-1)));
    Az = (z(2:end,2:end)-z(1:(end-1),1:(end-1)));
    Bx = (x(2:end,1:(end-1))-x(1:(end-1),2:end));
    By = (y(2:end,1:(end-1))-y(1:(end-1),2:end));
    Bz = (z(2:end,1:(end-1))-z(1:(end-1),2:end));
    
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

    % calculation of speeds, pressures and forces
    P_dyn = 0.5*paramFSI.analysisParams.rho*V^2;

    % cp
    X_to_be_analysed = XXX(1:608);
    k = 1;
    for i = 1:N
        for j = 1:M
            cp(i,j) = X_to_be_analysed(k);
            k = k+1;
        end
    end
    
    dX = -cp.*S(1:N,:).*n1(1:N,:)*P_dyn;
    dY = -cp.*S(1:N,:).*n2(1:N,:)*P_dyn;
    dZ = -cp.*S(1:N,:).*n3(1:N,:)*P_dyn;
    
    FX = sum(sum(dX));
    FZ = sum(sum(dZ));
    FL = sum(sum(-dY.*cz(1:N,:)+dZ.*cy(1:N,:)));
    
    FM = sum(sum( dX.*cz(1:N,:)-dZ.*(cx(1:N,:)+paramFSI.cog(1))));
     
    Sw = paramFSI.wingParams.chord*paramFSI.wingParams.b;
    FZA = FZ*cos(AlphaShape)-FX*sin(AlphaShape);
    
    cL = FZA/P_dyn/Sw;
    cRoll = FL/P_dyn/Sw/paramFSI.wingParams.b;
    cPitch = FM/P_dyn/Sw/paramFSI.wingParams.chord;
    
    
    %% recalculate induced drag with ellt
    CX = sum(dX,1)./sum(S(1:N,:),1)/P_dyn*2;
    CZ = sum(dZ,1)./sum(S(1:N,:),1)/P_dyn*2;
    
    cLdistr_val = CZ*cos(AlphaShape)-CX*sin(AlphaShape);
    cLdistr_val = cLdistr_val';
    
    cLdistr_val_interp = interp1(paramFSI.cLdistr_yPos, cLdistr_val, paramFSI.y_segMidpoint, 'linear', 'extrap');
   
    K_alpha_ind_ELLT = paramFSI.K_alpha_ind_ELLT;
    [cDi,cYaw,~] = eLLT_PM_cDi(paramFSI.y_segLimit, paramFSI.wingParams.chord, V, cLdistr_val_interp, K_alpha_ind_ELLT,P_dyn,Sw,paramFSI.wingParams.b);
    
    
    %% calculate viscous drag
    mu = 1.81*10^-5; % kg/(m·s)

    Re = paramFSI.viscPre.chordMean*paramFSI.analysisParams.rho*V/mu;
    try
        if Re > max(max(paramFSI.viscPre.ReINT))
            cDvisc = interp1(paramFSI.viscPre.alpha(2,:),paramFSI.viscPre.cd(2,:),AlphaShape*180/pi);
        elseif Re < max(max(paramFSI.viscPre.ReINT))
            cDvisc = interp1(paramFSI.viscPre.alpha(1,:),paramFSI.viscPre.cd(1,:),AlphaShape*180/pi);
        else
            cDvisc = interp2(paramFSI.viscPre.alpha,paramFSI.viscPre.ReINT,paramFSI.viscPre.cd,AlphaShape*180/pi,Re);
        end
    catch
        cDvisc = max(max(paramFSI.viscPre.cd));
    end
    if isnan(cDvisc)
        cDvisc = max(max(paramFSI.viscPre.cd));
    end
    
    
    cD = cDi + cDvisc;
    
    out = [cL, cD, cRoll, cPitch, cYaw];
end