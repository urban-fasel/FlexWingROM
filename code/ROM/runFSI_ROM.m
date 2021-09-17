% Function that simulates any of the 3 ROMs models (specified with "type_string")

function simOutROM = runFSI_ROM(Velocity, ROMs, dT, iTest, simInput, type_string)


paramFSI = simInput.paramFSI;

nModes = paramFSI.nModes;

%% Simulation parameters

paramFSI.HHT.xi = zeros(nModes,1);
paramFSI.HHT.xdi = zeros(nModes,1); % x dot i
paramFSI.HHT.xddi = zeros(nModes,1); % x dot dot i
paramFSI.HHT.fi = zeros(nModes,1); % f i
paramFSI.HHT.act = [0,0];

% These are the outputs calculated with the nonlinear output function
% This funciton was initially introduced in the aDMDc work to calculate
% accurate nonlinear outputs, as we only linearized the singularity
% strengths and modal amplitudes of the model.
simOutROM.cL = zeros(iTest,1);
simOutROM.cD = zeros(iTest,1);
simOutROM.cRoll = zeros(iTest,1);
simOutROM.cYaw = zeros(iTest,1);
simOutROM.cPitch = zeros(iTest,1);
simOutROM.bendingModeAmplitude = zeros(iTest,1);

% These are the outputs corresponding to the linear input-output models.
% They are only available for the IODMD and the BMD model.
simOutROM.cL_states = zeros(iTest,1);
simOutROM.cD_states = zeros(iTest,1);
simOutROM.cRoll_states = zeros(iTest,1);
simOutROM.cYaw_states = zeros(iTest,1);
simOutROM.cPitch_states = zeros(iTest,1);
simOutROM.bendingModeAmplitude_states = zeros(iTest,1);

% set inputs
meanRotationRates = simInput.meanRotationRates;

inputMorphingSymmetric = simInput.inputMorphingSymmetric;
force_input = simInput.force_input;
Omega = simInput.Omega;
deltaAlpha = simInput.deltaAlpha;

rotOmega = simInput.rotOmega;
AlphaShape = simInput.AlphaShape;
forceShape1 = simInput.forceShape1;
forceShape2 = simInput.forceShape2;

VelocitySWEEP = linspace(min(Velocity),max(Velocity),iTest);

Input_ROMs0 = [0;0;0;0;0;simInput.alpha*pi/180] + simInput.inputSnapshot0';


%% run reduced order models
% Preprocessing
if(strcmp(type_string,'aDMDc'))
    Atil_mat = ROMs.Atil_mat; % r x r x n_g
    Btil_mat = ROMs.Btil_mat;
    Ftil_mat = ROMs.Ftil_mat;
    Uhat_mat = ROMs.Uhat_mat;
    Xmean = ROMs.Xmean;
    vv = ROMs.velocity_vec;
    n_g = length(vv);
    ROMs_aDMDc = cell(1,n_g);
    r = size(Atil_mat,1);
    for k = 1:n_g
        ROMin.Atil = Atil_mat(:,:,k);
        ROMin.Btil = Btil_mat(:,:,k);
        ROMin.Ftil = Ftil_mat(:,:,k);
        ROMin.Uhat = Uhat_mat(:,:,k);
        ROMin.Xmean = Xmean(:,k);
        ROMin.Velocity = vv(k);
        ROMin.State = zeros(r,1);
        ROMin.Input = zeros(6,1);
        ROMin.cpxiOld = ROMin.Xmean;

        ROMs_aDMDc{k} = ROMin; 
    end
    Input_ROMs0 = [0;0;0;0;0;simInput.alpha*pi/180]+simInput.inputSnapshot0';
    
elseif(strcmp(type_string,'aIODMD')||strcmp(type_string,'BMD'))
    [n_x, ~] = size(ROMs.Xmean);
    [r, n_u, n_g] = size(ROMs.G_mat);
    [n_y,~,~] = size(ROMs.H_mat);

    rho_vec = ROMs.velocity_vec;

    % Interpolation preparation
    V_F = zeros(length(rho_vec), r*r);
    V_G = zeros(length(rho_vec), r*n_u);
    V_H = zeros(length(rho_vec), n_y*r);
    V_D = zeros(length(rho_vec), n_y*n_u);
    V_L = zeros(length(rho_vec), r*n_u);
    V_E = zeros(length(rho_vec), n_y*n_u);
    V_Xmean = zeros(length(rho_vec), n_x);
    if(strcmp(type_string,'BMD'))
        V_W = zeros(length(rho_vec), n_x*r);
    end
    for k = 1:n_g
        V_F(k,:) = reshape(ROMs.F_mat(:,:,k),1,[]);
        V_G(k,:) = reshape(ROMs.G_mat(:,:,k),1,[]);
        V_H(k,:) = reshape(ROMs.H_mat(:,:,k),1,[]);
        V_D(k,:) = reshape(ROMs.D_mat(:,:,k),1,[]);
        V_L(k,:) = reshape(ROMs.L_mat(:,:,k),1,[]);
        V_E(k,:) = reshape(ROMs.E_mat(:,:,k),1,[]);
        V_Xmean(k,:) = ROMs.Xmean(:,k).';
        if(strcmp(type_string,'BMD'))
            V_W(k,:) = reshape(ROMs.W_mat(:,:,k),1,[]);
        end
    end

    if n_g == 1
        Xold = V_Xmean';
    else
        Xold = interp1(rho_vec.', V_Xmean, VelocitySWEEP(1)).'; % mean corresponding initial parameter state
    end
    Z = zeros(size(ROMs.F_mat,1),1); 
    Y = zeros(1,1);
    Input = zeros(6,1);
end


% run ROMs
VelocitySWEEP(iTest+1) = VelocitySWEEP(iTest);                                        
for i = 1:iTest
    if(strcmp(type_string,'aDMDc'))
        Input_ROMs = [force_input*forceShape1(i); inputMorphingSymmetric*forceShape2(i); meanRotationRates'+Omega'*rotOmega(i); deltaAlpha*AlphaShape(i)];
        ROMsInterp = [];
        
        for index = 1:length(ROMs_aDMDc)
            ROMs_aDMDc{index}.State = ROMs_aDMDc{index}.Atil*ROMs_aDMDc{index}.State + ROMs_aDMDc{index}.Btil*ROMs_aDMDc{index}.Input + ROMs_aDMDc{index}.Ftil*Input_ROMs;
            ROMs_aDMDc{index}.Input = Input_ROMs;
            
            cpxiInterp = ROMs_aDMDc{index}.Uhat*ROMs_aDMDc{index}.State + ROMs_aDMDc{index}.Xmean;
            cpxiPrime = (cpxiInterp - ROMs_aDMDc{index}.cpxiOld)/dT;
            ROMs_aDMDc{index}.cpxiOld = cpxiInterp;
            ROMsInterp = [ROMsInterp,[cpxiInterp;cpxiPrime(609:618)]]; % states: [cp, xi, xdoti] -> they are used in output function ObtainOutputROM
        end
        
        interpolation_method = 'spline';
        if n_g == 1
            XXX = ROMsInterp';
        else
            XXX = interp1(vv',ROMsInterp',VelocitySWEEP(i),interpolation_method, 'extrap');
        end
        simOutROM.bendingModeAmplitude(i) = XXX(610); % extract first bending mode 610 from X
        simOutROM.bendingModeAmplitude_states(i) = XXX(610); % just use same, since output doesn't exist

        Input_ROMsI = Input_ROMs + Input_ROMs0;
        out = ObtainOutputROM(paramFSI,Input_ROMsI,VelocitySWEEP(i),XXX');

        simOutROM.cL(i) = out(1);
        simOutROM.cD(i) = out(2);
        simOutROM.cRoll(i) = out(3);
        simOutROM.cPitch(i) = out(4);
        simOutROM.cYaw(i) = out(5);
        
    elseif(strcmp(type_string,'aIODMD')||strcmp(type_string,'BMD'))  
        % Interpolation
        rho_k = VelocitySWEEP(i);
        if n_g == 1
            V_F_q = V_F; 
            V_G_q = V_G;
            V_H_q = V_H;
            V_D_q = V_D;
            V_L_q = V_L;
            V_E_q = V_E;
            V_Xmean_q = V_Xmean;
            V_Xmean_q_plus = V_Xmean;
        else
            V_F_q = interp1(rho_vec.', V_F, rho_k); % , interpolation_method,'extrap'
            V_G_q = interp1(rho_vec.', V_G, rho_k);
            V_H_q = interp1(rho_vec.', V_H, rho_k);
            V_D_q = interp1(rho_vec.', V_D, rho_k);
            V_L_q = interp1(rho_vec.', V_L, rho_k);
            V_E_q = interp1(rho_vec.', V_E, rho_k);
            V_Xmean_q = interp1(rho_vec.', V_Xmean, rho_k);
            V_Xmean_q_plus = interp1(rho_vec.', V_Xmean, VelocitySWEEP(i+1));
        end
        
        if(strcmp(type_string,'BMD'))
            if n_g == 1
                V_W_q = V_W; 
                V_W_q_plus = V_W;
            else
                V_W_q = interp1(rho_vec.',V_W, rho_k); 
                V_W_q_plus = interp1(rho_vec.', V_W, VelocitySWEEP(i+1));
            end
        end

        if n_g == 1
            cL_bar = ROMs.cL_bar_vec;
            cD_bar = ROMs.cD_bar_vec;
            cRoll_bar = ROMs.cRoll_bar_vec;
            cPitch_bar = ROMs.cPitch_bar_vec;
            cYaw_bar = ROMs.cYaw_bar_vec;
            bendingModeAmplitude_bar = ROMs.bendingModeAmplitude_bar_vec;
        else
            cL_bar = interp1(rho_vec.', ROMs.cL_bar_vec, rho_k);
            cD_bar = interp1(rho_vec.', ROMs.cD_bar_vec, rho_k);
            cRoll_bar = interp1(rho_vec.', ROMs.cRoll_bar_vec, rho_k);
            cPitch_bar = interp1(rho_vec.', ROMs.cPitch_bar_vec, rho_k);
            cYaw_bar = interp1(rho_vec.', ROMs.cYaw_bar_vec, rho_k);
            bendingModeAmplitude_bar = interp1(rho_vec.', ROMs.bendingModeAmplitude_bar_vec, rho_k);
        end
        
        F = reshape(V_F_q, [r, r]);
        G = reshape(V_G_q, [r, n_u]);
        H = reshape(V_H_q, [n_y, r]);
        D = reshape(V_D_q, [n_y, n_u]);
        L = reshape(V_L_q, [r, n_u]);
        E = reshape(V_E_q, [n_y, n_u]);
        Xmean = V_Xmean_q.';
        Xmean_plus = V_Xmean_q_plus.';                                       

        Input_ROMs = [force_input*forceShape1(i); inputMorphingSymmetric*forceShape2(i); meanRotationRates'+Omega'*rotOmega(i); deltaAlpha*AlphaShape(i)];
        Y = H*Z + D*Input + E*Input_ROMs; % y_{k} = ... maybe change order
                                                        %     norm(ROMs.Q.'*Xmean - ROMs.Q.'*Xmean_plus,2)/norm(ROMs.Q.'*Xmean,2)
        if(strcmp(type_string,'BMD'))
            W_k = reshape(V_W_q, [n_x r]);
            W_k_plus = reshape(V_W_q_plus, [n_x r]);
            Z = F*Z + G*Input + L*Input_ROMs + W_k.'*Xmean - W_k_plus.'*Xmean_plus; % z_{k+1} = .. + ROMs.Q.'*Xmean - ROMs.Q.'*Xmean_plus
        elseif(strcmp(type_string,'aIODMD'))
            Z = F*Z + G*Input + L*Input_ROMs + ROMs.Q.'*Xmean - ROMs.Q.'*Xmean_plus;
        else
            disp('error 4!');
        end

        Input = Input_ROMs;
    
        if(strcmp(type_string,'BMD'))
            X = ROMs.V*Z + Xmean_plus;          
        elseif(strcmp(type_string,'aIODMD'))
            X = ROMs.Q*Z + Xmean_plus;
        end
        
        simOutROM.cL(i) = Y(1) + cL_bar;
        simOutROM.cD(i) = Y(2) + cD_bar;
        simOutROM.cRoll(i) = Y(3) + cRoll_bar;
        simOutROM.cPitch(i) = Y(4) + cPitch_bar;
        simOutROM.cYaw(i) = Y(5) + cYaw_bar;
        simOutROM.bendingModeAmplitude(i) = Y(6) + bendingModeAmplitude_bar; 

        cpxiPrime = (X - Xold)/dT;
        Xold = X;
        XXX = [X; cpxiPrime(609:618)]; % states: [cp, xi, xdoti] -> they are used in output function ObtainOutputROM
        
        Input_ROMsI = Input_ROMs + Input_ROMs0; 
        out = ObtainOutputROM(paramFSI,Input_ROMsI,VelocitySWEEP(i),XXX);

        simOutROM.cL_states(i) = out(1);
        simOutROM.cD_states(i) = out(2);
        simOutROM.cRoll_states(i) = out(3);
        simOutROM.cPitch_states(i) = out(4);
        simOutROM.cYaw_states(i) = out(5);
        simOutROM.bendingModeAmplitude_states(i) = X(610);

    else
        disp('error 5!');
    end
end


end