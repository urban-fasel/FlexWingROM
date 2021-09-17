% Function that generates the BMD model for varying speed

function[BMD_Matrices] = Obtain_BMD_model(Snapshot, Snapshot_input, grid_velocity, type_string, rankMax, simInput)

[n_x, Nt, n_g] = size(Snapshot);
[n_u,~,~] = size(Snapshot_input);
n_s = Nt - 1;
n_y = 6; % number of outputs: cL, cD, cRoll, cPitch, cYaw, bending mode amplitude

% Reshaping the snapshots to be used for aIODMD
X_0 = zeros(n_x, n_s, n_g);
X_1 = zeros(n_x, n_s, n_g);
Y_0 = zeros(n_y, n_s, n_g);
U_0 = zeros(n_u, n_s, n_g);
U_1 = zeros(n_u, n_s, n_g);
Xmean = zeros(n_x, n_g);

cL_bar_vec = zeros(1, n_g);
cD_bar_vec = zeros(1, n_g);
cRoll_bar_vec = zeros(1, n_g);
cYaw_bar_vec = zeros(1, n_g);
cPitch_bar_vec = zeros(1, n_g);
bendingModeAmplitude_bar_vec = zeros(1, n_g);

% (empirical) gramian
L_c_mat = cell(1,n_g);
L_o_mat = cell(1,n_g);
smallest_size = inf;
for k = 1:n_g
    X_0(:,:,k) = Snapshot(:,1:end-1,k);
    X_1(:,:,k) = Snapshot(:,2:end,k);
    U_0(:,:,k) = Snapshot_input(:,1:end-1,k); % choose all snapshots or inputs
    U_1(:,:,k) = Snapshot_input(:,2:end,k); % choose all snapshots or inputs
    Xmean(:,k) = mean(X_0(:,:,k),2);
    
    V0 = grid_velocity(k);
   
    % mean values loaded from initialization simOut
    load(strcat(pwd,filesep,'data',filesep,simInput.paramFSI.wingParams.airfoil,filesep,'ROM',...
        filesep,'V',num2str(V0),filesep,'simOut.mat'), 'simOut');
    
    cL_bar_vec(k) = simOut.cL(1);
    cD_bar_vec(k) = simOut.cD(1);
    cRoll_bar_vec(k) = simOut.cRoll(1);
    cYaw_bar_vec(k) = simOut.cYaw(1);
    cPitch_bar_vec(k) = simOut.cPitch(1);
    bendingModeAmplitude_bar_vec(k) = simOut.Snapshot(610,1);
    
    Y_0(1,:,k) = simOut.cL(2:end) - cL_bar_vec(k);
    Y_0(2,:,k) = simOut.cD(2:end) - cD_bar_vec(k);
    Y_0(3,:,k) = simOut.cRoll(2:end) - cRoll_bar_vec(k);
    Y_0(4,:,k) = simOut.cPitch(2:end) - cPitch_bar_vec(k);
    Y_0(5,:,k) = simOut.cYaw(2:end) - cYaw_bar_vec(k);
    Y_0(6,:,k) = simOut.bendingModeAmplitude(2:end) - bendingModeAmplitude_bar_vec(k);
    
    % Subtracting mean or not
    X_0(:,:,k) = X_0(:,:,k) - Xmean(:,k); 
    X_1(:,:,k) = X_1(:,:,k) - Xmean(:,k);
    
    % GRAMIANS
    % 'gram' -> gramian calulated using linearized dynamics (linearized dynamics: IODMD) and solving Lyapunov equation
    % 'empGlin' -> empirical gramian calulated using linearized dynamics
    % (linearized dynamics: IODMD), using Lall 2002
    % 'empGnlin' -> empirical gramians calculated using nonlinear FSI 
    type_string_gram = 'gram';
    loadGramians = false; % load gramians from files or generate them
    [W_o, W_c] = generate_gramian(Snapshot(:,:,k), Snapshot_input(:,:,k), V0, type_string_gram, loadGramians, rankMax, simInput);
        
    L_c_mat{k} = cholcov(W_c).';
    L_o_mat{k} = cholcov(W_o).';
    if(size(L_o_mat{k},2) < smallest_size) % for the u-loop
        smallest_size = size(L_o_mat{k},2);
    end
    if(size(L_c_mat{k},2) < smallest_size) % for the u-loop
        smallest_size = size(L_c_mat{k},2);
    end
end



for r = 1:rankMax
    [F_list{r}, G_list{r}, H_list{r}, D_list{r}, L_list{r}, E_list{r}, W_list{r},V_list{r}, ~] = BMD_hf(r, X_0, X_1, U_0, Y_0, U_1, L_c_mat, L_o_mat); % rel_fit_error_vec_mat(r)
end


BMD_Matrices.F_list = F_list;
BMD_Matrices.G_list = G_list;
BMD_Matrices.H_list = H_list;
BMD_Matrices.D_list = D_list;
BMD_Matrices.L_list = L_list;
BMD_Matrices.E_list = E_list;
BMD_Matrices.W_list = W_list;
BMD_Matrices.V_list = V_list;

BMD_Matrices.Xmean = Xmean;
BMD_Matrices.rank = rankMax;
BMD_Matrices.type_string = type_string;
BMD_Matrices.velocity_vec = grid_velocity;

BMD_Matrices.cL_bar_vec = cL_bar_vec;
BMD_Matrices.cD_bar_vec = cD_bar_vec;
BMD_Matrices.cRoll_bar_vec = cRoll_bar_vec;
BMD_Matrices.cYaw_bar_vec = cYaw_bar_vec;
BMD_Matrices.cPitch_bar_vec = cPitch_bar_vec;
BMD_Matrices.bendingModeAmplitude_bar_vec = bendingModeAmplitude_bar_vec;

