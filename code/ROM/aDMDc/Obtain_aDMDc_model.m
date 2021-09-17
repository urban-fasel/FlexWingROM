% Function that generates the aDMDc model 

function [aDMDc_Matrices] = Obtain_aDMDc_model(Snapshot,Snapshot_input,grid_velocity,type_string,rankROM)

[n_x, Nt, n_g] = size(Snapshot);
[n_u,~,~] = size(Snapshot_input);
n_s = Nt -1;
n = n_x;
q = n_u;

rank_X_0 = inf;
Util_full_mat = zeros(n_x + 2*n_u, n_s, n_g);
Sigtil_full_mat = zeros(n_s, n_s, n_g);
Vtil_full_mat = zeros(n_s, n_s, n_g);
U_mat = zeros(n_x, n_s, n_g);
Sig_mat = zeros(n_s, n_s, n_g);
V_mat = zeros(n_s, n_s, n_g);
Xmean = zeros(n_x,n_g);
for k = 1:n_g
    X = Snapshot(:,1:end-1,k);
    Xp = Snapshot(:,2:end,k);
    
    Xmean(:,k) = mean(X,2);
    X = X - Xmean(:,k);
    Xp = Xp - Xmean(:,k);
    
    UPS = Snapshot_input(:,1:end-1,k);
    UPSp = Snapshot_input(:,2:end,k);
    XModified = [X;UPS;UPSp];
    [Util_full,Sigtil_full,Vtil_full] = svd(XModified,'econ');
    [U,Sig,V] = svd(Xp,'econ');

    Util_full_mat(:,:,k) = Util_full;
    Sigtil_full_mat(:,:,k) = Sigtil_full;
    Vtil_full_mat(:,:,k) = Vtil_full;
    U_mat(:,:,k) = U;
    Sig_mat(:,:,k) = Sig;
    V_mat(:,:,k) = V;
    
    rank_1 = rank(X);
    rank_2 = rank(Xp);
    if(rank_1 < rank_X_0)
        rank_X_0 = rank_1;
    end
    if(rank_2 < rank_X_0)
        rank_X_0 = rank_2;
    end
end

rank_X_0 = min(rankROM.rankMax, rank_X_0);


% Just build n_g models at each parameter
Atil_list = cell(1, rank_X_0);
Btil_list = cell(1, rank_X_0);
Ftil_list = cell(1, rank_X_0);
Uhat_list = cell(1, rank_X_0);
for r = 1:rank_X_0
    Atil_mat = zeros(r, r, n_g); % saving all models
    Btil_mat = zeros(r, n_u, n_g);
    Ftil_mat = zeros(r, n_u, n_g);
    Uhat_mat = zeros(n_x, r, n_g);

    % this first truncation, does not affect the rank of the final system,
    % however, also affects the accuracy of the model.
    rtil = max([r, rankROM.rankDMDcTruncation1]); 

	for k = 1:n_g
        Xp = Snapshot(:,2:end,k);
        Xp = Xp - Xmean(:,k);
        
        Util = Util_full_mat(:,1:rtil,k);
        Sigtil = Sigtil_full_mat(1:rtil,1:rtil,k);
        Vtil = Vtil_full_mat(:,1:rtil,k);
        
        Uhat = U_mat(:,1:r,k);
        Sighat = Sig_mat(1:r,1:r,k);
        Vhat = V_mat(:,1:r,k);
        
        U_1 = Util(1:n,:);
        U_2 = Util(n+1:n+q,:);
        U_3 = Util(n+q+1:n+q+q,:);
        
        Atil_mat(:,:,k) = Uhat'*(Xp)*Vtil*(Sigtil\U_1')*Uhat;
        Btil_mat(:,:,k) = Uhat'*(Xp)*Vtil*(Sigtil\U_2');
        Ftil_mat(:,:,k) = Uhat'*(Xp)*Vtil*(Sigtil\U_3');
        Uhat_mat(:,:,k) = Uhat;
    end
    Atil_list{r} = Atil_mat;
    Btil_list{r} = Btil_mat;
    Ftil_list{r} = Ftil_mat;
    Uhat_list{r} = Uhat_mat;
end

aDMDc_Matrices.Atil_list = Atil_list;
aDMDc_Matrices.Btil_list = Btil_list;
aDMDc_Matrices.Ftil_list = Ftil_list;
aDMDc_Matrices.Uhat_list = Uhat_list;
aDMDc_Matrices.Xmean = Xmean;

aDMDc_Matrices.rank = rank_X_0;
aDMDc_Matrices.type_string = type_string;
aDMDc_Matrices.velocity_vec = grid_velocity;
end