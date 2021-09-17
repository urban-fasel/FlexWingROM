%% Script, which calculates the empirical gramians MISO for different speeds
% Note: here the data computed in the other data generating routines of
% this folder are used. Therefore, those simulations must be run first (for
% MISO)

function generate_empirical_gramian(V, nInputs, InpAmpO, InpAmpC, maxIterationO, maxIterationC, simInput)

n_x = 618;
r = 1;
% s = length(c_aero_vec);
maxIteration = maxIterationO;%300;

% Loading data with trim point
% load(['data/simOutFULL_0_c_m_8deg_V',num2str(V),'_maxIteration_300.mat'],'simOutFULL_0'); % loads trim point
% load(['data/simOut0_ego_V',num2str(V),'.mat'],'simOut')
load(strcat('data',filesep,simInput.paramFSI.wingParams.airfoil,filesep,'ROM',...
	filesep,'V',num2str(V),filesep,'simOut0_ego_V',num2str(V),'.mat'),'simOut'); 

simOutFULL_0 = simOut;

x_bar = simOutFULL_0.Snapshot(:,30); 

% Observability Gramian loop
X = zeros(n_x, n_x);
Psi_temp = zeros(n_x, n_x, maxIteration);
Psi = zeros(n_x, n_x);

% maybe need to loop over several input amplitudes
c_aero = InpAmpO(1);%c_aero_vec(m);
c_struct = InpAmpO(2);%c_struct_vec(m);
absolute_amplitudes_aero = c_aero*x_bar - x_bar;
absolute_amplitudes_struct = c_struct*x_bar - x_bar;
% load(['data/simOutEG_two_relative_ImpAmp_2_05_nStates_618_V',num2str(V),'_maxIteration_',num2str(maxIteration),'.mat'],'simOutEG');
% load(['data/simOut_ego_V',num2str(V),'.mat'],'simOutEG')
load(strcat('data',filesep,simInput.paramFSI.wingParams.airfoil,filesep,'ROM',...
	filesep,'V',num2str(V),filesep,'simOut_ego_V',num2str(V),'.mat'),'simOutEG'); 

for i = 1:n_x
    if(i < 609)
        c_i = absolute_amplitudes_aero(i);
    else
        c_i = absolute_amplitudes_struct(i);
    end
    dz_im = simOutEG(i).Snapshot(610,:);
    idx_im = find(dz_im == 0);
    if(~isempty(idx_im))
        dz_im(1:idx_im(1)-1) = dz_im(1:idx_im(1)-1) - x_bar(610); % only subtract values until simulation was run & ignore the zero initialization
    end
    dz_im_bar = mean(dz_im, 2);
    temp_i = (dz_im - dz_im_bar)/c_i;

    for j = i:n_x % symmetric & single output start at i:n
        if(j < 609)
            c_j = absolute_amplitudes_aero(j);
        else
            c_j = absolute_amplitudes_struct(j);
        end
        dz_jm = simOutEG(j).Snapshot(610,:);
        idx_jm = find(dz_jm == 0);
        if(~isempty(idx_jm))
            dz_jm(1:idx_jm(1)-1) = dz_jm(1:idx_jm(1)-1) - x_bar(610);
        end
        dz_jm_bar = mean(dz_jm, 2);
        temp_j = (dz_jm - dz_jm_bar)/c_j;

        for k = 1:maxIteration
            temp = temp_i(:,k)'*temp_j(:,k); % Psi_ij(t) at time k
            Psi_temp(i,j,k) = temp; % third dimension element-wise time
            Psi_temp(j,i,k) = temp; % symmetric matrix, assuming no complex vectors
        end
    end
end % end of building Psi

% Integration
Psi(:,:) = sum(Psi_temp,3);
integral = Psi;

% X = X + 1/(r*s)*integral;
X = X + 1/(r)*integral;


%% Controllability Gramian
% Preparation
amplitudes = InpAmpC;
s = length(amplitudes);
r = 1;
maxIteration = maxIterationC; 

% Controllability Gramian loop
Y = zeros(n_x, n_x);
for i = 1:nInputs
    for m = 1:s
        c_m = amplitudes(m);
%         load(['data/simOut',num2str(i),'_egc_V',num2str(V),'.mat'],'simOut')
        load(strcat('data',filesep,simInput.paramFSI.wingParams.airfoil,filesep,'ROM',...
            filesep,'V',num2str(V),filesep,'simOut',num2str(i),'_egc_V',num2str(V),'.mat'),'simOut'); 

        simOutFULL = simOut;

        dx_im = simOutFULL.Snapshot(:,1:maxIteration) - x_bar;
        dx_im_bar = mean(dx_im, 2);
        temp = dx_im - dx_im_bar;

        % integral
        integral = (temp*temp');
        Y = Y + 1/(r*s*c_m^2)*integral;
    end
end

%% Saving
W_o = X;
W_c = Y;

%     save(['testing/EmpiricalG/two_relative_empirical_gramian_V',num2str(V),'.mat'],'W_o','W_c');
% save(['data/empirical_gramians_V',num2str(V),'.mat'],'W_o','W_c');
save(strcat('data',filesep,simInput.paramFSI.wingParams.airfoil,filesep,'ROM',...
            filesep,'V',num2str(V),filesep,'empirical_gramians_V',num2str(V),'.mat'),'W_o','W_c'); 


