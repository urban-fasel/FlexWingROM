% Generate the ROM and then perform simulation

function ROMsim(V0,rankROM,type_string,dT,iTest,simInput,simOutFULL)

n_g = length(V0); % number of velocities -> number of ROMs

% load snapshots matrices for each velocity V0
for k = 1:n_g
    Snapshot_temp = load(strcat(pwd,filesep,'data',filesep,simInput.paramFSI.wingParams.airfoil,filesep,'ROM',filesep,...
        'V',num2str(V0(k)),filesep,'Snapshot',filesep,'Snapshot.csv'));
    Snapshot_temp = transpose(Snapshot_temp);
    Snapshot_input_temp = load(strcat(pwd,filesep,'data',filesep,simInput.paramFSI.wingParams.airfoil,filesep,'ROM',filesep,...
        'V',num2str(V0(k)),filesep,'Snapshot',filesep,'Snapshot_input.csv'));
    Snapshot_input_temp = transpose(Snapshot_input_temp);
    
    Snapshot(:,:,k) = Snapshot_temp; % cropping, since extended sometimes
    Snapshot_input(:,:,k) = Snapshot_input_temp;
end


% generate ROM
switch type_string
    case 'aDMDc' % aDMDc model
        disp('Generate aDMDc models');
        ROMs = Obtain_aDMDc_model(Snapshot, Snapshot_input, V0, type_string, rankROM);
        rankMax = rankROM.rankMax;
    case 'aIODMD' % IODMD model
        disp('Generate aIODMD models');
        ROMs = Obtain_aIODMD_model(Snapshot, Snapshot_input,V0, type_string, rankROM.rankMax, simInput);
        rankMax = rankROM.rankMax;
    case 'BMD' % BMD model   
        disp('Generate BMD models');
        ROMs = Obtain_BMD_model(Snapshot, Snapshot_input, V0, type_string, rankROM.rankMaxBMD, simInput);
        rankMax = rankROM.rankMaxBMD;
end


% simulate ROM
simOutROM = cell(1, rankMax);

for r = 1:rankMax
    switch type_string
        case 'aDMDc' % aDMDc model
            disp(['Run aDMDc model with rank = ',num2str(r)]);
            ROMs.Atil_mat=ROMs.Atil_list{r};
            ROMs.Btil_mat=ROMs.Btil_list{r};
            ROMs.Ftil_mat=ROMs.Ftil_list{r};
            ROMs.Uhat_mat=ROMs.Uhat_list{r};

        case 'aIODMD' % IODMD model
            disp(['Run aIODMD model with rank = ',num2str(r)]);
            ROMs.F_mat=ROMs.F_list{r};
            ROMs.G_mat=ROMs.G_list{r};
            ROMs.H_mat=ROMs.H_list{r};
            ROMs.D_mat=ROMs.D_list{r};
            ROMs.L_mat=ROMs.L_list{r};
            ROMs.E_mat=ROMs.E_list{r};
            ROMs.Q=ROMs.Q_list{r};

        case 'BMD' % BMD model   
            disp(['Run BMD model with rank = ',num2str(r)]);
            ROMs.F_mat=ROMs.F_list{r};
            ROMs.G_mat=ROMs.G_list{r};
            ROMs.H_mat=ROMs.H_list{r};
            ROMs.D_mat=ROMs.D_list{r};
            ROMs.L_mat=ROMs.L_list{r};
            ROMs.E_mat=ROMs.E_list{r};
            ROMs.W_mat=ROMs.W_list{r};
            ROMs.V=ROMs.V_list{r};
            
        end

    simOutROM_temp = runFSI_ROM(V0,ROMs, dT, iTest, simInput, type_string); 
    simOutROM{r} = simOutROM_temp;
end

% save for comparison and plotting
save_ushape(simOutFULL, simOutROM, type_string, V0, simInput);
save_trajectory(simOutROM, type_string, V0, simInput)


