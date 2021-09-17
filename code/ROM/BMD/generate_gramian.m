function [W_o, W_c] = generate_gramian(Snapshot,Snapshot_input,V0,type_string_gram,loadGramians,rankMax,simInput)
    
if loadGramians
    load(strcat('data',filesep,paramFSI.wingParams.airfoil,filesep,'ROM',...
                filesep,'V',num2str(V0),filesep,'gramians_V',num2str(V0),'.mat'),'W_c','W_o'); 
else

    switch type_string_gram
        case 'gram' 
            % using aIODMD linear system
            type_string = 'aIODMD';
            rank_aIODMD = rankMax + 10; % need higher rank, otherwise not possible to generate BMD model with rankMax
            ROMsIODMD = Obtain_aIODMD_model(Snapshot, Snapshot_input,V0, type_string, rank_aIODMD, simInput);
            Q_mat = ROMsIODMD.Q_list{1,40};
            A_aIODMDc = Q_mat*ROMsIODMD.F_list{1,40}*Q_mat';
            B_aIODMDc = Q_mat*ROMsIODMD.G_list{1,40};
            C_aIODMDc = ROMsIODMD.H_list{1,40}*Q_mat'; 
            D_aIODMDc = ROMsIODMD.D_list{1,40};
            sys = ss(A_aIODMDc,B_aIODMDc,C_aIODMDc,D_aIODMDc,0.006);
            W_c = gram(sys,'c');
            W_o = gram(sys,'o');

        case 'empGlin' 
            % using aIODMD linear system
            type_string = 'aIODMD';
            rank_aIODMD = rankMax + 10; % need higher rank, otherwise not possible to generate BMD model with rankMax
            ROMsIODMD = Obtain_aIODMD_model(Snapshot, Snapshot_input,V0, type_string, rank_aIODMD, simInput);
            Q_mat = ROMsIODMD.Q_list{1,40};
            A_aIODMDc = Q_mat*ROMsIODMD.F_list{1,40}*Q_mat';
            B_aIODMDc = Q_mat*ROMsIODMD.G_list{1,40};
            C_aIODMDc = ROMsIODMD.H_list{1,40}*Q_mat'; 
            D_aIODMDc = ROMsIODMD.D_list{1,40};
            sys = ss(A_aIODMDc,B_aIODMDc,C_aIODMDc,D_aIODMDc,0.006);

            % generate empirical gramians: Lall 2002
            n_x = 618;
            n_u = 6;
            n_y = 6;
            M_gram = 0.001;
            dt_int = 0.006;
            T_c = 6.0;
            T_o = 6.0;
            [W_c, W_o] = empirical_gramian_Lall2002(sys, n_x, n_u, n_y, M_gram, dt_int, T_c, T_o);
            
       case 'empGnlin' 
            % generate data for empirical controllability gramian
            nIterationC = 300; % number of FSI iterations
            nInputs = 6; % number of inputs
            InpAmpC = 8/180*pi;
            empirical_gramian_controllability_loop(paramFSI,V0,dT,nIterationC,nInputs,InpAmpC)

            % generate data for empirical observability gramian
            nIterationO = 300; % number of FSI iterations
            InpAmpO = [1.2 1.5]; % different input amplitudes for aerodynamic and structural states: InpAmp(1) -> struct, InpAmp(2) -> aero
            empirical_gramian_observability_loop(paramFSI,V0,dT,nIterationO,InpAmpO)

            % generate empirical gramians
            [W_c, W_o] = generate_empirical_gramian(V0, nInputs, InpAmpO, InpAmpC, nIterationO, nIterationC, simInput);
    end
    
    save(strcat('data',filesep,simInput.paramFSI.wingParams.airfoil,filesep,'ROM',...
    	filesep,'V',num2str(V0),filesep,'gramians_V',num2str(V0),'.mat'),'W_c','W_o'); 

end
