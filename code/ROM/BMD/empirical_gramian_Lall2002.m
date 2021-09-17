%% implementation of the method to compute empirical Gramians from
% "Lall - 2002 - A subspace approach to balanced truncation for model reduction of nonlinear control systems - IJRNC"
%% Input
% G_lin: linear discrete-time system 
% n_x
% n_u
% M: Vector of the amplitude
% dt: discretization of the integral
% T_c: simulation end for contrallability gramian
% T_o: simulation end for observability gramian

%% Output
% Y is the Empirical Controllability Gramian
% X is the Empirical Observability Gramian
function [Y, X] = empirical_gramian_Lall2002(G_lin, n_x, n_u, n_y, M, dt, T_c, T_o)
%% Empirical Controllability Gramian
r = 2;
s = length(M);
p = n_u;

F = zeros(p, p, r);
F(:,:,1) = eye(p);
F(:,:,2) = -eye(p);
Eps = eye(p);

Y = zeros(n_x, n_x);
t_c = 0:dt:T_c;
N = length(t_c);
dx0 = zeros(n_x,1);
for l = 1:r
    T_l = F(:,:,l);
    
    for m = 1:s
        c_m = M(1,m);
        
        for i = 1:p
            e_i = Eps(:,i);

            du = zeros(n_u,N);
            du(:,1) = c_m*T_l*e_i; % divide by dt to get dirac-area of 1
            
            [~, ~, dx_ilm] = lsim(G_lin, du.', t_c, dx0);
            dx_ilm = dx_ilm.'; % n_x x t
            
%             assert(norm(dx_ilm(:,end),2) < 1e-6);
            if ~(norm(dx_ilm(:,end),2) < 1e-6)
                warning(['!!!! empirical gramian tolerance ',num2str(norm(dx_ilm(:,end),2))])
            end
            % [rel, ~] = max(vecnorm(dx_ilm, 2, 2)); % maximum of 2-norm at some time
            % assert(norm(dx_ilm(:,end),2)/rel < 1e-3 && norm(dx_ilm(:,end),2) < 1e-3);
            
            dx_ilm_bar = mean(dx_ilm, 2);
            
            temp = dx_ilm - dx_ilm_bar; 
            
            integral = (temp*temp');
            
            Y = Y + 1/(r*s*c_m^2)*integral;
        end
    end
end


%% Empirical Observability Gramian
F = zeros(n_x, n_x, r);
F(:,:,1) = eye(n_x);
F(:,:,2) = -eye(n_x);
Eps = eye(n_x);

X = zeros(n_x, n_x);
t_o = 0:dt:T_o;
N = length(t_o);
du = zeros(n_u, N);
Psi_temp = zeros(n_x,n_x,N); % thrid dimension time (element-wise)
Psi = zeros(n_x,n_x);
for l = 1:r
    T_l = F(:,:,l);
    
    for m = 1:s
        c_m = M(:,m);
        
        % Building Psi
        dz_ilm_temp = zeros(N,n_y,n_x); % cL, cD, cRoll, cPitch, cYaw, bendingModeAmplitude
        for i = 1:n_x
            e_i = Eps(:,i);
            dx0_i = c_m*T_l*e_i;
            [dz_ilm_temp(:,:,i),~, ~] = lsim(G_lin, du.', t_o, dx0_i);
        end
        
        for i = 1:n_x
            dz_ilm = dz_ilm_temp(:,:,i);
            dz_ilm = dz_ilm.'; % 1 x N
            
%             assert(norm(dz_ilm(:,end),2) < 1e-6);
            if ~(norm(dz_ilm(:,end),2) < 1e-6)
                warning(['!!!! empirical gramian tolerance ',num2str(norm(dz_ilm(:,end),2))])
            end
            % [rel, ~] = max(vecnorm(dz_ilm, 2, 2)); % maximum of 2-norm at some time   
            % assert(norm(dz_ilm(:,end),2)/rel < 1e-3 && norm(dz_ilm(:,end),2) < 1e-3);
            
            dz_ilm_bar = mean(dz_ilm, 2);
            
            temp_i = dz_ilm - dz_ilm_bar; 
            
            for j = i:n_x
                if(j == i)
                    temp_j = temp_i;
                else
                    dz_jlm = dz_ilm_temp(:,:,j);
                    dz_jlm = dz_jlm.';
                    
%                     assert(norm(dz_jlm(:,end),2) < 1e-6);
                    if ~(norm(dz_jlm(:,end),2) < 1e-6)
                        warning(['!!!! empirical gramian tolerance ',num2str(norm(dz_jlm(:,end),2))])
                    end
                    % [rel, ~] = max(vecnorm(dz_jlm, 2, 2)); % maximum of 2-norm at some time
                    % assert(norm(dz_jlm(:,end),2)/rel < 1e-3 && norm(dz_jlm(:,end),2) < 1e-3);
                    
                    dz_jlm_bar = mean(dz_jlm, 2);

                    temp_j = dz_jlm - dz_jlm_bar; 
                end
                
                for k = 1:N
                    temp = temp_i(:,k)'*temp_j(:,k); % Psi_ij(t) at time k
                    Psi_temp(i,j,k) = temp; % third dimension element-wise time
                    Psi_temp(j,i,k) = temp; % symmetric matrix, assuming no complex vectors
                end
            end
        end % end of building Psi
        
        % Integration
        Psi(:,:) = sum(Psi_temp, 3);
        integral = T_l*Psi*T_l'; % Linearity of integral
        
        X = X + 1/(r*s*c_m^2)*integral;
    end
end

end