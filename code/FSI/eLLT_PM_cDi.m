function [cDi_glob,cYaw,FN] = eLLT_PM_cDi(y_segLimit, c_ind, V, cLdistr_val_interp, K_alpha_ind_ELLT,q,Sw,b)
% y_segLimit, c_ind, V, cLdistr_val_interp, K_alpha_ind_ELLT

n_seg = length(y_segLimit)-1;
dy_seg = (y_segLimit(2:end) - y_segLimit(1:end-1)); % span of each segment
y_centerpoint = (y_segLimit(2:end) + y_segLimit(1:end-1))/2;

WingArea = sum(c_ind*dy_seg);

%% First step: get the cP on all sections, then Gamma
Gamma = zeros(n_seg, 1);
for i = 1:n_seg/2
	Gamma(i) = 0.5 * V * c_ind * cLdistr_val_interp(i);
end
for i = (n_seg/2+1):n_seg
	Gamma(i) = 0.5 * V * c_ind * cLdistr_val_interp(i);
end

%% Second step: create AIC matrix
K_alpha_ind = K_alpha_ind_ELLT/(4*pi*V);

%% Drag
alpha_ind = -K_alpha_ind * Gamma; % -K·Gamma returns the induced angle of attack
cDi_sect = 2*(Gamma .* (alpha_ind) ./ c_ind)./V; 
cDi_glob = sum(cDi_sect .* c_ind .* dy_seg(:))/WingArea;
cYaw = -dot(cDi_sect .* c_ind .* dy_seg(:),y_centerpoint)/WingArea/b;
FN = cYaw*q*Sw*b;