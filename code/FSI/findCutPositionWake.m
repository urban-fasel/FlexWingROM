function cutWakeAt = findCutPositionWake(paramWing,viscPre,Dt)
% U = 20;       % ESTIMATED MIN SPEED IN THE SIMULATION
% warning('fix min speed wake')
% U = 50;       % ESTIMATED MIN SPEED IN THE SIMULATION
% S = paramWing.S_mw;
% Cl = max(viscPre.maxCl);  % MAX Cl
% % Steps = 25:5:200;
% Steps = 25:5:400;
% r1 = sqrt((paramWing.b/2)^2+(Steps*U*Dt).^2);
% V_induced = ( S * Cl ) ./ ( 8 * pi * Dt * Steps .*r1 );
% requiredSteps = Steps(find(V_induced < 0.001*U,1,'first')); % We look for the first number of steps that produce an induced angle less than 0.001 (three orders of magnitude less than when the wake is cut immediately)
% if requiredSteps == Steps(end)
%     disp('Wake not long enough')
% end
% cutWakeAt = requiredSteps;
warning(' !!!!!!!!!!!!! fixed wake length !!!!!!!!!!!!!!!')
cutWakeAt = 100;
end