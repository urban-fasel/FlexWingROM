function [t, u] = segment_intersect_2d(p, dp, q, dq)
% Finds the intersections between segment (p - p+dp) and segment (q - q+dq) such that p+t*dp = q+u*dq

% (p+t*dp) = (q+u*dq)
% (p+t*dp) × dq = (q+u*dq) × dq
% t*(dp × dq) = (q - p) × dq
% OR
% (p+t*dp) × dp = (q+u*dq) × dp
% (p - q) × dp = u*(dq × dp)
% --> t = (q - p) × dq / (dp × dq)
% --> u = (p - q) × dp / (dq × dp)

t = twodcross((q - p), dq) / twodcross(dp, dq);
u = twodcross((p - q), dp) / twodcross(dq, dp);

end