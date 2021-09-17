function q = vortex_semiinf(r, v, l)
% Q = VORTEX_SEMIINF(R, V, L) calculates the speed induced by a semi-infinite vortex
%
% 	computes the velocity induced at the points in the columns of R by
% 	the unit-strength semi-infinite rectilinear vortex segments beginning
% 	from the points in the columns of V and extending indefinitely
% 	in the directions in the columns of L.

	%q = vortex_segment (r, v, v + 1e9*l);
	lhat = normalize(l);
	h = crossC(lhat, r - v);
	normh2 = dotC(h, h);
	%q = repmat (dot (lhat / 4 / pi, ...
	%	normalize(r - v))  + 1 ...
	%	./ normh2, [3, 1, 1] ) .* h;
	%q = repmat (dot (lhat / 4 / pi, ...
	%	normalize(r - v))  + 1 ...
	%	./ normh2, [3, 1] ) .* h;
	q = [	dotC(lhat / 4 / pi, normalize(r - v)) + 1 ./ normh2; ...
			dotC(lhat / 4 / pi, normalize(r - v)) + 1 ./ normh2; ...
			dotC(lhat / 4 / pi, normalize(r - v)) + 1 ./ normh2] .* h;
	%q(repmat (normh2 < eps, [3, 1, 1])) = 0; % r collinear
	%q(repmat (normh2 < eps, [3, 1])) = 0; % r collinear
	q([normh2 < eps; normh2 < eps; normh2 < eps]) = 0; % r collinear
end