function [t1a, t2a, t1b, t2b] = circleCircleTangent(X1, X2, r1, r2)
% [t1a t2a t1b t2b] = circleCircleTangent(X1, X2, R1, R2)
% Calculates external tangents to circles
%
% #### WARNING: FAILS IF THE CIRCLES ARE LOCATED AT THE SAME X or Y COORDINATE
%
%   Usage Example:
%
%   circleCircleTangent([1 3], [5 3], 3, 1); 
%
%   Kaushik Kandakatla <kkaushik@ymail.com>
%   Version 1.00
%   July, 2011

p1=X1(1);
q1=X1(2);
p2=X2(1);
q2=X2(2);

if p1 == p2 %###### NASTY HACK
	p1 = p1 + eps;
end
if q1 == q2
	p2 = p2 + eps;
end


d2 = (p2-p1)^2+(q2-q1)^2;
r = sqrt(d2-(r2-r1)^2);
s = ((q2-q1)*r+(p2-p1)*(r2-r1))/d2;
c = ((p2-p1)*r-(q2-q1)*(r2-r1))/d2;
x1 = p1-r1*s;
y1 = q1+r1*c;
x2 = p2-r2*s;
y2 = q2+r2*c;
m=(q2-q1)/(p2-p1);
xc1=(p1*m^2-m*q1+m*y1+x1)/(m^2+1);
yc1=(y1*m^2+(x1-p1)*m+q1)/(m^2+1);
xp1=2*xc1-x1;
yp1=2*yc1-y1;
xc2=(p1*m^2-m*q1+m*y2+x2)/(m^2+1);
yc2=(y2*m^2+(x2-p1)*m+q1)/(m^2+1);
xp2=2*xc2-x2;
yp2=2*yc2-y2;

t1a = [x1; y1];
t2a = [x2; y2];
t1b = [xp1; yp1];
t2b = [xp2; yp2];