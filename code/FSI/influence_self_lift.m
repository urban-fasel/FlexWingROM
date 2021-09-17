% Copyright (C) 2008  Daniel Filkovic

% This file is part of APAME-ADS.

% APAME-ADS is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% APAME-ADS is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with APAME-ADS.  If not, see <http://www.gnu.org/licenses/>.

% file influence_self_lift.m

function [atotal,b]=influence_self_lift_new(FF,ep,S,M,N,N1,O,l1,l2,l3,t1,t2,t3,n1,n2,n3,cx,cy,cz,x1,y1,x2,y2,x3,y3,x4,y4,d1,d2,d3,d4)

% calculation of influence coefficients
a = zeros(M*N,M*N);
b = zeros(M*N,M*N);
I = 0;
for k = 1:N
    for l = 1:M
        for i = 1:N
            for j = 1:M
                I = I+1;
                cxijcxkl = cx(i,j)-cx(k,l);
                cyijcykl = cy(i,j)-cy(k,l);
                czijczkl = cz(i,j)-cz(k,l);
                cpx = cxijcxkl*l1(k,l) + cyijcykl*l2(k,l) + czijczkl*l3(k,l);
                cpy = cxijcxkl*t1(k,l) + cyijcykl*t2(k,l) + czijczkl*t3(k,l);
                cpz = cxijcxkl*n1(k,l) + cyijcykl*n2(k,l) + czijczkl*n3(k,l);
                % if distance of panel from influenced point is greater
                % then product of longer diagonal and "far field" coefficient
                if realsqrt(cpx^2+cpy^2+cpz^2) > FF(k,l)
                    rad = cpx^2 + cpy^2 + cpz^2;
                    a(I) = -1/(4*pi)*S(k,l)*cpz*rad^(-1.5);
                    b(I) = -1/(4*pi)*S(k,l)/realsqrt(rad);
                else
                    cpx1 = cpx - x1(k,l);
                    cpx2 = cpx - x2(k,l);
                    cpx3 = cpx - x3(k,l);
                    cpx4 = cpx - x4(k,l);
                    cpy1 = cpy - y1(k,l);
                    cpy2 = cpy - y2(k,l);
                    cpy3 = cpy - y3(k,l);
                    cpy4 = cpy - y4(k,l);
                    e1 = cpx1^2+cpz^2;
                    e2 = cpx2^2+cpz^2;
                    e3 = cpx3^2+cpz^2;
                    e4 = cpx4^2+cpz^2;
                    r1 = realsqrt(e1 + cpy1^2);
                    r2 = realsqrt(e2 + cpy2^2);
                    r3 = realsqrt(e3 + cpy3^2);
                    r4 = realsqrt(e4 + cpy4^2);
                    h1 = cpx1*cpy1;
                    h2 = cpx2*cpy2;
                    h3 = cpx3*cpy3;
                    h4 = cpx4*cpy4;
                    x21 = x2(k,l)-x1(k,l);
                    x32 = x3(k,l)-x2(k,l);
                    x43 = x4(k,l)-x3(k,l);
                    x14 = x1(k,l)-x4(k,l);
                    y21 = y2(k,l)-y1(k,l);
                    y32 = y3(k,l)-y2(k,l);
                    y43 = y4(k,l)-y3(k,l);
                    y14 = y1(k,l)-y4(k,l);
                    if d1(k,l) < ep
                        a1 = 0;
                        b1 = 0;
                    else
                        F = y21*e1 - x21*h1;
                        G = y21*e2 - x21*h2;
                        a1 = atan2(cpz*x21*(F*r2-G*r1), cpz^2*x21^2*r1*r2+F*G);
                        b1 = (cpx1*y21-cpy1*x21)/d1(k,l)*reallog((r1+r2+d1(k,l))/(r1+r2-d1(k,l)));
%                         a1 = atan((y21*(cpx1^2+cpz^2)-x21*cpx1*cpy1)/(r1*cpz*x21)) - atan((y21*(cpx2^2+cpz^2)-x21*cpx2*cpy2)/(r2*cpz*x21));
                    end
                    if d2(k,l) < ep
                        a2 = 0;
                        b2 = 0;
                    else
                        F = y32*e2 - x32*h2;
                        G = y32*e3 - x32*h3;
                        a2 = atan2(cpz*x32*(F*r3-G*r2), cpz^2*x32^2*r2*r3+F*G);
                        b2 = (cpx2*y32-cpy2*x32)/d2(k,l)*reallog((r2+r3+d2(k,l))/(r2+r3-d2(k,l)));
%                         a2 = atan((y32*(cpx2^2+cpz^2)-x32*cpx2*cpy2)/(r2*cpz*x32)) - atan((y32*(cpx3^2+cpz^2)-x32*cpx3*cpy3)/(r3*cpz*x32));
                    end
                    if d3(k,l) < ep
                        a3 = 0;
                        b3 = 0;
                    else
                        F = y43*e3 - x43*h3;
                        G = y43*e4 - x43*h4;
                        a3 = atan2(cpz*x43*(F*r4-G*r3), cpz^2*x43^2*r3*r4+F*G);
                        b3 = (cpx3*y43-cpy3*x43)/d3(k,l)*reallog((r3+r4+d3(k,l))/(r3+r4-d3(k,l)));
%                         a3 = atan((y43*(cpx3^2+cpz^2)-x43*cpx3*cpy3)/(r3*cpz*x43)) - atan((y43*(cpx4^2+cpz^2)-x43*cpx4*cpy4)/(r4*cpz*x43));
                    end
                    if d4(k,l) < ep
                        a4 = 0;
                        b4 = 0;
                    else
                        F = y14*e4 - x14*h4;
                        G = y14*e1 - x14*h1;
                        a4 = atan2(cpz*x14*(F*r1-G*r4), cpz^2*x14^2*r4*r1+F*G);
                        b4 = (cpx4*y14-cpy4*x14)/d4(k,l)*reallog((r4+r1+d4(k,l))/(r4+r1-d4(k,l)));
%                         a4 = atan((y14*(cpx4^2+cpz^2)-x14*cpx4*cpy4)/(r4*cpz*x14)) - atan((y14*(cpx1^2+cpz^2)-x14*cpx1*cpy1)/(r1*cpz*x14));
                    end
                    dipol = -(a1+a2+a3+a4)*1/(4*pi);
                    if abs(cpz) < ep
                        dipol = 0;
                    end
                    if i == k && j == l
                        dipol = 0.5;
                    end
                    a(I) = dipol;                                                     % Doublet influence coefficient
                    b(I) = -(b1+b2+b3+b4)*1/(4*pi)-cpz*dipol;                % Source influence coefficient
                end
            end
        end
    end
end

% calculation of trail wake influence coefficients
I = 0;
atrag = zeros(M*N,M*N);
for k = 1:N
    for l = 1:M
        for i = 1:N
            for j = 1:M
                I = I+1;
                if (k == N)  || (k == 1)                                   % We are here at the trailing edge
                    cxijcxN1l = cx(i,j)-cx(N1,l);
                    cyijcyN1l = cy(i,j)-cy(N1,l);
                    czijczN1l = cz(i,j)-cz(N1,l);
                    cpx = cxijcxN1l*l1(N1,l) + cyijcyN1l*l2(N1,l) + czijczN1l*l3(N1,l);
                    cpy = cxijcxN1l*t1(N1,l) + cyijcyN1l*t2(N1,l) + czijczN1l*t3(N1,l);
                    cpz = cxijcxN1l*n1(N1,l) + cyijcyN1l*n2(N1,l) + czijczN1l*n3(N1,l);
                    if realsqrt(cpx^2+cpy^2+cpz^2) > FF(N1,l)
                        rad = cpx^2 + cpy^2 + cpz^2;
                        atrag(I) = -1/(4*pi)*S(N1,l)*cpz*rad^(-1.5);
                    else
                        cpx1 = cpx - x1(N1,l);
                        cpx2 = cpx - x2(N1,l);
                        cpx3 = cpx - x3(N1,l);
                        cpx4 = cpx - x4(N1,l);
                        cpy1 = cpy - y1(N1,l);
                        cpy2 = cpy - y2(N1,l);
                        cpy3 = cpy - y3(N1,l);
                        cpy4 = cpy - y4(N1,l);
                        e1 = cpx1^2+cpz^2;
                        e2 = cpx2^2+cpz^2;
                        e3 = cpx3^2+cpz^2;
                        e4 = cpx4^2+cpz^2;
                        r1 = realsqrt(e1 + cpy1^2);
                        r2 = realsqrt(e2 + cpy2^2);
                        r3 = realsqrt(e3 + cpy3^2);
                        r4 = realsqrt(e4 + cpy4^2);
                        h1 = cpx1*cpy1;
                        h2 = cpx2*cpy2;
                        h3 = cpx3*cpy3;
                        h4 = cpx4*cpy4;
                        x21 = x2(N1,l)-x1(N1,l);
                        x32 = x3(N1,l)-x2(N1,l);
                        x43 = x4(N1,l)-x3(N1,l);
                        x14 = x1(N1,l)-x4(N1,l);
                        y21 = y2(N1,l)-y1(N1,l);
                        y32 = y3(N1,l)-y2(N1,l);
                        y43 = y4(N1,l)-y3(N1,l);
                        y14 = y1(N1,l)-y4(N1,l);
                        if d1(N1,l) < ep
                            a1 = 0;
                        else
                            F = y21*e1 - x21*h1;
                            G = y21*e2 - x21*h2;
                            a1 = atan2(cpz*x21*(F*r2-G*r1), cpz^2*x21^2*r1*r2+F*G);
%                             a1 = atan((y21*(cpx1^2+cpz^2)-x21*cpx1*cpy1)/(r1*cpz*x21)) - atan((y21*(cpx2^2+cpz^2)-x21*cpx2*cpy2)/(r2*cpz*x21));
                        end
                        if d2(N1,l) < ep
                            a2 = 0;
                        else
                            F = y32*e2 - x32*h2;
                            G = y32*e3 - x32*h3;
                            a2 = atan2(cpz*x32*(F*r3-G*r2), cpz^2*x32^2*r2*r3+F*G);
%                             a2 = atan((y32*(cpx2^2+cpz^2)-x32*cpx2*cpy2)/(r2*cpz*x32)) - atan((y32*(cpx3^2+cpz^2)-x32*cpx3*cpy3)/(r3*cpz*x32));
                        end
                        if d3(N1,l) < ep
                            a3 = 0;
                        else
                            F = y43*e3 - x43*h3;
                            G = y43*e4 - x43*h4;
                            a3 = atan2(cpz*x43*(F*r4-G*r3), cpz^2*x43^2*r3*r4+F*G);
%                             a3 = atan((y43*(cpx3^2+cpz^2)-x43*cpx3*cpy3)/(r3*cpz*x43)) - atan((y43*(cpx4^2+cpz^2)-x43*cpx4*cpy4)/(r4*cpz*x43));
                        end
                        if d4(N1,l) < ep
                            a4 = 0;
                        else
                            F = y14*e4 - x14*h4;
                            G = y14*e1 - x14*h1;
                            a4 = atan2(cpz*x14*(F*r1-G*r4), cpz^2*x14^2*r4*r1+F*G);
%                             a4 = atan((y14*(cpx4^2+cpz^2)-x14*cpx4*cpy4)/(r4*cpz*x14)) - atan((y14*(cpx1^2+cpz^2)-x14*cpx1*cpy1)/(r1*cpz*x14));
                        end
                        atrag(I) = -(a1+a2+a3+a4)*1/(4*pi);
                        if abs(cpz) < ep
                            atrag(I) = 0;
                        end
                    end
                    if k == 1   % lower trailing edge
                        atrag(I) = -atrag(I);
                    end
                else
                    atrag(I) = 0;
                end
            end
        end
    end
end

atotal = a+atrag;