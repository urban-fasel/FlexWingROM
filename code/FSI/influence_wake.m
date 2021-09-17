function [C]=influence_wake(FF,ep,S,M,N,Nwakes,l1,l2,l3,t1,t2,t3,n1,n2,n3,cx,cy,cz,x1,y1,x2,y2,x3,y3,x4,y4,d1,d2,d3,d4)

% calculation of influence coefficients
C = zeros(M*N,M*Nwakes);
I = 0;
for k = N+2:N+1+Nwakes
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
                    C(I) = -0.079577471545948*S(k,l)*cpz*rad^(-1.5);
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
                    else
                        F = y21*e1 - x21*h1;
                        G = y21*e2 - x21*h2;
                        a1 = atan2(cpz*x21*(F*r2-G*r1), cpz^2*x21^2*r1*r2+F*G);
                    end
                    if d2(k,l) < ep
                        a2 = 0;
                    else
                        F = y32*e2 - x32*h2;
                        G = y32*e3 - x32*h3;
                        a2 = atan2(cpz*x32*(F*r3-G*r2), cpz^2*x32^2*r2*r3+F*G);
                    end
                    if d3(k,l) < ep
                        a3 = 0;
                    else
                        F = y43*e3 - x43*h3;
                        G = y43*e4 - x43*h4;
                        a3 = atan2(cpz*x43*(F*r4-G*r3), cpz^2*x43^2*r3*r4+F*G);
                    end
                    if d4(k,l) < ep
                        a4 = 0;
                    else
                        F = y14*e4 - x14*h4;
                        G = y14*e1 - x14*h1;
                        a4 = atan2(cpz*x14*(F*r1-G*r4), cpz^2*x14^2*r4*r1+F*G);
                    end
                    dipol = -(a1+a2+a3+a4)*0.079577471545948;
                    if abs(cpz) < ep
                        dipol = 0;
                    end
                    if i == k && j == l
                        dipol = 0.5;
                    end
                    C(I) = dipol;                                                     % Doublet influence coefficient
                end
            end
        end
    end
end
end