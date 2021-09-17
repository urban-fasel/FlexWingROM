function [us,vs,ws,ud,vd,wd]=induced_velocity(FF,ep,S,M,N,Nwakes,l1,l2,l3,t1,t2,t3,n1,n2,n3,cx,cy,cz,x1,y1,x2,y2,x3,y3,x4,y4,d1,d2,d3,d4,x,y,z,TimeStep,Scope)
ep2 = 1e-4;
us = zeros((M+1)*(Nwakes+1),M*N);
vs=us;
ws=vs;
ud = zeros((M+1)*(Nwakes+1),M*(N+1+Nwakes));
vd = ud;
wd = vd;
I = 0;
%%% First the bounded sources + doublets
for k = 1:N
    for l = 1:M
        for i = N+2:N+2+Nwakes
            for j = 1:M+1
                I = I+1;
                xijcxkl = x(i,j)-cx(k,l);
                yijcykl = y(i,j)-cy(k,l);
                zijczkl = z(i,j)-cz(k,l);
                cpx = xijcxkl*l1(k,l) + yijcykl*l2(k,l) + zijczkl*l3(k,l);
                cpy = xijcxkl*t1(k,l) + yijcykl*t2(k,l) + zijczkl*t3(k,l);
                cpz = xijcxkl*n1(k,l) + yijcykl*n2(k,l) + zijczkl*n3(k,l);
                % if distance of panel from influenced point is greater
                % then product of longer diagonal and "far field" coefficient
                if realsqrt(cpx^2+cpy^2+cpz^2) > FF(k,l)
                    rad = cpx^2 + cpy^2 + cpz^2;
                    uus = 0.079577471545948*S(k,l)*cpx*rad^(-1.5);
                    vvs = 0.079577471545948*S(k,l)*cpy*rad^(-1.5);
                    wws = 0.079577471545948*S(k,l)*cpz*rad^(-1.5);
                    uud = 3*0.079577471545948*S(k,l)*cpx*cpz*rad^(-2.5);
                    vvd = 3*0.079577471545948*S(k,l)*cpy*cpz*rad^(-2.5);
                    wwd = -0.079577471545948*S(k,l)*(cpx^2+cpy^2-2*cpz^2)*rad^(-2.5);
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
                        K1 = 0;
                    else
                        F = y21*e1 - x21*h1;
                        G = y21*e2 - x21*h2;
                        a1 = atan2(cpz*x21*(F*r2-G*r1), cpz^2*x21^2*r1*r2+F*G);
                        r1r2x = cpy1*cpz-cpz*cpy2;
                        r1r2y = -cpx1*cpz+cpz*cpx2;
                        r1r2z = cpx1*cpy2-cpy1*cpx2;
                        module = r1r2x^2+r1r2y^2+r1r2z^2;
                        if (module<ep) || (r1<ep) || (r2<ep)
                            K1 = 0;
                        else
                            K1 = 1/module*(((-x21*cpx1-y21*cpy1)/r1)-((-x21*cpx2-y21*cpy2)/r2));
                        end
                    end
                    if d2(k,l) < ep
                        a2 = 0;
                        K2 = 0;
                    else
                        F = y32*e2 - x32*h2;
                        G = y32*e3 - x32*h3;
                        a2 = atan2(cpz*x32*(F*r3-G*r2), cpz^2*x32^2*r2*r3+F*G);
                        r2r3x = cpy2*cpz-cpz*cpy3;
                        r2r3y = -cpx2*cpz+cpz*cpx3;
                        r2r3z = cpx2*cpy3-cpy2*cpx3;
                        module = r2r3x^2+r2r3y^2+r2r3z^2;
                        if (module<ep) || (r2<ep) || (r3<ep)
                            K2 = 0;
                        else
                            K2 = 1/module*(((-x32*cpx2-y32*cpy2)/r2)-((-x32*cpx3-y32*cpy3)/r3));
                        end
                    end
                    if d3(k,l) < ep
                        a3 = 0;
                        K3 = 0;
                    else
                        F = y43*e3 - x43*h3;
                        G = y43*e4 - x43*h4;
                        a3 = atan2(cpz*x43*(F*r4-G*r3), cpz^2*x43^2*r3*r4+F*G);
                        r3r4x = cpy3*cpz-cpz*cpy4;
                        r3r4y = -cpx3*cpz+cpz*cpx4;
                        r3r4z = cpx3*cpy4-cpy3*cpx4;
                        module = r3r4x^2+r3r4y^2+r3r4z^2;
                        if (module<ep) || (r3<ep) || (r4<ep)
                            K3 = 0;
                        else
                            K3 = 1/module*(((-x43*cpx3-y43*cpy3)/r3)-((-x43*cpx4-y43*cpy4)/r4));
                        end
                    end
                    if d4(k,l) < ep
                        a4 = 0;
                        K4 = 0;
                    else
                        F = y14*e4 - x14*h4;
                        G = y14*e1 - x14*h1;
                        a4 = atan2(cpz*x14*(F*r1-G*r4), cpz^2*x14^2*r4*r1+F*G);
                        r4r1x = cpy4*cpz-cpz*cpy1;
                        r4r1y = -cpx4*cpz+cpz*cpx1;
                        r4r1z = cpx4*cpy1-cpy4*cpx1;
                        module = r4r1x^2+r4r1y^2+r4r1z^2;
                        if (module<ep) || (r4<ep) || (r1<ep)
                            K4 = 0;
                        else
                            K4 = 1/module*(((-x14*cpx4-y14*cpy4)/r4)-((-x14*cpx1-y14*cpy1)/r1));
                        end
                    end
                    uus = 0.079577471545948*(    ...
                        y21/d1(k,l)*reallog( (r1+r2-d1(k,l))/(r1+r2+d1(k,l)) )...
                        +   y32/d2(k,l)*reallog( (r2+r3-d2(k,l))/(r2+r3+d2(k,l)) )...
                        +   y43/d3(k,l)*reallog( (r3+r4-d3(k,l))/(r3+r4+d3(k,l)) )...
                        +   y14/d4(k,l)*reallog( (r4+r1-d4(k,l))/(r4+r1+d4(k,l)) ));
                    
                    vvs = -0.079577471545948*(    ...
                        x21/d1(k,l)*reallog( (r1+r2-d1(k,l))/(r1+r2+d1(k,l)) )...
                        +   x32/d2(k,l)*reallog( (r2+r3-d2(k,l))/(r2+r3+d2(k,l)) )...
                        +   x43/d3(k,l)*reallog( (r3+r4-d3(k,l))/(r3+r4+d3(k,l)) )...
                        +   x14/d4(k,l)*reallog( (r4+r1-d4(k,l))/(r4+r1+d4(k,l)) ));
                    
                    wws = 0.079577471545948*(a1+a2+a3+a4);
                    
                    uud = 0.079577471545948*(K1*r1r2x+K2*r2r3x+K3*r3r4x+K4*r4r1x);
                    
                    vvd = 0.079577471545948*(K1*r1r2y+K2*r2r3y+K3*r3r4y+K4*r4r1y);
                    
                    wwd = 0.079577471545948*(K1*r1r2z+K2*r2r3z+K3*r3r4z+K4*r4r1z);
  
                end
                % Project everything in the common axis
                us(I) = uus*l1(k,l) + vvs*t1(k,l) + wws*n1(k,l);
                vs(I) = uus*l2(k,l) + vvs*t2(k,l) + wws*n2(k,l);
                ws(I) = uus*l3(k,l) + vvs*t3(k,l) + wws*n3(k,l);
                ud(I) = uud*l1(k,l) + vvd*t1(k,l) + wwd*n1(k,l);
                vd(I) = uud*l2(k,l) + vvd*t2(k,l) + wwd*n2(k,l);
                wd(I) = uud*l3(k,l) + vvd*t3(k,l) + wwd*n3(k,l);
                
                
                if (isnan(us(I)) || isnan(vs(I)) || isnan(ws(I)) || isinf(us(I)) || isinf(vs(I)) || isinf(ws(I)))
                    disp('Found an issue in the calculation of induced velocities by the body on the wake')
                    pause(5)
                end
            end
        end
    end
end

% Now the effect of the wake on itself
for k = N+1:N+1+Nwakes
    for l = 1:M
        for i = N+2:N+2+Nwakes
            for j = 1:M+1
                I = I+1;
                xijcxkl = x(i,j)-cx(k,l);
                yijcykl = y(i,j)-cy(k,l);
                zijczkl = z(i,j)-cz(k,l);
                cpx = xijcxkl*l1(k,l) + yijcykl*l2(k,l) + zijczkl*l3(k,l);
                cpy = xijcxkl*t1(k,l) + yijcykl*t2(k,l) + zijczkl*t3(k,l);
                cpz = xijcxkl*n1(k,l) + yijcykl*n2(k,l) + zijczkl*n3(k,l);
                % if distance of panel from influenced point is greater
                % then product of longer diagonal and "far field" coefficient
                if realsqrt(cpx^2+cpy^2+cpz^2) > FF(k,l)
                    rad = cpx^2 + cpy^2 + cpz^2;
                    uud = 3*0.079577471545948*S(k,l)*cpx*cpz*rad^(-2.5);
                    vvd = 3*0.079577471545948*S(k,l)*cpy*cpz*rad^(-2.5);
                    wwd = -0.079577471545948*S(k,l)*(cpx^2+cpy^2-2*cpz^2)*rad^(-2.5);
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
                    x21 = x2(k,l)-x1(k,l);
                    x32 = x3(k,l)-x2(k,l);
                    x43 = x4(k,l)-x3(k,l);
                    x14 = x1(k,l)-x4(k,l);
                    y21 = y2(k,l)-y1(k,l);
                    y32 = y3(k,l)-y2(k,l);
                    y43 = y4(k,l)-y3(k,l);
                    y14 = y1(k,l)-y4(k,l);
                    if d1(k,l) < ep
                        K1 = 0;
                    else
                        r1r2x = cpy1*cpz-cpz*cpy2;
                        r1r2y = -cpx1*cpz+cpz*cpx2;
                        r1r2z = cpx1*cpy2-cpy1*cpx2;
                        module = r1r2x^2+r1r2y^2+r1r2z^2;
                        if (min(module*r1,module*r2)<ep2) || (r1<ep) || (r2<ep)
                            K1 = 0;
                            if (min(module*r1,module*r2)<ep2)
%                                 disp('WEI')
%                                 pause(5)
                            end
                        else
                            LL = x21^2+y21^2;
                            if Scope == 'predict'
                                rcore2 = 4*1.25643*1.488e-5*k*TimeStep;      % Vasistas model with nu of air computed at 18 °C
                            else
                                rcore2 = 4*1.25643*1.488e-5*(k+1)*TimeStep;
                            end    
                            Cv1 = ((LL*r1^2-(-x21*cpx1-y21*cpy1)^2)/LL)/(rcore2+(LL*r1^2-(-x21*cpx1-y21*cpy1)^2)/LL);
                            K1 = Cv1/module*(((-x21*cpx1-y21*cpy1)/r1)-((-x21*cpx2-y21*cpy2)/r2));
                        end
                    end
                    if d2(k,l) < ep
                        K2 = 0;
                    else
                        r2r3x = cpy2*cpz-cpz*cpy3;
                        r2r3y = -cpx2*cpz+cpz*cpx3;
                        r2r3z = cpx2*cpy3-cpy2*cpx3;
                        module = r2r3x^2+r2r3y^2+r2r3z^2;
                        if (min(module*r2,module*r3)<ep2) || (r2<ep) || (r3<ep)
                            K2 = 0;
                            if (min(module*r2,module*r3)<ep2)
%                                 disp('WEI')
%                                 pause(5)
                            end
                        else
                            LL = x32^2+y32^2;
                            if Scope == 'predict'
                                rcore2 = 4*1.25643*1.488e-5*k*TimeStep;      % Vasistas model with nu of air computed at 18 °C
                            else
                                rcore2 = 4*1.25643*1.488e-5*(k+1)*TimeStep;
                            end    
                            Cv2 = ((LL*r2^2-(-x32*cpx2-y32*cpy2)^2)/LL)/(rcore2+(LL*r2^2-(-x32*cpx2-y32*cpy2)^2)/LL);
                            K2 = Cv2/module*(((-x32*cpx2-y32*cpy2)/r2)-((-x32*cpx3-y32*cpy3)/r3));
                        end
                    end
                    if d3(k,l) < ep
                        K3 = 0;
                    else
                        r3r4x = cpy3*cpz-cpz*cpy4;
                        r3r4y = -cpx3*cpz+cpz*cpx4;
                        r3r4z = cpx3*cpy4-cpy3*cpx4;
                        module = r3r4x^2+r3r4y^2+r3r4z^2;
                        if (min(module*r3,module*r4)<ep2) || (r3<ep) || (r4<ep)
                            K3 = 0;
                            if (min(module*r3,module*r4)<ep2)
%                                 disp('WEI')
%                                 pause(5)
                            end
                        else
                            LL = x43^2+y43^2;
                            if Scope == 'predict'
                                rcore2 = 4*1.25643*1.488e-5*k*TimeStep;      % Vasistas model with nu of air computed at 18 °C
                            else
                                rcore2 = 4*1.25643*1.488e-5*(k+1)*TimeStep;
                            end    
                            Cv3 = ((LL*r3^2-(-x43*cpx3-y43*cpy3)^2)/LL)/(rcore2+(LL*r3^2-(-x43*cpx3-y43*cpy3)^2)/LL);
                            K3 = Cv3/module*(((-x43*cpx3-y43*cpy3)/r3)-((-x43*cpx4-y43*cpy4)/r4));
                        end
                    end
                    if d4(k,l) < ep
                        K4 = 0;
                    else
                        r4r1x = cpy4*cpz-cpz*cpy1;
                        r4r1y = -cpx4*cpz+cpz*cpx1;
                        r4r1z = cpx4*cpy1-cpy4*cpx1;
                        module = r4r1x^2+r4r1y^2+r4r1z^2;
                        if (min(module*r4,module*r1)<ep2) || (r4<ep) || (r1<ep)
                            K4 = 0;
                            if (min(module*r4,module*r1)<ep2)
%                                 disp('WEI')
%                                 pause(5)
                            end
                        else
                            LL = x14^2+y14^2;
                            if Scope == 'predict'
                                rcore2 = 4*1.25643*1.488e-5*k*TimeStep;      % Vasistas model with nu of air computed at 18 °C
                            else
                                rcore2 = 4*1.25643*1.488e-5*(k+1)*TimeStep;
                            end    
                            Cv4 = ((LL*r4^2-(-x14*cpx4-y14*cpy4)^2)/LL)/(rcore2+(LL*r4^2-(-x14*cpx4-y14*cpy4)^2)/LL);
                            K4 = Cv4/module*(((-x14*cpx4-y14*cpy4)/r4)-((-x14*cpx1-y14*cpy1)/r1));
                        end
                    end
                    
                    uud = 0.079577471545948*(K1*r1r2x+K2*r2r3x+K3*r3r4x+K4*r4r1x);
                    
                    vvd = 0.079577471545948*(K1*r1r2y+K2*r2r3y+K3*r3r4y+K4*r4r1y);
                    
                    wwd = 0.079577471545948*(K1*r1r2z+K2*r2r3z+K3*r3r4z+K4*r4r1z); 
                                 
                end
                % Project everything in the common axis
                ud(I) = uud*l1(k,l) + vvd*t1(k,l) + wwd*n1(k,l);
                vd(I) = uud*l2(k,l) + vvd*t2(k,l) + wwd*n2(k,l);
                wd(I) = uud*l3(k,l) + vvd*t3(k,l) + wwd*n3(k,l);
                
                if (isnan(ud(I)) || isnan(vd(I)) || isnan(wd(I)) || isinf(ud(I)) || isinf(vd(I)) || isinf(wd(I)))
                    disp('Found an issue in the calculation of induced velocities by the wake on the wake')
                    pause(5)
                end
                        
            end
        end
    end
end
end