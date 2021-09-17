function [G,Gs] = TPS(x,y,z,xA,yA,zA)

R = [ones(size(x,1),1),x,y,z];

x2 = repmat(x,1,length(x));
y2 = repmat(y,1,length(y));
z2 = repmat(z,1,length(z));
rSS = (x2-x2').^2+(y2-y2').^2+(z2-z2').^2; % rSS = r^2
K = rSS.*log(sqrt(rSS));
K(isnan(K)) = 0;

C = [zeros(4,4), R'; R, K];
Cinv = inv(C);

x3 = repmat(x',length(xA),1);
y3 = repmat(y',length(yA),1);
z3 = repmat(z',length(zA),1);
xA3 = repmat(xA,1,length(x));
yA3 = repmat(yA,1,length(y));
zA3 = repmat(zA,1,length(z));
rSA = (x3-xA3).^2+(y3-yA3).^2+(z3-zA3).^2; % rSA = r^2
Kk = rSA.*log(sqrt(rSA));
Kk(isnan(Kk)) = 0;
DKk = -2*(x3-xA3).*(1+log(rSA));
DKk(isnan(DKk)) = 0;

skI = [ones(length(xA),1),xA,yA,zA,Kk]*Cinv;
sk = skI(:,5:end);

DskI = [zeros(length(xA),1),-ones(length(xA),1),zeros(length(xA),1),zeros(length(xA),1),DKk]*Cinv;
Dsk = DskI(:,5:end);

sizeSK = size(x,1);
G = zeros(size(xA,1)*6,sizeSK*6);
Gs = zeros(size(xA,1),sizeSK);
for i = 1:size(xA,1)
    G(1+(i-1)*6,1:sizeSK) = sk(i,:);
    G(2+(i-1)*6,sizeSK+1:2*sizeSK) = sk(i,:);
    G(3+(i-1)*6,2*sizeSK+1:3*sizeSK) = sk(i,:);
    G(4+(i-1)*6,1:sizeSK) = Dsk(i,:);
    G(5+(i-1)*6,sizeSK+1:2*sizeSK) = Dsk(i,:);
    G(6+(i-1)*6,2*sizeSK+1:3*sizeSK) = Dsk(i,:);
    Gs(i,:) = sk(i,:);
end
