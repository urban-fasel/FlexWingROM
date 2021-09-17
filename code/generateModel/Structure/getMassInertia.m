function [m,I,cog] = getMassInertia(Nodes,Elements,thickness,density)

m = 0;
cog = zeros(1,3);
I = zeros(3,3);

for ii = 1:length(Elements)
    P1 = Nodes(Elements(ii,1),:);
    P2 = Nodes(Elements(ii,2),:);
    P3 = Nodes(Elements(ii,3),:);
    
    cogii = mean([P1;P2;P3]);
    A = 1/2*norm(cross(P2-P1,P3-P1));
    mii = A*thickness*density;
    I(1,1) = I(1,1) + mii*norm(cogii([2 3]))^2;
    I(2,2) = I(2,2) + mii*norm(cogii([1 3]))^2;
    I(3,3) = I(3,3) + mii*norm(cogii([1 2]))^2;
    cog = (cog*m + cogii*mii)/(m+mii);
    m = m + mii;
end

end

