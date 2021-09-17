function IDWout = IDW(x,y,z,xA,yA,zA)

% IDWout = zeros(length(x),length(xA));

x2 = repmat(x,1,length(xA));
y2 = repmat(y,1,length(yA));
z2 = repmat(z,1,length(zA));
xA2 = repmat(xA',length(x),1);
yA2 = repmat(yA',length(y),1);
zA2 = repmat(zA',length(z),1);
d = 1./((x2-xA2).^2+(y2-yA2).^2+(z2-zA2).^2); % 1/distance^2

dSum = sum(d,1);
dSumA = repmat(dSum,length(x),1);

IDWout = d./dSumA;