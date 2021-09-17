function plotModel(wingDesign, simParam, wingModelAero, wingModelStructure)

%% plot FE model

Nodes = wingModelStructure.wP.Nodes;
Element_IS = wingModelStructure.wP.Element_IS;
Element_FSp = wingModelStructure.wP.Element_FSp;
Element_SPARp = wingModelStructure.wP.Element_SPARp;
Element_Ap = wingModelStructure.wP.Element_Ap;
Element_Sp = wingModelStructure.wP.Element_Sp;
Element_RRp = wingModelStructure.wP.Element_RRp;

if wingDesign.nRibsC > 0
    figure
    PlotMesh(Nodes,Element_IS,0); % plot internal structure
    title('Mesh: compliant internal structure')
end

figure
PlotMesh(Nodes,Element_FSp,0); % plot flexible skin
title('Mesh: flexible skin')

figure
PlotMesh(Nodes,Element_SPARp,0); % plot spars
title('Mesh: spars')

figure
PlotMesh(Nodes,Element_RRp,0); % plot rigid ribs
title('Mesh: rigid ribs')

if wingDesign.nRibsC > 0
    figure
    PlotMesh(Nodes,Element_Ap,0); % plot actuation mounting
    title('Mesh: actuation mounting')
end

figure
PlotMesh(Nodes,Element_Sp,0); % plot wing skin
title('Mesh: wing skin')

figure
PlotMesh(Nodes,Element_SPARp,0); 
PlotMesh(Nodes,Element_RRp,0); 
if wingDesign.nRibsC > 0
    PlotMesh(Nodes,Element_Ap,0); 
    PlotMesh(Nodes,Element_IS,0); 
end
PlotMesh(Nodes,Element_FSp,0); 
title('Mesh: wing internal structure')
view([0, 89])

figure
offset = 0.2;
PlotMesh(Nodes + repmat([0 offset*0 0], length(Nodes), 1),Element_Sp,0); 
PlotMesh(Nodes + repmat([0 offset*2 0], length(Nodes), 1),Element_SPARp,0); 
PlotMesh(Nodes + repmat([0 offset*2 0], length(Nodes), 1),Element_RRp,0); 
if wingDesign.nRibsC > 0
    PlotMesh(Nodes + repmat([0 offset*2 0], length(Nodes), 1),Element_Ap,0); 
    PlotMesh(Nodes + repmat([0 offset*2 0], length(Nodes), 1),Element_IS,0); 
end
PlotMesh(Nodes + repmat([0 offset*1 0], length(Nodes), 1),Element_FSp,0); 
title('Mesh: flexible wing')




%% plot vibration and morphing deformation mode shapes
modePlotScale = 20;
for i = 1:(simParam.nmodes + 2)
    q_ASET = zeros(simParam.nmodes + 2,1);
    q_ASET(i) = 1;
    plotModeShape(wingDesign, simParam, wingModelAero, q_ASET, modePlotScale, i);
end

