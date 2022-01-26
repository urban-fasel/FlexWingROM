function plotModel(wingDesign, simParam, wingModelAero, wingModelStructure)

%% plot FE model

Nodes = wingModelStructure.wP.Nodes;
Element_IS = wingModelStructure.wP.Element_IS;
Element_FSp = wingModelStructure.wP.Element_FSp;
Element_SPARp = wingModelStructure.wP.Element_SPARp;
Element_Ap = wingModelStructure.wP.Element_Ap;
Element_Sp = wingModelStructure.wP.Element_Sp;
Element_RRp = wingModelStructure.wP.Element_RRp;


offset = 0.4;

figure

PlotMesh([Nodes(:,3),Nodes(:,1),Nodes(:,2)+offset],Element_Sp,0); % plot wing skin

PlotMesh([Nodes(:,3),Nodes(:,1),Nodes(:,2)],Element_SPARp,0); % plot spars
PlotMesh([Nodes(:,3),Nodes(:,1),Nodes(:,2)],Element_RRp,0); % plot rigid ribs
if wingDesign.nRibsC > 0
    PlotMesh([Nodes(:,3),Nodes(:,1),Nodes(:,2)],Element_Ap,0); % plot actuation mounting
    PlotMesh([Nodes(:,3),Nodes(:,1),Nodes(:,2)],Element_IS,0); % plot compliant ribs internal structure
    PlotMesh([Nodes(:,3),Nodes(:,1),Nodes(:,2)-offset],Element_FSp,0); % plot flexible skin, for morphing wing
else
    PlotMesh([Nodes(:,3),Nodes(:,1),Nodes(:,2)+offset],Element_FSp,0); % plot flexible skin, for non-morphing wing
end

set(gca, 'Projection','perspective')
% view([142 18])
% view([165 49])
view([68 11])
set(gca,'color','w')
title('Mesh: flexible wing')



%% plot vibration and morphing deformation mode shapes
modePlotScale = 20;
for i = 1:(simParam.nmodes + 2)
    q_ASET = zeros(simParam.nmodes + 2,1);
    q_ASET(i) = 1;
    plotModeShape(wingDesign, simParam, wingModelAero, q_ASET, modePlotScale, i);
end

