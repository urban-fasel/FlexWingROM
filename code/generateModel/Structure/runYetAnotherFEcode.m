function [M,K,wP] = runYetAnotherFEcode(wingProperties, wingDesign)

%% Split quad Elements to have only triangular elements
Elements_all = wingProperties.triaElements(:,2:end);
quadElements = wingProperties.quadElements(:,2:end); 
split_elements = zeros(2*size(quadElements,1),4);

j = 1;
i = 1;
while ( i <= size(quadElements,1))
    split_elements(j,:) = quadElements(i,1:4);
    split_elements(j+1,:) = quadElements(i,[1 2 4 5]);
    i = i+1;
    j = j+2;
end
Elements_all = [Elements_all; split_elements];


%% Assign nodes and elements
Nodes = wingProperties.FE_Grid;

% Assign elements, according to which body their nodes belong to
Elements = zeros(0,4);
for i = 1 : size(Elements_all,1)
    if (sum(Nodes (:,1) == Elements_all(i,2)) && sum(Nodes (:,1) == Elements_all(i,3)) && sum(Nodes (:,1) == Elements_all(i,4)))
        Elements = [Elements; Elements_all(i,:)];
    end
end

% Deleting the first row of Nodes{j} (which gave the node_id). Now the node at
% the i-th row of the array Nodes{j} has node_id = i. Since the node_ids have changed, 
% the node_ids in Elements{j} are changed accordingly
Element_grid = Elements(:,2:end);
for i = 1 : size(Nodes,1)
    Element_grid (Element_grid == Nodes(i,1)) = i;
end
Elements(:,2:end) = Element_grid;
Nodes = Nodes(:,2:end);


%% Define material and assign properties to mesh

%% material
for i = 1:length(wingDesign.matProp.E)
    myMaterial{i}  = KirchoffMaterial();
    set(myMaterial{i},'YOUNGS_MODULUS',wingDesign.matProp.E(i),'DENSITY',wingDesign.matProp.rho(i),'POISSONS_RATIO',wingDesign.matProp.nu(i));
end


%% mesh
myMesh = Mesh(Nodes);

mTot = 0; % wing mass
Itot = zeros(3,3); % wing inertia
cogTot = zeros(1,3); % wing center of gravity

% internal morphing structure
if wingDesign.nRibsC > 0
    EID = find(Elements(:,1) == wingProperties.isPIDout(1,1));
    Element_IS = Elements(EID,2:4);
    myElementConstructor =  @()TriShellElement(wingProperties.isPIDout(1,3), myMaterial{2}); 
    myMesh.create_elements_table(Element_IS,myElementConstructor);
    
    [m,I,cog] = getMassInertia(Nodes,Element_IS,wingProperties.isPIDout(1,3),myMaterial{2}.DENSITY);
    cogTot = (cogTot*mTot + cog*m)/(mTot+m);
    mTot = mTot + m;
    Itot = Itot + I;
end

% flexible skin
EID = find(Elements(:,1) == wingProperties.fsPIDout(1,1));
Element_FSp = Elements(EID,2:4);
myElementConstructor =  @()TriShellElement(wingProperties.fsPIDout(1,3), myMaterial{3});
myMesh.create_elements_table(Element_FSp,myElementConstructor);

[m,I,cog] = getMassInertia(Nodes,Element_FSp,wingProperties.fsPIDout(1,3),myMaterial{3}.DENSITY);
cogTot = (cogTot*mTot + cog*m)/(mTot+m);
mTot = mTot + m;
Itot = Itot + I;

% spars
EID = find(Elements(:,1) == wingProperties.sparPIDout(1,1));
Element_SPARp = Elements(EID,2:4);
myElementConstructor =  @()TriShellElement(wingProperties.sparPIDout(1,3), myMaterial{1});
myMesh.create_elements_table(Element_SPARp,myElementConstructor);

[m,I,cog] = getMassInertia(Nodes,Element_SPARp,wingProperties.sparPIDout(1,3),myMaterial{1}.DENSITY);
cogTot = (cogTot*mTot + cog*m)/(mTot+m);
mTot = mTot + m;
Itot = Itot + I;

% rigid ribs
EID = find(Elements(:,1) == wingProperties.rigidRibPIDout(1,1));
Element_RRp = Elements(EID,2:4);
myElementConstructor =  @()TriShellElement(wingProperties.rigidRibPIDout(1,3), myMaterial{1});
myMesh.create_elements_table(Element_RRp,myElementConstructor);

[m,I,cog] = getMassInertia(Nodes,Element_RRp,wingProperties.rigidRibPIDout(1,3),myMaterial{1}.DENSITY);
cogTot = (cogTot*mTot + cog*m)/(mTot+m);
mTot = mTot + m;
Itot = Itot + I;

% actuation
EID = find(Elements(:,1) == wingProperties.actuationPIDout(1,1));
Element_Ap = Elements(EID,2:4);
myElementConstructor =  @()TriShellElement(wingProperties.actuationPIDout(1,3), myMaterial{1});
myMesh.create_elements_table(Element_Ap,myElementConstructor);

[m,I,cog] = getMassInertia(Nodes,Element_Ap,wingProperties.actuationPIDout(1,3),myMaterial{1}.DENSITY);
cogTot = (cogTot*mTot + cog*m)/(mTot+m);
mTot = mTot + m;
Itot = Itot + I;

% skin
EID = find(Elements(:,1) == wingProperties.skinPIDout(1,1));
Element_Sp = Elements(EID,2:4);
myElementConstructor =  @()TriShellElement(wingProperties.skinPIDout(1,3), myMaterial{1});
myMesh.create_elements_table(Element_Sp,myElementConstructor);

[m,I,cog] = getMassInertia(Nodes,Element_Sp,wingProperties.skinPIDout(1,3),myMaterial{1}.DENSITY);
cogTot = (cogTot*mTot + cog*m)/(mTot+m);
mTot = mTot + m;
Itot = Itot + I;

% coordinate system located in wing leading edge, x points in front, y points in right wing, z points down  
wP.cog = [-cogTot(1), cogTot(3), -cogTot(2)]; % rotate cos
wP.m = mTot;
wP.I = [Itot(1,1) 0 0; 0 Itot(3,3) 0; 0 0 Itot(2,2)]; % rotate cos


%% Save nodes and elements for plotting mesh
if wingDesign.plot
    wP.Nodes = Nodes;
    wP.Element_IS = Element_IS;
    wP.Element_FSp = Element_FSp;
    wP.Element_SPARp = Element_SPARp;
    wP.Element_Ap = Element_Ap;
    wP.Element_Sp = Element_Sp;
    wP.Element_RRp = Element_RRp;
end

%% Assemble linear stiffness and mass
PlateAssembly = Assembly(myMesh);
K = PlateAssembly.stiffness_matrix();
M = PlateAssembly.mass_matrix();

K = sparse(K);
M = sparse(M);
