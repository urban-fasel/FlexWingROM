function wingModelStructure = generateWingModelStructure(wingDesign,simParam)

%% define material and property IDs
MID_Skin = 1; % wing skin material ID
MID_IS = 2; % morphing ribs internal structure material ID
MID_FS = 3; % flexible skin material ID
MID_Spar = MID_Skin;

PID_Skin = 1; % wing skin property ID
PID_FS = 2; % flexible skin property ID
PID_IS = 3; % morphing ribs internal structure property ID
PID_Spar = 4; % spar property ID


%% Create planform

% General settings
numRibs = 2*wingDesign.nRibs;
extraRibOffset = 0.005; % offset needed to generate mesh
ribPos = linspace(0, 1, wingDesign.nRibs); % position as a percentage of the non-rib space
span = wingDesign.span - wingDesign.ribWidth*(wingDesign.nRibs-2)*2 - 4*extraRibOffset; % intermediate wing span 
semispan = span/2; % semispan 

% Define planform parameters
ribsDR = struct('LEPos', mat2cell([zeros(1, wingDesign.nRibs); zeros(1, wingDesign.nRibs); ribPos*semispan], 3, ones(1, wingDesign.nRibs)));

for i = 1:wingDesign.nRibs-1
    x1 = ribsDR(i).LEPos([1 3]);
    x2 = ribsDR(i+1).LEPos([1 3]);
    r1 = wingDesign.chord;
    r2 = wingDesign.chord;
    
    [t1a, t2a, t1b, t2b] = circleCircleTangent(x1, x2, r1, r2);
    
    % choose between the two solutions
    if (t1a(1) < t1b(1)) && (t2a(1) < t2b(1))
        t1 = t1a;
        t2 = t2a;
    elseif (t1a(1) > t1b(1)) && (t2a(1) > t2b(1))
        t1 = t1b;
        t2 = t2b;
    else
        error('runIndividual_NASPan:PlanformCreation', 'Cannot create suitable planform - circle tangency failed!');
    end
    
    ribsDR(i).TEPos_tmp1(:,2) = [t1(1); 0; t1(2)]; % right side of rib i
    ribsDR(i+1).TEPos_tmp1(:,1) = [t2(1); 0; t2(2)]; % left side of rib i+1
end

ribsDR(1).TEPos = ribsDR(1).TEPos_tmp1(:,2);
ribsDR(end).TEPos = ribsDR(end).TEPos_tmp1(:,1);

% calculate the intersection between the various tangents
for i = 2:wingDesign.nRibs-1
    x1 = ribsDR(i-1).TEPos_tmp1([1 3], 2);
    x2 = ribsDR(i).TEPos_tmp1([1 3], 1);
    x3 = ribsDR(i).TEPos_tmp1([1 3], 2);
    x4 = ribsDR(i+1).TEPos_tmp1([1 3], 1);
    
    if abs(x2-x3) < eps(10) % points are basically coincident
        point = x2/2 + x3/2;
    else
        point = lineLineIntersection(x1, x2, x3, x4);
    end
    
    ribsDR(i).TEPos = [point(1); 0; point(2)]; % right side of rib i
end

% Calculate the position of the various ribs
for i = 1:wingDesign.nRibs

    ribDRTotalWidth = wingDesign.ribActiveWidth+2*wingDesign.ribWidth;
    rib1Ox1 = [ribsDR(i).LEPos(1); ribsDR(i).LEPos(3)] + ribDRTotalWidth/2*[0; -1];
    rib1Ox2 = [ribsDR(i).TEPos(1); ribsDR(i).TEPos(3)] + ribDRTotalWidth/2*[0; -1];
    rib2Ox1 = [ribsDR(i).LEPos(1); ribsDR(i).LEPos(3)] + ribDRTotalWidth/2*[0; 1];
    rib2Ox2 = [ribsDR(i).TEPos(1); ribsDR(i).TEPos(3)] + ribDRTotalWidth/2*[0; 1];
  
    rib1Ix1 = [ribsDR(i).LEPos(1); ribsDR(i).LEPos(3)] + wingDesign.ribActiveWidth/2*[0; -1];
    rib1Ix2 = [ribsDR(i).TEPos(1); ribsDR(i).TEPos(3)] + wingDesign.ribActiveWidth/2*[0; -1];
    rib2Ix1 = [ribsDR(i).LEPos(1); ribsDR(i).LEPos(3)] + wingDesign.ribActiveWidth/2*[0; 1];
    rib2Ix2 = [ribsDR(i).TEPos(1); ribsDR(i).TEPos(3)] + wingDesign.ribActiveWidth/2*[0; 1];

    % intersect the outermost part of the double rib with the LE and TE
    if i == 1
        LE1x1 = [ribsDR(i).LEPos(1); ribsDR(i).LEPos(3)];
        LE1x2 = [ribsDR(i+1).LEPos(1); ribsDR(i+1).LEPos(3)];
        TE1x1 = [ribsDR(i).TEPos(1); ribsDR(i).TEPos(3)];
        TE1x2 = [ribsDR(i+1).TEPos(1); ribsDR(i+1).TEPos(3)];
    else
        LE1x1 = [ribsDR(i-1).LEPos(1); ribsDR(i-1).LEPos(3)];
        LE1x2 = [ribsDR(i).LEPos(1); ribsDR(i).LEPos(3)];
        TE1x1 = [ribsDR(i-1).TEPos(1); ribsDR(i-1).TEPos(3)];
        TE1x2 = [ribsDR(i).TEPos(1); ribsDR(i).TEPos(3)];
    end
    if i == wingDesign.nRibs
        LE2x1 = [ribsDR(i-1).LEPos(1); ribsDR(i-1).LEPos(3)];
        LE2x2 = [ribsDR(i).LEPos(1); ribsDR(i).LEPos(3)];
        TE2x1 = [ribsDR(i-1).TEPos(1); ribsDR(i-1).TEPos(3)];
        TE2x2 = [ribsDR(i).TEPos(1); ribsDR(i).TEPos(3)];
    else
        LE2x1 = [ribsDR(i).LEPos(1); ribsDR(i).LEPos(3)];
        LE2x2 = [ribsDR(i+1).LEPos(1); ribsDR(i+1).LEPos(3)];
        TE2x1 = [ribsDR(i).TEPos(1); ribsDR(i).TEPos(3)];
        TE2x2 = [ribsDR(i+1).TEPos(1); ribsDR(i+1).TEPos(3)];
    end
    rib1OintLE = lineLineIntersection(rib1Ox1, rib1Ox2, LE1x1, LE1x2);
    rib1OintTE = lineLineIntersection(rib1Ox1, rib1Ox2, TE1x1, TE1x2);
    rib2OintLE = lineLineIntersection(rib2Ox1, rib2Ox2, LE2x1, LE2x2);
    rib2OintTE = lineLineIntersection(rib2Ox1, rib2Ox2, TE2x1, TE2x2);
    
    % the leading and trailing edge will be formed in addition by the segments between the newly found points
    
    % intersect the innermost part of the double rib with the LE and TE
    rib1IintLE = lineLineIntersection(rib1Ix1, rib1Ix2, rib1OintLE, rib2OintLE);
    rib1IintTE = lineLineIntersection(rib1Ix1, rib1Ix2, rib1OintTE, rib2OintTE);
    rib2IintLE = lineLineIntersection(rib2Ix1, rib2Ix2, rib1OintLE, rib2OintLE);
    rib2IintTE = lineLineIntersection(rib2Ix1, rib2Ix2, rib1OintTE, rib2OintTE);
    
    ribs(2*(i-1)+1).LE1 = [rib1OintLE(1); 0; rib1OintLE(2)];
    ribs(2*(i-1)+1).LE2 = [rib1IintLE(1); 0; rib1IintLE(2)];
    ribs(2*(i-1)+2).LE1 = [rib2IintLE(1); 0; rib2IintLE(2)];
    ribs(2*(i-1)+2).LE2 = [rib2OintLE(1); 0; rib2OintLE(2)];
    
    ribs(2*(i-1)+1).TE1 = [rib1OintTE(1); 0; rib1OintTE(2)];
    ribs(2*(i-1)+1).TE2 = [rib1IintTE(1); 0; rib1IintTE(2)];
    ribs(2*(i-1)+2).TE1 = [rib2IintTE(1); 0; rib2IintTE(2)];
    ribs(2*(i-1)+2).TE2 = [rib2OintTE(1); 0; rib2OintTE(2)];
end

numRibs = numRibs + 1;

% generate an additional rib at the root
LE1x1 = [ribs(1).LE1(1); ribs(1).LE1(3)];
LE1x2 = [ribs(1).LE2(1); ribs(1).LE2(3)];
TE1x1 = [ribs(1).TE1(1); ribs(1).TE1(3)];
TE1x2 = [ribs(1).TE2(1); ribs(1).TE2(3)];

% get the wing farthest point to the root
zClosestToRoot = min([ribs(1).LE1(3) ribs(1).TE1(3)]);

rib0x1 = [0; zClosestToRoot-extraRibOffset];
rib0x2 = [1; zClosestToRoot-extraRibOffset];
rib0intLE = lineLineIntersection(rib0x1, rib0x2, LE1x1, LE1x2);
rib0intTE = lineLineIntersection(rib0x1, rib0x2, TE1x1, TE1x2);


ribs(2:end+1) = ribs(1:end);
ribs(1).LE2 = [rib0intLE(1); 0; rib0intLE(2)];
ribs(1).LE1 = ribs(1).LE2;
ribs(1).TE2 = [rib0intTE(1); 0; rib0intTE(2)];
ribs(1).TE1 = ribs(1).TE2;

numRibs = numRibs + 1;

% generate an additional rib at the root
LE1x1 = [ribs(end).LE1(1); ribs(end).LE1(3)];
LE1x2 = [ribs(end).LE2(1); ribs(end).LE2(3)];
TE1x1 = [ribs(end).TE1(1); ribs(end).TE1(3)];
TE1x2 = [ribs(end).TE2(1); ribs(end).TE2(3)];

% get the wing farthest point to the tip
zClosestToTip = max([ribs(end).LE2(3) ribs(end).TE2(3)]);

ribNx1 = [0; zClosestToTip+extraRibOffset];
ribNx2 = [1; zClosestToTip+extraRibOffset];
ribNintLE = lineLineIntersection(ribNx1, ribNx2, LE1x1, LE1x2);
ribNintTE = lineLineIntersection(ribNx1, ribNx2, TE1x1, TE1x2);

nRibsp1 = numel(ribs)+1;
ribs(nRibsp1).LE2 = [ribNintLE(1); 0; ribNintLE(2)];
ribs(nRibsp1).LE1 = ribs(nRibsp1).LE2;
ribs(nRibsp1).TE2 = [ribNintTE(1); 0; ribNintTE(2)];
ribs(nRibsp1).TE1 = ribs(nRibsp1).TE2;

% Shift all ribs so that the root LE is at 0, 0
ribsShift = zeros(3, 1);

if ribs(2).LE1(3) < ribs(2).TE1(3)
    ribsShift(1) = ribs(1).LE1(1);
    ribsShift(3) = ribs(2).LE1(3);
else
    ribsShift(1) = ribs(1).LE1(1);
    ribsShift(3) = ribs(2).TE1(3);
end

for i = 1:numel(ribs)
    ribs(i).LE1 = ribs(i).LE1 - ribsShift;
    ribs(i).TE1 = ribs(i).TE1 - ribsShift;
    ribs(i).LE2 = ribs(i).LE2 - ribsShift;
    ribs(i).TE2 = ribs(i).TE2 - ribsShift;
end

% Calculate/update resulting parameters
Z_ribPos = zeros(size(ribs));
for i = 1:numel(ribs)
    Z_ribPos(i) = .5*ribs(i).LE1(3) + .5*ribs(i).LE2(3);
end


%% Create FE Models
rib_width = zeros(numRibs, 1);

%% FE Settings
nodes_mult = 5000; 
elems_mult = 5000; 
props_mult = 100;

FE_SkinNodeListU = [];
FE_SkinNodeListD = [];


%% Material / properties
domain_segments_PID_extra = zeros(1, numRibs); 

PID_Skin_up = PID_Skin*ones(9,(2*wingDesign.nRibs+1));
PID_Skin_dn = PID_Skin*ones(9,(2*wingDesign.nRibs+1));
PID_Skin_dn(6,:) = PID_FS;

domain_segments_PID_layup = PID_Skin*ones(numRibs,7);
domain_segments_PID_layup(:,5) = PID_Spar;
domain_segments_PID_layup(:,1) = 10; % define lower skin no corrugion, will be set to PID_Skin inside of createRibStructure

num_sdv_inner = 8;
isdv_spars = [4];
isdv_spar2 = [3];

%% Create structure
domain_settings.elem_min_number = 1;
domain_settings.elem_min_number_ext = 1;
domain_settings.elem_min_number_spanwise = 2;
domain_settings.elem_rib_min_number_spanwise = 2; % minimum number of elements on the rib
domain_settings.elem_max_size = simParam.elem_max_sizeC;
domain_settings.elem_max_size_ext = simParam.elem_max_sizeC;
domain_settings.elem_max_size_LE = simParam.elem_max_sizeC_LE;
domain_settings.elem_max_size_spanwise = simParam.elem_max_sizeS;
domain_settings.elem_rib_max_size_spanwise = simParam.elem_max_sizeS; % maximum element size on the rib

rib_numElemThickness = zeros(numRibs, 1);

meshOut.cquad = [];

for iRib = 1:numRibs
    % actuation motor actuation parameters and necessary nodes and positions
    rib(iRib).actuation.param_spar_y = wingDesign.actuatorFrontY;
    rib(iRib).actuation.param_connector_x = wingDesign.actuatorRearX;
    rib(iRib).actuation.param_connector_y = wingDesign.actuatorRearY;
    
    rib(iRib).actuation.sparNode_ID = [];
    rib(iRib).actuation.sparNode_Up_ID = [];
    rib(iRib).actuation.sparNode_Dn_ID = [];
    rib(iRib).actuation.sparNode_Pos = [];
    rib(iRib).actuation.sparNode_Up_Pos = [];
    rib(iRib).actuation.sparNode_Dn_Pos = [];
    
    rib(iRib).actuation.sdvNode_Up_ID = [];
    rib(iRib).actuation.sdvNode_Dn_ID = [];
    rib(iRib).actuation.sdvNode_Up_Pos = [];
    rib(iRib).actuation.sdvNode_Dn_Pos = [];
    
    rib(iRib).actuation.location_spar_ID = [];
    rib(iRib).actuation.location_spar_Pos = [];
	rib(iRib).actuation.sdvNode_Up_Bk_Pos = [];
	rib(iRib).actuation.sdvNode_Dn_Bk_Pos = [];
    rib(iRib).actuation.location_spar_Pos_prov = [];
    rib(iRib).actuation.location_sdv_ID = [];
    rib(iRib).actuation.location_sdv_Pos = [];
    rib(iRib).actuation.location_sdv_Pos_prov = [];
    
    isdv_actuation_1 = 0;
    isdv_actuation_2 = 0;
    sdvs_dist_act_connection = 0.01;
    

	domain_settings.domain_segments_thickness = repmat(0.002, 1, 7); % values are unused, but a 7-elements vector is needed!
	
	rib(iRib).targetShape.X_start = wingDesign.rearSpar*wingDesign.chord;
    rib(iRib).targetShape.X_trans = rib(iRib).targetShape.X_start + (rib(iRib).actuation.param_connector_x*(wingDesign.chord-rib(iRib).targetShape.X_start))*wingDesign.corrugationStart;
    rib(iRib).targetShape.X_trans2 = rib(iRib).targetShape.X_start + (rib(iRib).actuation.param_connector_x*(wingDesign.chord-rib(iRib).targetShape.X_start))*wingDesign.corrugationEnd;
    rib(iRib).targetShape.X_stop = wingDesign.chord;
    actuation_spar_dist = wingDesign.actuatorFrontX*(rib(iRib).actuation.param_connector_x*(wingDesign.chord-rib(iRib).targetShape.X_start));
    
	domain_settings.bounds = [	0, 0; 0.75, 0.1; 1, 0.2];
	
	rib_width(iRib) = wingDesign.ribWidth;
	rib_numElemThickness(iRib) = ceil(rib_width(iRib)/domain_settings.elem_rib_max_size_spanwise);
	if rib_numElemThickness(iRib) < domain_settings.elem_min_number_spanwise
        rib_numElemThickness(iRib) = domain_settings.elem_rib_min_number_spanwise; 
    end
	
	if iRib == 1
		rib_width(iRib) = 0;
		rib_numElemThickness(iRib) = 0;
	end
	if iRib == numRibs
		rib_width(iRib) = 0;
		rib_numElemThickness(iRib) = 0;
	end
     
	voronoi_settings.sites_X = wingDesign.voronoiXY(1:12);
	voronoi_settings.sites_Y = wingDesign.voronoiXY(13:24);
	voronoi_settings.sites_thickness = wingDesign.voronoiT*ones(12,1);
 
    airfoilShape = load(sprintf('airfoils/%s.dat',wingDesign.airfoil));
    rib(iRib).targetShape.nodes = [airfoilShape(1:100,1), airfoilShape(1:100,2), airfoilShape(100:199,1), airfoilShape(100:199,2)];

	rib(iRib).targetShape.nodes = rib(iRib).targetShape.nodes * wingDesign.chord;
	
	rib(iRib).sparPos = rib(iRib).targetShape.X_start;
    
	nodeOffset = nodes_mult*iRib;
	elemOffset = elems_mult*iRib;

    if iRib == 1 || iRib == numRibs || iRib <= (2*(wingDesign.nRibs-wingDesign.nRibsC)+1)
        createOnlyOuterNodes = true;
    else
        createOnlyOuterNodes = false;
    end
    
    % define specific PID
    domain_segments_PID = domain_segments_PID_layup(iRib,:);

	[rib(iRib).FE_nodes_3D, rib(iRib).FE_nodes_ext_up_sortedCCW, rib(iRib).FE_nodes_ext_dn_sortedCCW, rib(iRib).FE_nodes_spar, rib(iRib).FE_delta_iNodes_EdgeToEdge, rib(iRib).corner_nodes, rib(iRib).FE_Grid, rib(iRib).spar2, meshOut] ...
        = createRibStructure(meshOut, wingDesign, rib(iRib).targetShape, voronoi_settings, domain_settings, rib_width(iRib), nodeOffset, elemOffset, Z_ribPos(iRib), rib_numElemThickness(iRib), domain_segments_PID + domain_segments_PID_extra(iRib), createOnlyOuterNodes, PID_FS, PID_IS);
	
end


%% Create 3D structure
nodeOffset = nodes_mult*(numRibs+1);
elemOffset = elems_mult*(numRibs+1);

num_elem_span = zeros(1, numRibs-1);
for iRib = 1:numRibs-1
	num_elem_span(iRib) = ceil((Z_ribPos(iRib+1) - Z_ribPos(iRib))/domain_settings.elem_max_size_spanwise);
	if num_elem_span(iRib) < domain_settings.elem_min_number_spanwise 
        num_elem_span(iRib) = domain_settings.elem_min_number_spanwise; 
    end
end

% chordwise subdivision of wing: sets spar positions and can be used to define chordwise changing layup thickness
for iRib = 1:numRibs
	nsdv = 10; % number of subdivisions
	
	sdvPos = zeros(nsdv, 9); % 3 rows for the point on the left side of the rib, 3 rows for the one on the right
    sdvPos_rot = sdvPos;
	
    % sdv pos upscaled wing
    sdvPos(1,1) = .10*rib(iRib).sparPos; % wingbox sdv 1
    sdvPos(2,1) = .15*rib(iRib).sparPos; % wingbox sdv 2
	sdvPos(3,1) = wingDesign.frontSparX*rib(iRib).sparPos; % front spar
	sdvPos(4,1) = rib(iRib).sparPos; % <- spar
	sdvPos(5,1) = rib(iRib).targetShape.X_trans; % first sdv between spar and actuation    
	sdvPos(6,1) = rib(iRib).targetShape.X_trans2; % second sdv between spar and actuation 
    
    
    sdvPos(9,1) = sdvPos(4,1) + rib(iRib).actuation.param_connector_x*(wingDesign.chord-rib(iRib).sparPos); % sdv at actuation
    sdvPos(8,1) = sdvPos(9,1) - sdvs_dist_act_connection; % sdv in front of actuation    
    sdvPos(7,1) = sdvPos(8,1) - 0.01; % sdv in front of actuation  
    
    sdvPos(10,1) = inf; % <- Fake sdv

    isdv_actuation_1 = 7;
    isdv_actuation_2 = 8;
   
	sdvPos(:,3) = -0.5*rib_width(iRib);
	
	sdvPos(:,4) = sdvPos(:,1);
	sdvPos(:,6) = 0.5*rib_width(iRib);
	
	sdvPos_rot(:,[1:3]) = [sdvPos(:,1), sdvPos(:,2), sdvPos(:,3)];
	sdvPos_rot(:,[4:6]) = [sdvPos(:,4), sdvPos(:,5), sdvPos(:,6)];

	sdvPos_rot(:,[1:3]) = sdvPos_rot(:,[1:3]);
	sdvPos_rot(:,[4:6]) = sdvPos_rot(:,[4:6]);

    sdvPos_rot(:,1) = sdvPos_rot(:,1);
	sdvPos_rot(:,4) = sdvPos_rot(:,4);
	
	sdvPos_rot(isinf(sdvPos(:,1)),:) = inf;
	
	for j = 1:nsdv
		sdv(j).pos(iRib, 1) = sdvPos_rot(j, 1);
		sdv(j).pos(iRib, 2) = sdvPos_rot(j, 4);
	end
	
end

skinPanels_FE_Grid = [];
skinPanels_nodes.Lbound_up = NaN(numRibs-1, num_sdv_inner+2);
skinPanels_nodes.Lbound_dn = NaN(numRibs-1, num_sdv_inner+2);
skinPanels_nodes.Ubound_up = NaN(numRibs-1, num_sdv_inner+2);
skinPanels_nodes.Ubound_dn = NaN(numRibs-1, num_sdv_inner+2);

meshOut.ctria = [];

for iRib = 1:numRibs-1
	% create subdivisions along the chord for each station along the wingspan (in case we want sdvs, we need more zones)
	num_subdiv = num_sdv_inner+1;
	for iSubdiv = 1:num_subdiv+1
		dNodesEE = rib(iRib).FE_delta_iNodes_EdgeToEdge;
		switch iSubdiv
			case 1
				spanDiv(iRib).edgeSpanwise(iSubdiv).top.node1 = rib(iRib).FE_nodes_ext_up_sortedCCW(end);
				spanDiv(iRib).edgeSpanwise(iSubdiv).top.iNode1 = length(rib(iRib).FE_nodes_ext_up_sortedCCW);
				spanDiv(iRib).edgeSpanwise(iSubdiv).top.node2 = rib(iRib+1).FE_nodes_ext_up_sortedCCW(end);
				spanDiv(iRib).edgeSpanwise(iSubdiv).top.iNode2 = length(rib(iRib+1).FE_nodes_ext_up_sortedCCW);
				
				spanDiv(iRib).edgeSpanwise(iSubdiv).bot.node1 = rib(iRib).FE_nodes_ext_dn_sortedCCW(1);
				spanDiv(iRib).edgeSpanwise(iSubdiv).bot.iNode1 = 1;
				spanDiv(iRib).edgeSpanwise(iSubdiv).bot.node2 = rib(iRib+1).FE_nodes_ext_dn_sortedCCW(1);
				spanDiv(iRib).edgeSpanwise(iSubdiv).bot.iNode2 = 1;
			case num_subdiv+1
				spanDiv(iRib).edgeSpanwise(iSubdiv).top.node1 = rib(iRib).FE_nodes_ext_up_sortedCCW(1);
				spanDiv(iRib).edgeSpanwise(iSubdiv).top.iNode1 = 1;
				spanDiv(iRib).edgeSpanwise(iSubdiv).top.node2 = rib(iRib+1).FE_nodes_ext_up_sortedCCW(1);
				spanDiv(iRib).edgeSpanwise(iSubdiv).top.iNode2 = 1;
				spanDiv(iRib).edgeSpanwise(iSubdiv).bot.node1 = rib(iRib).FE_nodes_ext_dn_sortedCCW(end);
				spanDiv(iRib).edgeSpanwise(iSubdiv).bot.iNode1 = length(rib(iRib).FE_nodes_ext_dn_sortedCCW);
				spanDiv(iRib).edgeSpanwise(iSubdiv).bot.node2 = rib(iRib+1).FE_nodes_ext_dn_sortedCCW(end);
				spanDiv(iRib).edgeSpanwise(iSubdiv).bot.iNode2 = length(rib(iRib+1).FE_nodes_ext_dn_sortedCCW);
			otherwise
				% find the closest node
				[~, iNode] = min(abs(rib(iRib).FE_nodes_3D(rib(iRib).FE_nodes_ext_up_sortedCCW(:) + dNodesEE, 1) - sdv(iSubdiv-1).pos(iRib, 2)));
				spanDiv(iRib).edgeSpanwise(iSubdiv).top.iNode1 = iNode;
				% CHECK: are the starting and ending nodes on the short edges coincident? Yes -> move by one
				if spanDiv(iRib).edgeSpanwise(iSubdiv).top.iNode1 == spanDiv(iRib).edgeSpanwise(iSubdiv-1).top.iNode1
					spanDiv(iRib).edgeSpanwise(iSubdiv).top.iNode1 = spanDiv(iRib).edgeSpanwise(iSubdiv).top.iNode1 - 1;
				end
				spanDiv(iRib).edgeSpanwise(iSubdiv).top.node1 = rib(iRib).FE_nodes_ext_up_sortedCCW(iNode);
				
				[~, iNode] = min(abs(rib(iRib+1).FE_nodes_3D(rib(iRib+1).FE_nodes_ext_up_sortedCCW(:), 1) - sdv(iSubdiv-1).pos(iRib+1, 1)));
				spanDiv(iRib).edgeSpanwise(iSubdiv).top.iNode2 = iNode;
				% CHECK: are the starting and ending nodes on the short edges coincident? Yes -> move by one
				if spanDiv(iRib).edgeSpanwise(iSubdiv).top.iNode2 == spanDiv(iRib).edgeSpanwise(iSubdiv-1).top.iNode2
					spanDiv(iRib).edgeSpanwise(iSubdiv).top.iNode2 = spanDiv(iRib).edgeSpanwise(iSubdiv).top.iNode2 - 1;
				end
				spanDiv(iRib).edgeSpanwise(iSubdiv).top.node2 = rib(iRib+1).FE_nodes_ext_up_sortedCCW(iNode);
				
				[~, iNode] = min(abs(rib(iRib).FE_nodes_3D(rib(iRib).FE_nodes_ext_dn_sortedCCW(:) + dNodesEE, 1) - sdv(iSubdiv-1).pos(iRib, 2)));
				spanDiv(iRib).edgeSpanwise(iSubdiv).bot.iNode1 = iNode;
				% CHECK: are the starting and ending nodes on the short edges coincident? Yes -> move by one
				if spanDiv(iRib).edgeSpanwise(iSubdiv).bot.iNode1 == spanDiv(iRib).edgeSpanwise(iSubdiv-1).bot.iNode1
					spanDiv(iRib).edgeSpanwise(iSubdiv).bot.iNode1 = spanDiv(iRib).edgeSpanwise(iSubdiv).bot.iNode1 + 1;
				end
				spanDiv(iRib).edgeSpanwise(iSubdiv).bot.node1 = rib(iRib).FE_nodes_ext_dn_sortedCCW(iNode);
				
				[~, iNode] = min(abs(rib(iRib+1).FE_nodes_3D(rib(iRib+1).FE_nodes_ext_dn_sortedCCW(:), 1) - sdv(iSubdiv-1).pos(iRib+1, 1)));
				spanDiv(iRib).edgeSpanwise(iSubdiv).bot.iNode2 = iNode;
				% CHECK: are the starting and ending nodes on the short edges coincident? Yes -> move by one
				if spanDiv(iRib).edgeSpanwise(iSubdiv).bot.iNode2 == spanDiv(iRib).edgeSpanwise(iSubdiv-1).bot.iNode2
					spanDiv(iRib).edgeSpanwise(iSubdiv).bot.iNode2 = spanDiv(iRib).edgeSpanwise(iSubdiv).bot.iNode2 + 1;
				end
				spanDiv(iRib).edgeSpanwise(iSubdiv).bot.node2 = rib(iRib+1).FE_nodes_ext_dn_sortedCCW(iNode);
		end
	end
	for iSubdiv = 1:num_subdiv
		dNodesEE = rib(iRib).FE_delta_iNodes_EdgeToEdge;
		if 1 == iSubdiv
			% we need to create the nodes on one additional edge
			perc_nodes_edge4 = linspace(0, 1, num_elem_span(iRib)+1);
			% NOTE: first and last node coincide with rib nodes
			node1Upos = rib(iRib).FE_nodes_3D(spanDiv(iRib).edgeSpanwise(iSubdiv).top.node1 + dNodesEE, :);
			node2Upos = rib(iRib+1).FE_nodes_3D(spanDiv(iRib).edgeSpanwise(iSubdiv).top.node2, :);
			node1Dpos = rib(iRib).FE_nodes_3D(spanDiv(iRib).edgeSpanwise(iSubdiv).bot.node1 + dNodesEE, :);
			node2Dpos = rib(iRib+1).FE_nodes_3D(spanDiv(iRib).edgeSpanwise(iSubdiv).bot.node2, :);
			
			% check: are the upper and lower panel connected on edge4? If so, don't create the lower points...
			%%%%%%%%%%%% ################# NEED TO DO THE SAME FOR EDGE2
			if spanDiv(iRib).edgeSpanwise(iSubdiv).top.node1 == spanDiv(iRib).edgeSpanwise(iSubdiv).bot.node1
				spanDiv(iRib).wingSubdiv(iSubdiv).edge4U_nodes = (repmat(node1Upos, num_elem_span(iRib)-1, 1) + ((node2Upos-node1Upos)'*perc_nodes_edge4(2:end-1))');
				spanDiv(iRib).wingSubdiv(iSubdiv).edge4D_nodes = spanDiv(iRib).wingSubdiv(iSubdiv).edge4U_nodes;

				nodeID_edge4U = nodeOffset + (1:length(perc_nodes_edge4));
				nodeID_edge4D = nodeID_edge4U;
				nodeOffset = nodeOffset + length(perc_nodes_edge4);
				
				skinPanels_FE_Grid = [skinPanels_FE_Grid; [nodeID_edge4U(2:end-1)', spanDiv(iRib).wingSubdiv(iSubdiv).edge4U_nodes(:,1), spanDiv(iRib).wingSubdiv(iSubdiv).edge4U_nodes(:,2), spanDiv(iRib).wingSubdiv(iSubdiv).edge4U_nodes(:,3)]];
				
				spanDiv(iRib).wingSubdiv(iSubdiv).sdv4U = nodeID_edge4U(2:end-1);
				spanDiv(iRib).wingSubdiv(iSubdiv).sdv4D = NaN*nodeID_edge4U(2:end-1);
			else
				spanDiv(iRib).wingSubdiv(iSubdiv).edge4U_nodes = (repmat(node1Upos, num_elem_span(iRib)-1, 1) + ((node2Upos-node1Upos)'*perc_nodes_edge4(2:end-1))');
				spanDiv(iRib).wingSubdiv(iSubdiv).edge4D_nodes = (repmat(node1Dpos, num_elem_span(iRib)-1, 1) + ((node2Dpos-node1Dpos)'*perc_nodes_edge4(2:end-1))');

				nodeID_edge4U = nodeOffset + (1:length(perc_nodes_edge4));
				nodeID_edge4D = nodeOffset + length(perc_nodes_edge4) + (1:length(perc_nodes_edge4));
				nodeOffset = nodeOffset + 2*length(perc_nodes_edge4);
				
				skinPanels_FE_Grid = [skinPanels_FE_Grid; ...
					[nodeID_edge4U(2:end-1)', spanDiv(iRib).wingSubdiv(iSubdiv).edge4U_nodes(:,1), spanDiv(iRib).wingSubdiv(iSubdiv).edge4U_nodes(:,2), spanDiv(iRib).wingSubdiv(iSubdiv).edge4U_nodes(:,3)]; ...
					[nodeID_edge4D(2:end-1)', spanDiv(iRib).wingSubdiv(iSubdiv).edge4D_nodes(:,1), spanDiv(iRib).wingSubdiv(iSubdiv).edge4D_nodes(:,2), spanDiv(iRib).wingSubdiv(iSubdiv).edge4D_nodes(:,3)]];
				
				spanDiv(iRib).wingSubdiv(iSubdiv).sdv4U = nodeID_edge4U(2:end-1);
				spanDiv(iRib).wingSubdiv(iSubdiv).sdv4D = nodeID_edge4D(2:end-1);
			end
		else
			% copy from previous run
			perc_nodes_edge4 = perc_nodes_edge2;
			
			nodeID_edge4U = (nodeID_edge2U);
			nodeID_edge4D = (nodeID_edge2D);
			
			spanDiv(iRib).wingSubdiv(iSubdiv).edge4U_nodes = spanDiv(iRib).wingSubdiv(iSubdiv-1).edge2U_nodes;
			spanDiv(iRib).wingSubdiv(iSubdiv).edge4D_nodes = spanDiv(iRib).wingSubdiv(iSubdiv-1).edge2D_nodes;
			
			spanDiv(iRib).wingSubdiv(iSubdiv).sdv4U = NaN*nodeID_edge4U(2:end-1);
			spanDiv(iRib).wingSubdiv(iSubdiv).sdv4D = NaN*nodeID_edge4D(2:end-1);
		end
		
		perc_nodes_edge2 = linspace(0, 1, num_elem_span(iRib)+1);
		% NOTE: first and last node coincide with rib nodes
		node1Upos = rib(iRib).FE_nodes_3D(spanDiv(iRib).edgeSpanwise(iSubdiv+1).top.node1 + dNodesEE, :);
		node2Upos = rib(iRib+1).FE_nodes_3D(spanDiv(iRib).edgeSpanwise(iSubdiv+1).top.node2, :);
		node1Dpos = rib(iRib).FE_nodes_3D(spanDiv(iRib).edgeSpanwise(iSubdiv+1).bot.node1 + dNodesEE, :);
		node2Dpos = rib(iRib+1).FE_nodes_3D(spanDiv(iRib).edgeSpanwise(iSubdiv+1).bot.node2, :);
		
		% check: are the upper and lower panel connected on edge2? If so, don't create the lower points...
		if spanDiv(iRib).edgeSpanwise(iSubdiv+1).top.node1 == spanDiv(iRib).edgeSpanwise(iSubdiv+1).bot.node1
			spanDiv(iRib).wingSubdiv(iSubdiv).edge2U_nodes = (repmat(node1Upos, num_elem_span(iRib)-1, 1) + ((node2Upos-node1Upos)'*perc_nodes_edge2(2:end-1))');
			spanDiv(iRib).wingSubdiv(iSubdiv).edge2D_nodes = spanDiv(iRib).wingSubdiv(iSubdiv).edge2U_nodes;

			nodeID_edge2U = nodeOffset + (1:length(perc_nodes_edge2));
			nodeID_edge2D = nodeID_edge2U;
			nodeOffset = nodeOffset + length(perc_nodes_edge2);
			
			skinPanels_FE_Grid = [skinPanels_FE_Grid; [nodeID_edge2U(2:end-1)', spanDiv(iRib).wingSubdiv(iSubdiv).edge2U_nodes(:,1), spanDiv(iRib).wingSubdiv(iSubdiv).edge2U_nodes(:,2), spanDiv(iRib).wingSubdiv(iSubdiv).edge2U_nodes(:,3)]];
			
			spanDiv(iRib).wingSubdiv(iSubdiv).sdv2U = nodeID_edge2U(2:end-1);
			spanDiv(iRib).wingSubdiv(iSubdiv).sdv2D = NaN*nodeID_edge2U(2:end-1);
		else
			spanDiv(iRib).wingSubdiv(iSubdiv).edge2U_nodes = (repmat(node1Upos, num_elem_span(iRib)-1, 1) + ((node2Upos-node1Upos)'*perc_nodes_edge2(2:end-1))');
			spanDiv(iRib).wingSubdiv(iSubdiv).edge2D_nodes = (repmat(node1Dpos, num_elem_span(iRib)-1, 1) + ((node2Dpos-node1Dpos)'*perc_nodes_edge2(2:end-1))');

			nodeID_edge2U = nodeOffset + (1:length(perc_nodes_edge2));
			nodeID_edge2D = nodeOffset + length(perc_nodes_edge2) + (1:length(perc_nodes_edge2));
			nodeOffset = nodeOffset + 2*length(perc_nodes_edge2);
			
			skinPanels_FE_Grid = [skinPanels_FE_Grid; ...
				[nodeID_edge2U(2:end-1)', spanDiv(iRib).wingSubdiv(iSubdiv).edge2U_nodes(:,1), spanDiv(iRib).wingSubdiv(iSubdiv).edge2U_nodes(:,2), spanDiv(iRib).wingSubdiv(iSubdiv).edge2U_nodes(:,3)]; ...
				[nodeID_edge2D(2:end-1)', spanDiv(iRib).wingSubdiv(iSubdiv).edge2D_nodes(:,1), spanDiv(iRib).wingSubdiv(iSubdiv).edge2D_nodes(:,2), spanDiv(iRib).wingSubdiv(iSubdiv).edge2D_nodes(:,3)]];
			
			spanDiv(iRib).wingSubdiv(iSubdiv).sdv2U = nodeID_edge2U(2:end-1);
			spanDiv(iRib).wingSubdiv(iSubdiv).sdv2D = nodeID_edge2D(2:end-1);
		end
		
		% Adjust node IDs at the ends of the segment to reflect GLOBAL node ID
		nodeID_edge4U(1) = spanDiv(iRib).edgeSpanwise(iSubdiv).top.node1 + dNodesEE + nodes_mult*iRib;
		nodeID_edge4U(end) = spanDiv(iRib).edgeSpanwise(iSubdiv).top.node2 + nodes_mult*(iRib+1);
		nodeID_edge4D(1) = spanDiv(iRib).edgeSpanwise(iSubdiv).bot.node1 + dNodesEE + nodes_mult*iRib;
		nodeID_edge4D(end) = spanDiv(iRib).edgeSpanwise(iSubdiv).bot.node2 + nodes_mult*(iRib+1);
		
		nodeID_edge2U(1) = spanDiv(iRib).edgeSpanwise(iSubdiv+1).top.node1 + dNodesEE + nodes_mult*iRib;
		nodeID_edge2U(end) = spanDiv(iRib).edgeSpanwise(iSubdiv+1).top.node2 + nodes_mult*(iRib+1);
		nodeID_edge2D(1) = spanDiv(iRib).edgeSpanwise(iSubdiv+1).bot.node1 + dNodesEE + nodes_mult*iRib;
		nodeID_edge2D(end) = spanDiv(iRib).edgeSpanwise(iSubdiv+1).bot.node2 + nodes_mult*(iRib+1);
		
		spanDiv(iRib).wingSubdiv(iSubdiv).edge1U_nodes = rib(iRib).FE_nodes_3D(dNodesEE + (rib(iRib).FE_nodes_ext_up_sortedCCW(spanDiv(iRib).edgeSpanwise(iSubdiv).top.iNode1:-1:spanDiv(iRib).edgeSpanwise(iSubdiv+1).top.iNode1)), :);
		spanDiv(iRib).wingSubdiv(iSubdiv).edge1D_nodes = rib(iRib).FE_nodes_3D(dNodesEE + (rib(iRib).FE_nodes_ext_dn_sortedCCW(spanDiv(iRib).edgeSpanwise(iSubdiv).bot.iNode1:+1:spanDiv(iRib).edgeSpanwise(iSubdiv+1).bot.iNode1)), :);
		
		nodeID_edge1U = nodes_mult*iRib + rib(iRib).FE_nodes_ext_up_sortedCCW(spanDiv(iRib).edgeSpanwise(iSubdiv).top.iNode1:-1:spanDiv(iRib).edgeSpanwise(iSubdiv+1).top.iNode1);
		nodeID_edge1D = nodes_mult*iRib + rib(iRib).FE_nodes_ext_dn_sortedCCW(spanDiv(iRib).edgeSpanwise(iSubdiv).bot.iNode1:+1:spanDiv(iRib).edgeSpanwise(iSubdiv+1).bot.iNode1);
		
		spanDiv(iRib).wingSubdiv(iSubdiv).edge3U_nodes = rib(iRib+1).FE_nodes_3D(rib(iRib+1).FE_nodes_ext_up_sortedCCW(spanDiv(iRib).edgeSpanwise(iSubdiv).top.iNode2:-1:spanDiv(iRib).edgeSpanwise(iSubdiv+1).top.iNode2), :);
		spanDiv(iRib).wingSubdiv(iSubdiv).edge3D_nodes = rib(iRib+1).FE_nodes_3D(rib(iRib+1).FE_nodes_ext_dn_sortedCCW(spanDiv(iRib).edgeSpanwise(iSubdiv).bot.iNode2:+1:spanDiv(iRib).edgeSpanwise(iSubdiv+1).bot.iNode2), :);
		
		nodeID_edge3U = nodes_mult*(iRib+1) + rib(iRib+1).FE_nodes_ext_up_sortedCCW(spanDiv(iRib).edgeSpanwise(iSubdiv).top.iNode2:-1:spanDiv(iRib).edgeSpanwise(iSubdiv+1).top.iNode2);
		nodeID_edge3D = nodes_mult*(iRib+1) + rib(iRib+1).FE_nodes_ext_dn_sortedCCW(spanDiv(iRib).edgeSpanwise(iSubdiv).bot.iNode2:+1:spanDiv(iRib).edgeSpanwise(iSubdiv+1).bot.iNode2);
		
		% Save nodes for later
		skinPanel_div(iRib).chordDiv(iSubdiv).nodeID_edge1U = nodeID_edge1U;
		skinPanel_div(iRib).chordDiv(iSubdiv).nodeID_edge1D = nodeID_edge1D;
		skinPanel_div(iRib).chordDiv(iSubdiv).nodeID_edge2U = nodeID_edge2U;
		skinPanel_div(iRib).chordDiv(iSubdiv).nodeID_edge2D = nodeID_edge2D;
		skinPanel_div(iRib).chordDiv(iSubdiv).nodeID_edge3U = nodeID_edge3U;
		skinPanel_div(iRib).chordDiv(iSubdiv).nodeID_edge3D = nodeID_edge3D;
		skinPanel_div(iRib).chordDiv(iSubdiv).nodeID_edge4U = nodeID_edge4U;
		skinPanel_div(iRib).chordDiv(iSubdiv).nodeID_edge4D = nodeID_edge4D;
		        
		xVerticesU = [rib(iRib).FE_nodes_3D(dNodesEE + spanDiv(iRib).edgeSpanwise(iSubdiv).top.node1, 1); rib(iRib).FE_nodes_3D(dNodesEE + spanDiv(iRib).edgeSpanwise(iSubdiv+1).top.node1, 1); rib(iRib+1).FE_nodes_3D(spanDiv(iRib).edgeSpanwise(iSubdiv+1).top.node2, 1); rib(iRib+1).FE_nodes_3D(spanDiv(iRib).edgeSpanwise(iSubdiv).top.node2, 1)];
		yVerticesU = [rib(iRib).FE_nodes_3D(dNodesEE + spanDiv(iRib).edgeSpanwise(iSubdiv).top.node1, 3); rib(iRib).FE_nodes_3D(dNodesEE + spanDiv(iRib).edgeSpanwise(iSubdiv+1).top.node1, 3); rib(iRib+1).FE_nodes_3D(spanDiv(iRib).edgeSpanwise(iSubdiv+1).top.node2, 3); rib(iRib+1).FE_nodes_3D(spanDiv(iRib).edgeSpanwise(iSubdiv).top.node2, 3)];
		xVerticesD = [rib(iRib).FE_nodes_3D(dNodesEE + spanDiv(iRib).edgeSpanwise(iSubdiv).bot.node1, 1); rib(iRib).FE_nodes_3D(dNodesEE + spanDiv(iRib).edgeSpanwise(iSubdiv+1).bot.node1, 1); rib(iRib+1).FE_nodes_3D(spanDiv(iRib).edgeSpanwise(iSubdiv+1).bot.node2, 1); rib(iRib+1).FE_nodes_3D(spanDiv(iRib).edgeSpanwise(iSubdiv).bot.node2, 1)];
		yVerticesD = [rib(iRib).FE_nodes_3D(dNodesEE + spanDiv(iRib).edgeSpanwise(iSubdiv).bot.node1, 3); rib(iRib).FE_nodes_3D(dNodesEE + spanDiv(iRib).edgeSpanwise(iSubdiv+1).bot.node1, 3); rib(iRib+1).FE_nodes_3D(spanDiv(iRib).edgeSpanwise(iSubdiv+1).bot.node2, 3); rib(iRib+1).FE_nodes_3D(spanDiv(iRib).edgeSpanwise(iSubdiv).bot.node2, 3)];
		
		edge1U = ((spanDiv(iRib).wingSubdiv(iSubdiv).edge1U_nodes(:,1) - spanDiv(iRib).wingSubdiv(iSubdiv).edge1U_nodes(1,1))./(spanDiv(iRib).wingSubdiv(iSubdiv).edge1U_nodes(end,1) - spanDiv(iRib).wingSubdiv(iSubdiv).edge1U_nodes(1,1)));
		edge1D = ((spanDiv(iRib).wingSubdiv(iSubdiv).edge1D_nodes(:,1) - spanDiv(iRib).wingSubdiv(iSubdiv).edge1D_nodes(1,1))./(spanDiv(iRib).wingSubdiv(iSubdiv).edge1D_nodes(end,1) - spanDiv(iRib).wingSubdiv(iSubdiv).edge1D_nodes(1,1)));
		edge2 = perc_nodes_edge2';
		edge3U = 1-flipud((spanDiv(iRib).wingSubdiv(iSubdiv).edge3U_nodes(:,1) - spanDiv(iRib).wingSubdiv(iSubdiv).edge3U_nodes(1,1))./(spanDiv(iRib).wingSubdiv(iSubdiv).edge3U_nodes(end,1) - spanDiv(iRib).wingSubdiv(iSubdiv).edge3U_nodes(1,1)));
		edge3D = 1-flipud((spanDiv(iRib).wingSubdiv(iSubdiv).edge3D_nodes(:,1) - spanDiv(iRib).wingSubdiv(iSubdiv).edge3D_nodes(1,1))./(spanDiv(iRib).wingSubdiv(iSubdiv).edge3D_nodes(end,1) - spanDiv(iRib).wingSubdiv(iSubdiv).edge3D_nodes(1,1)));
		edge4 = perc_nodes_edge4';
		
		[verticesU, origVerticesU, spanDiv(iRib).wingSubdiv(iSubdiv).sectionsU] = meshQuadrilateral(xVerticesU, yVerticesU, edge1U, edge2, edge3U, edge4);
		[verticesD, origVerticesD, spanDiv(iRib).wingSubdiv(iSubdiv).sectionsD] = meshQuadrilateral(xVerticesD, yVerticesD, edge1D, edge2, edge3D, edge4);
		
		nodeID_verticesU = nodeOffset + (1:size(verticesU, 1));
		spanDiv(iRib).wingSubdiv(iSubdiv).sectionsU = spanDiv(iRib).wingSubdiv(iSubdiv).sectionsU + nodeOffset;
		nodeOffset = nodeOffset + size(verticesU, 1);
		nodeID_verticesD = nodeOffset + (1:size(verticesD, 1));
		spanDiv(iRib).wingSubdiv(iSubdiv).sectionsD = spanDiv(iRib).wingSubdiv(iSubdiv).sectionsD + nodeOffset;
		nodeOffset = nodeOffset + size(verticesD, 1);
		
		nodeID_origVerticesU = [nodeID_edge1U(1:end) + dNodesEE, nodeID_edge2U(2:end-1), fliplr(nodeID_edge3U(1:end)), fliplr(nodeID_edge4U(2:end-1))];
		nodeID_origVerticesD = [nodeID_edge1D(1:end) + dNodesEE, nodeID_edge2D(2:end-1), fliplr(nodeID_edge3D(1:end)), fliplr(nodeID_edge4D(2:end-1))];
		
		% in VERTICES there are the coordinates in the X and Z plane. We have to get Y by means of bi-linear interpolation
		if numel(verticesU) > 0; [xiU, etaU] = QUAD_global2localCS(xVerticesU, yVerticesU, verticesU(:,1), verticesU(:,2)); end
		if numel(verticesD) > 0; [xiD, etaD] = QUAD_global2localCS(xVerticesD, yVerticesD, verticesD(:,1), verticesD(:,2)); end
		% OSS: the transformation distorts the coordinates. ############# FIXED
		[xi_1U, eta_1U] = QUAD_global2localCS(xVerticesU, yVerticesU, spanDiv(iRib).wingSubdiv(iSubdiv).edge1U_nodes(:,1), spanDiv(iRib).wingSubdiv(iSubdiv).edge1U_nodes(:,3));
		[xi_1D, eta_1D] = QUAD_global2localCS(xVerticesD, yVerticesD, spanDiv(iRib).wingSubdiv(iSubdiv).edge1D_nodes(:,1), spanDiv(iRib).wingSubdiv(iSubdiv).edge1D_nodes(:,3));
		[xi_3U, eta_3U] = QUAD_global2localCS(xVerticesU, yVerticesU, spanDiv(iRib).wingSubdiv(iSubdiv).edge3U_nodes(:,1), spanDiv(iRib).wingSubdiv(iSubdiv).edge3U_nodes(:,3));
		[xi_3D, eta_3D] = QUAD_global2localCS(xVerticesD, yVerticesD, spanDiv(iRib).wingSubdiv(iSubdiv).edge3D_nodes(:,1), spanDiv(iRib).wingSubdiv(iSubdiv).edge3D_nodes(:,3));
		
		z_edge1U = interp1(xi_1U, spanDiv(iRib).wingSubdiv(iSubdiv).edge1U_nodes(:,2), xiU);
		z_edge3U = interp1(xi_3U, spanDiv(iRib).wingSubdiv(iSubdiv).edge3U_nodes(:,2), xiU);
		
		z_edge1D = interp1(xi_1D, spanDiv(iRib).wingSubdiv(iSubdiv).edge1D_nodes(:,2), xiD);
		z_edge3D = interp1(xi_3D, spanDiv(iRib).wingSubdiv(iSubdiv).edge3D_nodes(:,2), xiD);
		
		y_verticesU = (1-(etaU./2+.5)).*(z_edge1U) + (etaU./2+.5).*z_edge3U;
		y_verticesD = (1-(etaD./2+.5)).*(z_edge1D) + (etaD./2+.5).*z_edge3D;
		
		if numel(verticesU) > 0
			vertices3DU = [verticesU(:, 1), y_verticesU, verticesU(:, 2)];
		else
			vertices3DU = zeros(0,3);
		end
		
		vertices3DU_all_nodeID = [nodeID_verticesU, nodeID_origVerticesU];
		
		if numel(verticesD) > 0
			vertices3DD = [verticesD(:, 1), y_verticesD, verticesD(:, 2)];
		else
			vertices3DD = zeros(0,3);
		end
		
		vertices3DD_all_nodeID = [nodeID_verticesD, nodeID_origVerticesD];
		
		skinPanels_nodes.Lbound_up(iRib, iSubdiv) = numel(FE_SkinNodeListU)+1;
		skinPanels_nodes.Lbound_dn(iRib, iSubdiv) = numel(FE_SkinNodeListD)+1;
		
		FE_SkinNodeListU = [FE_SkinNodeListU; vertices3DU_all_nodeID'];
		FE_SkinNodeListD = [FE_SkinNodeListD; vertices3DD_all_nodeID'];
		
		skinPanels_nodes.Ubound_up(iRib, iSubdiv) = numel(FE_SkinNodeListU);
		skinPanels_nodes.Ubound_dn(iRib, iSubdiv) = numel(FE_SkinNodeListD);
		
		skinPanels_FE_Grid = [skinPanels_FE_Grid; ...
			[nodeID_verticesU', vertices3DU(:,1), vertices3DU(:,2), vertices3DU(:,3)]; ...
			[nodeID_verticesD', vertices3DD(:,1), vertices3DD(:,2), vertices3DD(:,3)]];
		
		verticesU_all = [verticesU; origVerticesU];
		triU = delaunay(verticesU_all(:,1), verticesU_all(:,2));
		% need normals pointing outwards -> flip first and second row
		triU = [triU(:,2), triU(:,1), triU(:,3)];
		
		verticesD_all = [verticesD; origVerticesD];
		triD = delaunay(verticesD_all(:,1), verticesD_all(:,2));
        
        if iSubdiv == isdv_actuation_1 || iSubdiv == isdv_actuation_2
            rib(iRib).actuation.sdvNode_Up_Pos = [rib(iRib).actuation.sdvNode_Up_Pos, spanDiv(iRib).wingSubdiv(iSubdiv).edge1U_nodes(end,:)'];
            rib(iRib).actuation.sdvNode_Dn_Pos = [rib(iRib).actuation.sdvNode_Dn_Pos, spanDiv(iRib).wingSubdiv(iSubdiv).edge1D_nodes(end,:)'];
            rib(iRib).actuation.sdvNode_Up_Pos = [rib(iRib).actuation.sdvNode_Up_Pos, spanDiv(iRib).wingSubdiv(iSubdiv).edge2U_nodes'];
            rib(iRib).actuation.sdvNode_Dn_Pos = [rib(iRib).actuation.sdvNode_Dn_Pos, spanDiv(iRib).wingSubdiv(iSubdiv).edge2D_nodes'];
            rib(iRib).actuation.sdvNode_Up_Pos = [rib(iRib).actuation.sdvNode_Up_Pos, spanDiv(iRib).wingSubdiv(iSubdiv).edge3U_nodes(end,:)'];
            rib(iRib).actuation.sdvNode_Dn_Pos = [rib(iRib).actuation.sdvNode_Dn_Pos, spanDiv(iRib).wingSubdiv(iSubdiv).edge3D_nodes(end,:)'];
        end
        
		if iSubdiv == isdv_actuation_2
            rib(iRib).actuation.sdvNode_Up_Bk_Pos = [rib(iRib).actuation.sdvNode_Up_Bk_Pos, spanDiv(iRib).wingSubdiv(iSubdiv).edge1U_nodes(end,:)'];
            rib(iRib).actuation.sdvNode_Dn_Bk_Pos = [rib(iRib).actuation.sdvNode_Dn_Bk_Pos, spanDiv(iRib).wingSubdiv(iSubdiv).edge1D_nodes(end,:)'];
            rib(iRib).actuation.sdvNode_Up_Bk_Pos = [rib(iRib).actuation.sdvNode_Up_Bk_Pos, spanDiv(iRib).wingSubdiv(iSubdiv).edge2U_nodes'];
            rib(iRib).actuation.sdvNode_Dn_Bk_Pos = [rib(iRib).actuation.sdvNode_Dn_Bk_Pos, spanDiv(iRib).wingSubdiv(iSubdiv).edge2D_nodes'];
            rib(iRib).actuation.sdvNode_Up_Bk_Pos = [rib(iRib).actuation.sdvNode_Up_Bk_Pos, spanDiv(iRib).wingSubdiv(iSubdiv).edge3U_nodes(end,:)'];
            rib(iRib).actuation.sdvNode_Dn_Bk_Pos = [rib(iRib).actuation.sdvNode_Dn_Bk_Pos, spanDiv(iRib).wingSubdiv(iSubdiv).edge3D_nodes(end,:)'];
        end
		
		% cleanup delaunay-generated triangles
		for i = 1:size(triU, 1)
			tri_angles_U = triangle_angles([verticesU_all(triU(i, 1), :); verticesU_all(triU(i, 2), :); verticesU_all(triU(i, 3), :)]);
			if min(tri_angles_U) < 1e-3
				triU(i, :) = [NaN, NaN, NaN];
			end
		end
		for i = 1:size(triD, 1)
			tri_angles_D = triangle_angles([verticesD_all(triD(i, 1), :); verticesD_all(triD(i, 2), :); verticesD_all(triD(i, 3), :)]);
			if min(tri_angles_D) < 1e-3
				triD(i, :) = [NaN, NaN, NaN];
			end
		end
		
		triU = triU(~isnan(triU));
		triU = reshape(triU, length(triU)/3, 3);
		triD = triD(~isnan(triD));
		triD = reshape(triD, length(triD)/3, 3);
		
		elemU = elemOffset + (1:size(triU,1));
		elemOffset = elemOffset + size(triU,1);
		
		elemD = elemOffset + (1:size(triD,1));
		elemOffset = elemOffset + size(triD,1);
		
        mirrorWingV = vertcat(rib(:).FE_Grid);
        mirrorWing = mirrorWingV(find(abs(mirrorWingV(:,4)+extraRibOffset)<0.00001),:);
		ctria_out = FE_OUT_ctria3(mirrorWing, elemU, PID_Skin_up(iSubdiv, iRib), vertices3DU_all_nodeID(triU(:,1)), vertices3DU_all_nodeID(triU(:,2)), vertices3DU_all_nodeID(triU(:,3)));
        meshOut.ctria = [meshOut.ctria; ctria_out];
        ctria_out = FE_OUT_ctria3(mirrorWing, elemD, PID_Skin_dn(iSubdiv, iRib), vertices3DD_all_nodeID(triD(:,1)), vertices3DD_all_nodeID(triD(:,2)), vertices3DD_all_nodeID(triD(:,3)));
		meshOut.ctria = [meshOut.ctria; ctria_out];
%         ctria_out = FE_OUT_ctria3(mirrorWing, elemU, PID_Skin_up, vertices3DU_all_nodeID(triU(:,1)), vertices3DU_all_nodeID(triU(:,2)), vertices3DU_all_nodeID(triU(:,3)));
%         meshOut.ctria = [meshOut.ctria; ctria_out];
%         ctria_out = FE_OUT_ctria3(mirrorWing, elemD, PID_Skin_dn, vertices3DD_all_nodeID(triD(:,1)), vertices3DD_all_nodeID(triD(:,2)), vertices3DD_all_nodeID(triD(:,3)));
% 		meshOut.ctria = [meshOut.ctria; ctria_out];
	end
end

%% CREATE MAIN SPAR (at sdv location)
% NEED TO DEFINE sdv POSITION FOR EACH RIB, (iSubdiv).POS(iRib)

sparPanels_FE_Grid = [];
FE_nodes_SD_Truss_1 = [];
FE_nodes_SD_Truss_2 = [];
for iRib = 1:numRibs-1
	dNodesEE = rib(iRib).FE_delta_iNodes_EdgeToEdge;
    for isdv_spar = isdv_spar2
		nodeID_edge2 = skinPanel_div(iRib).chordDiv(isdv_spar).nodeID_edge2U;
		nodeID_edge4 = skinPanel_div(iRib).chordDiv(isdv_spar).nodeID_edge2D;
        nodeID_edge1 = fliplr(rib(iRib).spar2.ID + dNodesEE); %nodeID_edge1 = rib(iRib).FE_nodes_spar + dNodesEE;
		nodeID_edge3 = fliplr(rib(iRib+1).spar2.ID); %		nodeID_edge3 = rib(iRib+1).FE_nodes_spar;
		
        xVertices = [rib(iRib).FE_nodes_3D(nodeID_edge1(1), 3); rib(iRib).FE_nodes_3D(nodeID_edge1(end), 3); ...
					 rib(iRib+1).FE_nodes_3D(nodeID_edge3(end), 3); rib(iRib+1).FE_nodes_3D(nodeID_edge3(1), 3)];
		yVertices = [rib(iRib).FE_nodes_3D(nodeID_edge1(1), 2); rib(iRib).FE_nodes_3D(nodeID_edge1(end), 2); ...
					 rib(iRib+1).FE_nodes_3D(nodeID_edge3(end), 2); rib(iRib+1).FE_nodes_3D(nodeID_edge3(1), 2)];

		edge1 = 1-fliplr(rib(iRib).FE_nodes_3D(rib(iRib).spar2.ID + dNodesEE, 2) - yVertices(1))./(yVertices(2) - yVertices(1));
		edge2 = [0; (spanDiv(iRib).wingSubdiv(isdv_spar).edge2U_nodes(:, 3) - xVertices(2))./(xVertices(3) - xVertices(2)); 1];
		edge3 = 1-fliplr(rib(iRib+1).FE_nodes_3D(rib(iRib+1).spar2.ID, 2) - yVertices(4))./(yVertices(3) - yVertices(4));
		edge4 = [0; (spanDiv(iRib).wingSubdiv(isdv_spar).edge2D_nodes(:, 3) - xVertices(1))./(xVertices(4) - xVertices(1)); 1];

		[vertices, origVertices] = meshQuadrilateral(xVertices, yVertices, edge1, edge2, edge3, edge4); % we don't need to take the section node info of the spar

        zVertices = [rib(iRib).FE_nodes_3D(nodeID_edge1(1), 1); rib(iRib).FE_nodes_3D(nodeID_edge1(end), 1); ...
					 rib(iRib+1).FE_nodes_3D(nodeID_edge3(end), 1); rib(iRib+1).FE_nodes_3D(nodeID_edge3(1), 1)];
        zVerticesInt = zVertices(1) + (zVertices(4)-zVertices(1))*edge2;
        xVerticesInt = xVertices(1) + (xVertices(4)-xVertices(1))*edge2;
        yVerticesInt = yVertices(1) + (yVertices(4)-yVertices(1))*edge2;  
        
        FE_nodes_SD_Truss_2 = [FE_nodes_SD_Truss_2; nodeID_edge4', zVerticesInt, yVerticesInt, xVerticesInt];
        
		nodeID_vertices = nodeOffset + (1:size(vertices, 1));
		nodeOffset = nodeOffset + size(vertices, 1);

		nodeID_origVertices = [nodes_mult*iRib + nodeID_edge1(1:end), nodeID_edge2(2:end-1), nodes_mult*(iRib+1) + fliplr(nodeID_edge3(1:end)), fliplr(nodeID_edge4(2:end-1))];
		rib(iRib).actuation.nodes_outline_spar.ID = nodeID_origVertices;
		
        if numel(vertices) > 0
            % Have new vertices been generated?

            % in VERTICES there are the coordinates in the X and Z plane. We have to get Y by means of bi-linear interpolation
            % COUNTER-CLOCKWISE!
            xVerticesA = [xVertices(1); xVertices(4); xVertices(3); xVertices(2)];
            yVerticesA = [yVertices(1); yVertices(4); yVertices(3); yVertices(2)];
            [eta, xi] = QUAD_global2localCS(xVerticesA, yVerticesA, vertices(:,1), vertices(:,2));
            % OSS: the transformation distorts the coordinates. ############# FIXME
            [eta_1, xi_1] = QUAD_global2localCS(xVerticesA, yVerticesA, rib(iRib).FE_nodes_3D(nodeID_edge1, 3), rib(iRib).FE_nodes_3D(nodeID_edge1, 2));
            [eta_3, xi_3] = QUAD_global2localCS(xVerticesA, yVerticesA, rib(iRib+1).FE_nodes_3D(nodeID_edge3, 3), rib(iRib+1).FE_nodes_3D(nodeID_edge3, 2));

            z_edge1 = interp1(xi_1, rib(iRib).FE_nodes_3D(nodeID_edge1, 1), xi);
            z_edge3 = interp1(xi_3, rib(iRib+1).FE_nodes_3D(nodeID_edge3, 1), xi);

            y_vertices = (1-(eta./2+.5)).*(z_edge1) + (eta./2+.5).*z_edge3;

            vertices3D = [y_vertices, vertices(:, 2), vertices(:, 1)];

            vertices3D_all_nodeID = [nodeID_vertices, nodeID_origVertices];

            sparPanels_FE_Grid = [sparPanels_FE_Grid; [nodeID_vertices', vertices3D(:,1), vertices3D(:,2), vertices3D(:,3)]];
        else
            vertices3D_all_nodeID = [nodeID_vertices, nodeID_origVertices];
        end

		vertices_all = [vertices; origVertices];
		tri = delaunay(vertices_all(:,1), vertices_all(:,2));

		% cleanup delaunay-generated triangles
        for i = 1:size(tri, 1)
            tri_angles = triangle_angles([vertices_all(tri(i, 1), :); vertices_all(tri(i, 2), :); vertices_all(tri(i, 3), :)]);
            if min(tri_angles) < 1e-3
                tri(i, :) = [NaN, NaN, NaN];
            end
        end
        
		tri = tri(~isnan(tri));
		tri = reshape(tri, length(tri)/3, 3);

		elem = elemOffset + (1:size(tri,1));
		elemOffset = elemOffset + size(tri,1);
        
%         PID_Spar_use = PID_Spar(iRib);
        PID_Spar_use = PID_Spar;
        
        mirrorWingV = vertcat(rib(:).FE_Grid);
        mirrorWing = mirrorWingV(find(abs(mirrorWingV(:,4)+extraRibOffset)<0.00001),:);
		ctria_out = FE_OUT_ctria3(mirrorWing, elem, PID_Spar_use, vertices3D_all_nodeID(tri(:,1)), vertices3D_all_nodeID(tri(:,2)), vertices3D_all_nodeID(tri(:,3)));
		meshOut.ctria = [meshOut.ctria; ctria_out];
    end

    for isdv_spar = isdv_spars
        nodeID_edge2 = skinPanel_div(iRib).chordDiv(isdv_spar).nodeID_edge2U;
        nodeID_edge4 = skinPanel_div(iRib).chordDiv(isdv_spar).nodeID_edge2D;
        nodeID_edge1 = rib(iRib).FE_nodes_spar + dNodesEE;
        nodeID_edge3 = rib(iRib+1).FE_nodes_spar;

        xVertices = [rib(iRib).FE_nodes_3D(nodeID_edge1(1), 3); rib(iRib).FE_nodes_3D(nodeID_edge1(end), 3); ...
                     rib(iRib+1).FE_nodes_3D(nodeID_edge3(end), 3); rib(iRib+1).FE_nodes_3D(nodeID_edge3(1), 3)];
        yVertices = [rib(iRib).FE_nodes_3D(nodeID_edge1(1), 2); rib(iRib).FE_nodes_3D(nodeID_edge1(end), 2); ...
                     rib(iRib+1).FE_nodes_3D(nodeID_edge3(end), 2); rib(iRib+1).FE_nodes_3D(nodeID_edge3(1), 2)];

        edge1 = (rib(iRib).FE_nodes_3D(dNodesEE + rib(iRib).FE_nodes_spar, 2) - yVertices(1))./(yVertices(2) - yVertices(1));
        edge2 = [0; (spanDiv(iRib).wingSubdiv(isdv_spar).edge2U_nodes(:, 3) - xVertices(2))./(xVertices(3) - xVertices(2)); 1];
        edge3 = 1-flipud(rib(iRib+1).FE_nodes_3D(rib(iRib+1).FE_nodes_spar, 2) - yVertices(4))./(yVertices(3) - yVertices(4));
        edge4 = [0; (spanDiv(iRib).wingSubdiv(isdv_spar).edge2D_nodes(:, 3) - xVertices(1))./(xVertices(4) - xVertices(1)); 1];

        [vertices, origVertices] = meshQuadrilateral(xVertices, yVertices, edge1, edge2, edge3, edge4); % we don't need to take the section node info of the spar

        zVertices = [rib(iRib).FE_nodes_3D(nodeID_edge1(1), 1); rib(iRib).FE_nodes_3D(nodeID_edge1(end), 1); ...
                     rib(iRib+1).FE_nodes_3D(nodeID_edge3(end), 1); rib(iRib+1).FE_nodes_3D(nodeID_edge3(1), 1)];
        zVerticesInt = zVertices(1) + (zVertices(4)-zVertices(1))*edge2;
        xVerticesInt = xVertices(1) + (xVertices(4)-xVertices(1))*edge2;
        yVerticesInt = yVertices(1) + (yVertices(4)-yVertices(1))*edge2;  

        FE_nodes_SD_Truss_1 = [FE_nodes_SD_Truss_1; nodeID_edge4', zVerticesInt, yVerticesInt, xVerticesInt];

        nodeID_vertices = nodeOffset + (1:size(vertices, 1));
        nodeOffset = nodeOffset + size(vertices, 1);

        nodeID_origVertices = [nodes_mult*iRib + nodeID_edge1(1:end), nodeID_edge2(2:end-1), nodes_mult*(iRib+1) + fliplr(nodeID_edge3(1:end)), fliplr(nodeID_edge4(2:end-1))];
        rib(iRib).actuation.nodes_outline_spar.ID = nodeID_origVertices;

        if numel(vertices) > 0
            % Have new vertices been generated?

            % in VERTICES there are the coordinates in the X and Z plane. We have to get Y by means of bi-linear interpolation
            % COUNTER-CLOCKWISE!
            xVerticesA = [xVertices(1); xVertices(4); xVertices(3); xVertices(2)];
            yVerticesA = [yVertices(1); yVertices(4); yVertices(3); yVertices(2)];
            [eta, xi] = QUAD_global2localCS(xVerticesA, yVerticesA, vertices(:,1), vertices(:,2));
            % OSS: the transformation distorts the coordinates. ############# FIXME
            [eta_1, xi_1] = QUAD_global2localCS(xVerticesA, yVerticesA, rib(iRib).FE_nodes_3D(nodeID_edge1, 3), rib(iRib).FE_nodes_3D(nodeID_edge1, 2));
            [eta_3, xi_3] = QUAD_global2localCS(xVerticesA, yVerticesA, rib(iRib+1).FE_nodes_3D(nodeID_edge3, 3), rib(iRib+1).FE_nodes_3D(nodeID_edge3, 2));

            z_edge1 = interp1(xi_1, rib(iRib).FE_nodes_3D(nodeID_edge1, 1), xi);
            z_edge3 = interp1(xi_3, rib(iRib+1).FE_nodes_3D(nodeID_edge3, 1), xi);

            y_vertices = (1-(eta./2+.5)).*(z_edge1) + (eta./2+.5).*z_edge3;

            vertices3D = [y_vertices, vertices(:, 2), vertices(:, 1)];

            vertices3D_all_nodeID = [nodeID_vertices, nodeID_origVertices];

            sparPanels_FE_Grid = [sparPanels_FE_Grid; [nodeID_vertices', vertices3D(:,1), vertices3D(:,2), vertices3D(:,3)]];
        else
            vertices3D_all_nodeID = [nodeID_vertices, nodeID_origVertices];
        end

        vertices_all = [vertices; origVertices];
        tri = delaunay(vertices_all(:,1), vertices_all(:,2));

        % cleanup delaunay-generated triangles
        for i = 1:size(tri, 1)
            tri_angles = triangle_angles([vertices_all(tri(i, 1), :); vertices_all(tri(i, 2), :); vertices_all(tri(i, 3), :)]);
            if min(tri_angles) < 1e-3
                tri(i, :) = [NaN, NaN, NaN];
            end
        end

        tri = tri(~isnan(tri));
        tri = reshape(tri, length(tri)/3, 3);

        elem = elemOffset + (1:size(tri,1));
        elemOffset = elemOffset + size(tri,1);

%         PID_Spar_use = PID_Spar(iRib);
        PID_Spar_use = PID_Spar;

        mirrorWingV = vertcat(rib(:).FE_Grid);
        mirrorWing = mirrorWingV(find(abs(mirrorWingV(:,4)+extraRibOffset)<0.00001),:);
        ctria_out = FE_OUT_ctria3(mirrorWing, elem, PID_Spar_use, vertices3D_all_nodeID(tri(:,1)), vertices3D_all_nodeID(tri(:,2)), vertices3D_all_nodeID(tri(:,3)));
        meshOut.ctria = [meshOut.ctria; ctria_out];

        if mod(iRib,2) == 0
            rib(iRib).actuation.sparNode_ID = vertices3D_all_nodeID;
            rib(iRib).actuation.sparNode_Pos = vertices3D';
            rib(iRib).actuation.sparNode_Pos = [rib(iRib).actuation.sparNode_Pos, rib(iRib).FE_nodes_3D(dNodesEE + rib(iRib).FE_nodes_spar, :)'];
            rib(iRib).actuation.sparNode_Pos = [rib(iRib).actuation.sparNode_Pos, spanDiv(iRib).wingSubdiv(isdv_spar).edge2U_nodes'];
            rib(iRib).actuation.sparNode_Pos = [rib(iRib).actuation.sparNode_Pos, flipud(rib(iRib+1).FE_nodes_3D(rib(iRib+1).FE_nodes_spar, :))'];
            rib(iRib).actuation.sparNode_Pos = [rib(iRib).actuation.sparNode_Pos, flipud(spanDiv(iRib).wingSubdiv(isdv_spar).edge2D_nodes)'];
            rib(iRib).actuation.nodes_outline_spar.ID = nodeID_origVertices;

            rib(iRib).actuation.sparNode_Up_ID = skinPanel_div(iRib).chordDiv(isdv_spar).nodeID_edge2U;
            rib(iRib).actuation.sparNode_Up_Pos = zeros(3,length(rib(iRib).actuation.sparNode_Up_ID));
            rib(iRib).actuation.sparNode_Up_Pos(:,1) = rib(iRib).FE_nodes_3D(nodeID_edge1(end), :)';
            rib(iRib).actuation.sparNode_Up_Pos(:,end) = rib(iRib+1).FE_nodes_3D(nodeID_edge3(end), :)';
            rib(iRib).actuation.sparNode_Up_Pos(:,2:end-1) = spanDiv(iRib).wingSubdiv(3).edge2U_nodes';

            rib(iRib).actuation.sparNode_Dn_ID = skinPanel_div(iRib).chordDiv(isdv_spar).nodeID_edge2D;
            rib(iRib).actuation.sparNode_Dn_Pos = zeros(3,length(rib(iRib).actuation.sparNode_Dn_ID));
            rib(iRib).actuation.sparNode_Dn_Pos(:,1) = rib(iRib).FE_nodes_3D(nodeID_edge1(1), :)';
            rib(iRib).actuation.sparNode_Dn_Pos(:,end) = rib(iRib+1).FE_nodes_3D(nodeID_edge3(1), :)';
            rib(iRib).actuation.sparNode_Dn_Pos(:,2:end-1) = spanDiv(iRib).wingSubdiv(3).edge2D_nodes';

            % define truss attachment points at the spar
            rib(iRib).truss.sparNode_Dn_Pos =  rib(iRib).actuation.sparNode_Dn_Pos;
            rib(iRib).truss.sparNode_Dn_ID =  rib(iRib).actuation.sparNode_Dn_ID; 
        end
    end
end


%% Create actuation
actuation_FE_Grid = [];
actuationID = 9;

for iRib = numRibs-2:-2:numRibs-2*wingDesign.nRibsC
    
    % sdv ID
    rib(iRib).actuation.sdvNode_Up_ID = skinPanel_div(iRib).chordDiv(isdv_actuation_1).nodeID_edge2U;
    rib(iRib).actuation.sdvNode_Dn_ID = skinPanel_div(iRib).chordDiv(isdv_actuation_1).nodeID_edge2D;
    rib(iRib).actuation.sdvNode_Up_Bk_ID = skinPanel_div(iRib).chordDiv(isdv_actuation_2).nodeID_edge2U;
    rib(iRib).actuation.sdvNode_Dn_Bk_ID = skinPanel_div(iRib).chordDiv(isdv_actuation_2).nodeID_edge2D;
    
    rib(iRib).actuation.location_spar_ID = elemOffset + 1;
    rib(iRib).actuation.location_sdv_ID = elemOffset + 2;
    rib(iRib).actuation.spar_rbe_ID =elemOffset + 3;
    rib(iRib).actuation.sdv_rbe_ID =elemOffset + 4;
    rib(iRib).actuation.rod_ID = elemOffset + 5;
    elem = elemOffset + 5;
    
    rib(iRib).actuation.location_spar_Pos_prov = [rib(iRib).actuation.sparNode_Up_Pos(:,1), rib(iRib).actuation.sparNode_Up_Pos(:,end), rib(iRib).actuation.sparNode_Dn_Pos(:,1), rib(iRib).actuation.sparNode_Dn_Pos(:,end)];
    rib(iRib).actuation.location_normal = (cross((rib(iRib).actuation.sparNode_Up_Pos(:,end)-rib(iRib).actuation.sparNode_Up_Pos(:,1)),[0;-1;0]));
    rib(iRib).actuation.location_spar_Pos = mean(rib(iRib).actuation.location_spar_Pos_prov')' + actuation_spar_dist*rib(iRib).actuation.location_normal/norm(rib(iRib).actuation.location_normal);
    rib(iRib).actuation.location_spar_Pos(2) = mean(rib(iRib).actuation.sparNode_Dn_Pos(2,:))+rib(iRib).actuation.param_spar_y*(mean(rib(iRib).actuation.sparNode_Up_Pos(2,:))-mean(rib(iRib).actuation.sparNode_Dn_Pos(2,:)));
	
	rib(iRib).actuation.location_sdv_Pos = mean(rib(iRib).actuation.sdvNode_Dn_Pos');
    rib(iRib).actuation.location_sdv_Pos(2) = mean(rib(iRib).actuation.sdvNode_Dn_Pos(2,:)) - rib(iRib).actuation.param_connector_y * (mean(rib(iRib).actuation.sdvNode_Dn_Pos(2,:)) - mean(rib(iRib).actuation.sdvNode_Up_Pos(2,:)));
    distance_rod_points = [rib(iRib).actuation.location_spar_Pos'; rib(iRib).actuation.location_sdv_Pos];
	rib(iRib).actuation.rodlength = pdist(distance_rod_points);
    
    mirrorWingV = vertcat(rib(:).FE_Grid);
    mirrorWing = mirrorWingV(find(abs(mirrorWingV(:,4)+extraRibOffset)<0.00001),:);  
    for pos = 1:length(rib(iRib).actuation.nodes_outline_spar.ID)-1
        element = elem + pos;
        ctria_out = FE_OUT_ctria3(mirrorWing, element, actuationID, rib(iRib).actuation.location_spar_ID, rib(iRib).actuation.nodes_outline_spar.ID(pos), rib(iRib).actuation.nodes_outline_spar.ID(pos+1));
        meshOut.ctria = [meshOut.ctria; ctria_out];
    end
    
    elem = element+1;
    
    ctria_out = FE_OUT_ctria3(mirrorWing, elem, actuationID, rib(iRib).actuation.location_spar_ID, rib(iRib).actuation.nodes_outline_spar.ID(end), rib(iRib).actuation.nodes_outline_spar.ID(1));
    meshOut.ctria = [meshOut.ctria; ctria_out];
    
    elem = element+2;
       
    for pos2 = 1:length(rib(iRib).actuation.sdvNode_Dn_ID)-1
        element2 = elem + pos2;
        ctria_out = FE_OUT_ctria3(mirrorWing, element2, actuationID, rib(iRib).actuation.location_sdv_ID, rib(iRib).actuation.sdvNode_Dn_ID(pos2), rib(iRib).actuation.sdvNode_Dn_ID(pos2+1));
        meshOut.ctria = [meshOut.ctria; ctria_out];
    end
    
    elem = element2;
    
     for pos3 = 1:length(rib(iRib).actuation.sdvNode_Dn_Bk_ID)-1
        element3 = elem + pos3;
        ctria_out = FE_OUT_ctria3(mirrorWing, element3, actuationID, rib(iRib).actuation.location_sdv_ID, rib(iRib).actuation.sdvNode_Dn_Bk_ID(pos3), rib(iRib).actuation.sdvNode_Dn_Bk_ID(pos3+1));
        meshOut.ctria = [meshOut.ctria; ctria_out];
     end   
    
    elem = element3;
            
    actuation_FE_Grid = [actuation_FE_Grid; ...
        [rib(iRib).actuation.location_spar_ID, rib(iRib).actuation.location_spar_Pos(1), rib(iRib).actuation.location_spar_Pos(2), rib(iRib).actuation.location_spar_Pos(3)];...
        [rib(iRib).actuation.location_sdv_ID, rib(iRib).actuation.location_sdv_Pos(1), rib(iRib).actuation.location_sdv_Pos(2), rib(iRib).actuation.location_sdv_Pos(3)]];
    elemOffset = elem;
end


%% Create SPCs 
iRibSPC = 2;
[~, iRear_wingbox_dn] = find(rib(iRibSPC).FE_nodes_ext_dn_sortedCCW == rib(iRibSPC).corner_nodes(5));
[~, iRear_wingbox_up] = find(rib(iRibSPC).FE_nodes_ext_up_sortedCCW == rib(iRibSPC).corner_nodes(4));

SPC_node_list = [rib(iRibSPC).FE_nodes_ext_up_sortedCCW((iRear_wingbox_up+1):end)+iRibSPC*nodes_mult*1, rib(iRibSPC).FE_nodes_ext_dn_sortedCCW(1:(iRear_wingbox_dn-1))+iRibSPC*nodes_mult*1, rib(iRibSPC).FE_nodes_spar+iRibSPC*nodes_mult*1];    

iRibSPC = iRibSPC + 1;
[~, iRear_wingbox_dn] = find(rib(iRibSPC).FE_nodes_ext_dn_sortedCCW == rib(iRibSPC).corner_nodes(5));
[~, iRear_wingbox_up] = find(rib(iRibSPC).FE_nodes_ext_up_sortedCCW == rib(iRibSPC).corner_nodes(4));
SPC_node_list = [SPC_node_list, rib(iRibSPC).FE_nodes_ext_up_sortedCCW((iRear_wingbox_up+1):end)+iRibSPC*nodes_mult*1, rib(iRibSPC).FE_nodes_ext_dn_sortedCCW(1:(iRear_wingbox_dn-1))+iRibSPC*nodes_mult*1, rib(iRibSPC).FE_nodes_spar+iRibSPC*nodes_mult*1];    


%% SKIN NODE LIST
FE_SkinNodeList_ribsU = [];
FE_SkinNodeList_ribsD = [];
for iRib = 1:numRibs
	dNodesEE = rib(iRib).FE_delta_iNodes_EdgeToEdge;
	if rib_numElemThickness(iRib) > 0
		FE_SkinNodeList_UP = nodes_mult*iRib + repmat(rib(iRib).FE_nodes_ext_up_sortedCCW', 1, rib_numElemThickness(iRib)+1) + repmat((0:rib_numElemThickness(iRib))*(dNodesEE/rib_numElemThickness(iRib)), length(rib(iRib).FE_nodes_ext_up_sortedCCW),1);
		FE_SkinNodeList_DN = nodes_mult*iRib + repmat(rib(iRib).FE_nodes_ext_dn_sortedCCW(1:end)', 1, rib_numElemThickness(iRib)+1) + repmat((0:rib_numElemThickness(iRib))*(dNodesEE/rib_numElemThickness(iRib)), length(rib(iRib).FE_nodes_ext_dn_sortedCCW),1);
	else
		FE_SkinNodeList_UP = [];
		FE_SkinNodeList_DN = [];
	end
	FE_SkinNodeList_ribsU = [FE_SkinNodeList_ribsU; FE_SkinNodeList_UP(:)];
	FE_SkinNodeList_ribsD = [FE_SkinNodeList_ribsD; FE_SkinNodeList_DN(:)];
end
FE_SkinNodeListU = [FE_SkinNodeListU; FE_SkinNodeList_ribsU];
FE_SkinNodeListD = [FE_SkinNodeListD; FE_SkinNodeList_ribsD];


%% SKIN NODES
FE_Grid_mat = [vertcat(rib(:).FE_Grid); skinPanels_FE_Grid; sparPanels_FE_Grid; actuation_FE_Grid];
FE_Grid = struct('ID', FE_Grid_mat(:,1), 'X', FE_Grid_mat(:,2), 'Y', FE_Grid_mat(:,3), 'Z', FE_Grid_mat(:,4));

%% define rigid ribs

% set if rigid ribs are added where no morphing structure is defined
addRearPlate = 1;
addFrontPlate = 1;
addMidPlate = 1;

rigidRibPID = 5;

for i = 2:numRibs-1
    
    % inner rib side sandwich
    coordpointsUP = rib(i).FE_nodes_3D(rib(i).FE_nodes_ext_up_sortedCCW,:);
    coordpointsDN = rib(i).FE_nodes_3D(rib(i).FE_nodes_ext_dn_sortedCCW,:);
    coordpointsSPARrear = rib(i).FE_nodes_3D(rib(i).FE_nodes_spar,:);
    coordpointsSPARfront = rib(i).FE_nodes_3D(rib(i).spar2.ID,:);
    sparRearPosPLATE = coordpointsSPARrear(1,1);
    sparFrontPosPLATE = coordpointsSPARfront(1,1);
    
    coordpointsFRONT = [coordpointsUP(coordpointsUP(:,1)<=sparFrontPosPLATE,:); coordpointsDN(coordpointsDN(:,1)<=sparFrontPosPLATE,:); coordpointsSPARfront];
    coordpointsREAR = [coordpointsUP(coordpointsUP(:,1)>sparRearPosPLATE,:);  coordpointsSPARrear(end:-1:1,:); coordpointsDN(coordpointsDN(:,1)>sparRearPosPLATE,:)];
    coordpointsMID = [coordpointsUP((coordpointsUP(:,1)>=sparFrontPosPLATE & coordpointsUP(:,1)<=sparRearPosPLATE),:);...
        coordpointsSPARfront(2:end-1,:);...
        coordpointsDN((coordpointsDN(:,1)>=sparFrontPosPLATE & coordpointsDN(:,1)<=sparRearPosPLATE),:);...
        coordpointsSPARrear(2:end-1,:)];
    
    trailing_edge_Z =  coordpointsUP(1,3);
    spar_Rear_Z = coordpointsSPARrear(1,3);
    spar_Front_Z = coordpointsSPARfront(1,3);
    leading_edge_Z = coordpointsUP(end,3);
    [p2Dplate,t2Dplate]=distmesh2d(@dpoly,@huniform,0.02,[-1,-1; 1,1],coordpointsREAR(:,1:2),coordpointsREAR(:,1:2));
    [p2Dplatefront,t2Dplatefront]=distmesh2d(@dpoly,@huniform,0.02,[-1,-1; 1,1],coordpointsFRONT(:,1:2),coordpointsFRONT(:,1:2));
    [p2DplateMID,t2DplateMID]=distmesh2d(@dpoly,@huniform,0.02,[-1,-1; 1,1],coordpointsMID(:,1:2),coordpointsMID(:,1:2));
    
    deltaIDplate = 200000+5020*i;
    t2Dplate = t2Dplate+deltaIDplate;
    FE_Grid_innerPlate.ID = [1:length(p2Dplate)]' + deltaIDplate;
    FE_Grid_innerPlate.X = p2Dplate(:,1);
    FE_Grid_innerPlate.Y = p2Dplate(:,2);
    FE_Grid_innerPlate.Z = interp1([coordpointsUP(1,1),sparRearPosPLATE],[trailing_edge_Z, spar_Rear_Z], FE_Grid_innerPlate.X);%ones(length(p2Dplate),1)*coordpointsUP(1,3);
    tolPlate = 0.00001;
    tipplateDoubleA = [];
    
    for j = 1:length(FE_Grid_innerPlate.X)
        tipplateDouble = ~isempty(FE_Grid.ID(abs(FE_Grid.X - FE_Grid_innerPlate.X(j))+abs(FE_Grid.Y - FE_Grid_innerPlate.Y(j))+abs(FE_Grid.Z - FE_Grid_innerPlate.Z(j)) < tolPlate));
        if tipplateDouble
            t2Dplate(t2Dplate==FE_Grid_innerPlate.ID(j)) = FE_Grid.ID(abs(FE_Grid.X - FE_Grid_innerPlate.X(j))+abs(FE_Grid.Y - FE_Grid_innerPlate.Y(j))+abs(FE_Grid.Z - FE_Grid_innerPlate.Z(j))<tolPlate);
            tipplateDoubleA = [tipplateDoubleA;j];
        end
    end
    FE_Grid_innerPlate.ID(tipplateDoubleA) = [];
    FE_Grid_innerPlate.X(tipplateDoubleA) = [];
    FE_Grid_innerPlate.Y(tipplateDoubleA) = [];
    FE_Grid_innerPlate.Z(tipplateDoubleA) = [];
    
    %%%
    if addRearPlate && i <= (2*(wingDesign.nRibs-wingDesign.nRibsC)+1)%11
        FE_Grid.ID = [FE_Grid.ID; FE_Grid_innerPlate.ID];
        FE_Grid.X = [FE_Grid.X; FE_Grid_innerPlate.X];
        FE_Grid.Y = [FE_Grid.Y; FE_Grid_innerPlate.Y];
        FE_Grid.Z = [FE_Grid.Z; FE_Grid_innerPlate.Z];
               
        IDroot_rearEL = [1:length(t2Dplate)]' + deltaIDplate;
        ctria_out = FE_OUT_ctria3([], IDroot_rearEL, rigidRibPID, t2Dplate(:,1), t2Dplate(:,2), t2Dplate(:,3));
        meshOut.ctria = [meshOut.ctria; ctria_out];
    end
    
    %%%
    if addFrontPlate 
        deltaIDplate = deltaIDplate+7600*i;
        t2Dplatefront = t2Dplatefront+deltaIDplate;
        FE_Grid_innerPlatefront.ID = [1:length(p2Dplatefront)]' + deltaIDplate;
        FE_Grid_innerPlatefront.X = p2Dplatefront(:,1);
        FE_Grid_innerPlatefront.Y = p2Dplatefront(:,2);
        FE_Grid_innerPlatefront.Z = interp1([coordpointsUP(end,1),sparFrontPosPLATE],[leading_edge_Z, spar_Front_Z], FE_Grid_innerPlatefront.X);%ones(length(p2Dplate),1)*coordpointsUP(1,3);
        tolPlate = 0.00001;
        tipplateDoubleA = [];
        
        for j = 1:length(FE_Grid_innerPlatefront.X)
            tipplateDouble = ~isempty(FE_Grid.ID(abs(FE_Grid.X - FE_Grid_innerPlatefront.X(j))+abs(FE_Grid.Y - FE_Grid_innerPlatefront.Y(j))+abs(FE_Grid.Z - FE_Grid_innerPlatefront.Z(j)) < tolPlate));
            if tipplateDouble
                t2Dplatefront(t2Dplatefront==FE_Grid_innerPlatefront.ID(j)) = FE_Grid.ID(abs(FE_Grid.X - FE_Grid_innerPlatefront.X(j))+abs(FE_Grid.Y - FE_Grid_innerPlatefront.Y(j))+abs(FE_Grid.Z - FE_Grid_innerPlatefront.Z(j))<tolPlate);
                tipplateDoubleA = [tipplateDoubleA;j];
            end
        end
        FE_Grid_innerPlatefront.ID(tipplateDoubleA) = [];
        FE_Grid_innerPlatefront.X(tipplateDoubleA) = [];
        FE_Grid_innerPlatefront.Y(tipplateDoubleA) = [];
        FE_Grid_innerPlatefront.Z(tipplateDoubleA) = [];
        
        FE_Grid.ID = [FE_Grid.ID; FE_Grid_innerPlatefront.ID];
        FE_Grid.X = [FE_Grid.X; FE_Grid_innerPlatefront.X];
        FE_Grid.Y = [FE_Grid.Y; FE_Grid_innerPlatefront.Y];
        FE_Grid.Z = [FE_Grid.Z; FE_Grid_innerPlatefront.Z];
        
        ID_frontEL = [1:length(t2Dplatefront)]' + deltaIDplate;
        ctria_out = FE_OUT_ctria3([], ID_frontEL, rigidRibPID, t2Dplatefront(:,1), t2Dplatefront(:,2), t2Dplatefront(:,3));
        meshOut.ctria = [meshOut.ctria; ctria_out];
    end
    
    %%%
    if addMidPlate
        
        deltaIDplate = deltaIDplate+7800*i;
        t2DplateMID = t2DplateMID+deltaIDplate;
        FE_Grid_innerPlatefront.ID = [1:length(p2DplateMID)]' + deltaIDplate;
        FE_Grid_innerPlatefront.X = p2DplateMID(:,1);
        FE_Grid_innerPlatefront.Y = p2DplateMID(:,2);
        FE_Grid_innerPlatefront.Z = ones(length(p2DplateMID),1)*coordpointsUP(1,3);
        tolPlate = 0.00001;
        tipplateDoubleA = [];
        
        for j = 1:length(FE_Grid_innerPlatefront.X)
            tipplateDouble = ~isempty(FE_Grid.ID(abs(FE_Grid.X - FE_Grid_innerPlatefront.X(j))+abs(FE_Grid.Y - FE_Grid_innerPlatefront.Y(j))+abs(FE_Grid.Z - FE_Grid_innerPlatefront.Z(j)) < tolPlate));
            if tipplateDouble
                t2DplateMID(t2DplateMID==FE_Grid_innerPlatefront.ID(j)) = FE_Grid.ID(abs(FE_Grid.X - FE_Grid_innerPlatefront.X(j))+abs(FE_Grid.Y - FE_Grid_innerPlatefront.Y(j))+abs(FE_Grid.Z - FE_Grid_innerPlatefront.Z(j))<tolPlate);
                tipplateDoubleA = [tipplateDoubleA;j];
            end
        end
        FE_Grid_innerPlatefront.ID(tipplateDoubleA) = [];
        FE_Grid_innerPlatefront.X(tipplateDoubleA) = [];
        FE_Grid_innerPlatefront.Y(tipplateDoubleA) = [];
        FE_Grid_innerPlatefront.Z(tipplateDoubleA) = [];
        
        FE_Grid.ID = [FE_Grid.ID; FE_Grid_innerPlatefront.ID];
        FE_Grid.X = [FE_Grid.X; FE_Grid_innerPlatefront.X];
        FE_Grid.Y = [FE_Grid.Y; FE_Grid_innerPlatefront.Y];
        FE_Grid.Z = [FE_Grid.Z; FE_Grid_innerPlatefront.Z];

        ID_frontEL = [1:size(t2DplateMID,1)]' + deltaIDplate;
        ctria_out = FE_OUT_ctria3([], ID_frontEL, rigidRibPID, t2DplateMID(:,1), t2DplateMID(:,2), t2DplateMID(:,3));
        meshOut.ctria = [meshOut.ctria; ctria_out];
    end
    
    
    % outer rib side sandwich
    coordpointsUP = rib(i).FE_nodes_3D(rib(i).FE_nodes_ext_up_sortedCCW+rib(i).FE_delta_iNodes_EdgeToEdge,:);
    coordpointsDN = rib(i).FE_nodes_3D(rib(i).FE_nodes_ext_dn_sortedCCW+rib(i).FE_delta_iNodes_EdgeToEdge,:);
    coordpointsSPARrear = rib(i).FE_nodes_3D(rib(i).FE_nodes_spar+rib(i).FE_delta_iNodes_EdgeToEdge,:);
    coordpointsSPARfront = rib(i).FE_nodes_3D(rib(i).spar2.ID+rib(i).FE_delta_iNodes_EdgeToEdge,:);
    sparRearPosPLATE = coordpointsSPARrear(1,1);
    sparFrontPosPLATE = coordpointsSPARfront(1,1);
    
    coordpointsREAR = [coordpointsUP(coordpointsUP(:,1)>sparRearPosPLATE,:);  coordpointsSPARrear(end:-1:1,:); coordpointsDN(coordpointsDN(:,1)>sparRearPosPLATE,:)];
    coordpointsFRONT = [coordpointsUP(coordpointsUP(:,1)<=sparFrontPosPLATE,:); coordpointsDN(coordpointsDN(:,1)<=sparFrontPosPLATE,:); coordpointsSPARfront];
    coordpointsMID = [coordpointsUP((coordpointsUP(:,1)>=sparFrontPosPLATE & coordpointsUP(:,1)<=sparRearPosPLATE),:);...
        coordpointsSPARfront(2:end-1,:);...
        coordpointsDN((coordpointsDN(:,1)>=sparFrontPosPLATE & coordpointsDN(:,1)<=sparRearPosPLATE),:);...
        coordpointsSPARrear(2:end-1,:)];
    
    trailing_edge_Z =  coordpointsUP(1,3);
    spar_Rear_Z = coordpointsSPARrear(1,3);
    spar_Front_Z = coordpointsSPARfront(1,3);
    leading_edge_Z = coordpointsUP(end,3);
    
    [p2Dplate,t2Dplate]=distmesh2d(@dpoly,@huniform,0.02,[-1,-1; 1,1],coordpointsREAR(:,1:2),coordpointsREAR(:,1:2));
    [p2Dplatefront,t2Dplatefront]=distmesh2d(@dpoly,@huniform,0.02,[-1,-1; 1,1],coordpointsFRONT(:,1:2),coordpointsFRONT(:,1:2));
    [p2DplateMID,t2DplateMID]=distmesh2d(@dpoly,@huniform,0.02,[-1,-1; 1,1],coordpointsMID(:,1:2),coordpointsMID(:,1:2));

    deltaIDplate = deltaIDplate+3300*i;
    t2Dplate = t2Dplate+deltaIDplate;
    FE_Grid_outerPlate.ID = [1:length(p2Dplate)]' + deltaIDplate;
    FE_Grid_outerPlate.X = p2Dplate(:,1);
    FE_Grid_outerPlate.Y = p2Dplate(:,2);
    FE_Grid_outerPlate.Z = interp1([coordpointsUP(1,1),sparRearPosPLATE],[trailing_edge_Z, spar_Rear_Z], FE_Grid_outerPlate.X);%ones(length(p2Dplate),1)*coordpointsUP(1,3);
    tolPlate = 0.00001;
    tipplateDoubleA = [];
    
    for j = 1:length(FE_Grid_outerPlate.X)
        tipplateDouble = ~isempty(FE_Grid.ID(abs(FE_Grid.X - FE_Grid_outerPlate.X(j))+abs(FE_Grid.Y - FE_Grid_outerPlate.Y(j))+abs(FE_Grid.Z - FE_Grid_outerPlate.Z(j)) < tolPlate));
        if tipplateDouble
            t2Dplate(t2Dplate==FE_Grid_outerPlate.ID(j)) = FE_Grid.ID(abs(FE_Grid.X - FE_Grid_outerPlate.X(j))+abs(FE_Grid.Y - FE_Grid_outerPlate.Y(j))+abs(FE_Grid.Z - FE_Grid_outerPlate.Z(j))<tolPlate);
            tipplateDoubleA = [tipplateDoubleA;j];
        end
    end
    FE_Grid_outerPlate.ID(tipplateDoubleA) = [];
    FE_Grid_outerPlate.X(tipplateDoubleA) = [];
    FE_Grid_outerPlate.Y(tipplateDoubleA) = [];
    FE_Grid_outerPlate.Z(tipplateDoubleA) = [];
    
    %%%
    if addRearPlate && i<= (2*(wingDesign.nRibs-wingDesign.nRibsC)+1)%11
        FE_Grid.ID = [FE_Grid.ID; FE_Grid_outerPlate.ID];
        FE_Grid.X = [FE_Grid.X; FE_Grid_outerPlate.X];
        FE_Grid.Y = [FE_Grid.Y; FE_Grid_outerPlate.Y];
        FE_Grid.Z = [FE_Grid.Z; FE_Grid_outerPlate.Z];
                
        IDroot_rearEL = [1:length(t2Dplate)]' + deltaIDplate;
        ctria_out = FE_OUT_ctria3([], IDroot_rearEL, rigidRibPID, t2Dplate(:,1), t2Dplate(:,2), t2Dplate(:,3));
        meshOut.ctria = [meshOut.ctria; ctria_out];
    end
    
    %%%
    if addFrontPlate
        deltaIDplate = deltaIDplate+4230*i;
        t2Dplatefront = t2Dplatefront+deltaIDplate;
        FE_Grid_outerPlatefront.ID = [1:length(p2Dplatefront)]' + deltaIDplate;
        FE_Grid_outerPlatefront.X = p2Dplatefront(:,1);
        FE_Grid_outerPlatefront.Y = p2Dplatefront(:,2);
        FE_Grid_outerPlatefront.Z = interp1([coordpointsUP(end,1),sparFrontPosPLATE],[leading_edge_Z, spar_Front_Z], FE_Grid_outerPlatefront.X);%ones(length(p2Dplate),1)*coordpointsUP(1,3);
        tolPlate = 0.00001;
        tipplateDoubleA = [];
        
        for j = 1:length(FE_Grid_outerPlatefront.X)
            tipplateDouble = ~isempty(FE_Grid.ID(abs(FE_Grid.X - FE_Grid_outerPlatefront.X(j))+abs(FE_Grid.Y - FE_Grid_outerPlatefront.Y(j))+abs(FE_Grid.Z - FE_Grid_outerPlatefront.Z(j)) < tolPlate));
            if tipplateDouble
                t2Dplatefront(t2Dplatefront==FE_Grid_outerPlatefront.ID(j)) = FE_Grid.ID(abs(FE_Grid.X - FE_Grid_outerPlatefront.X(j))+abs(FE_Grid.Y - FE_Grid_outerPlatefront.Y(j))+abs(FE_Grid.Z - FE_Grid_outerPlatefront.Z(j))<tolPlate);
                tipplateDoubleA = [tipplateDoubleA;j];
            end
        end
        FE_Grid_outerPlatefront.ID(tipplateDoubleA) = [];
        FE_Grid_outerPlatefront.X(tipplateDoubleA) = [];
        FE_Grid_outerPlatefront.Y(tipplateDoubleA) = [];
        FE_Grid_outerPlatefront.Z(tipplateDoubleA) = [];
        
        FE_Grid.ID = [FE_Grid.ID; FE_Grid_outerPlatefront.ID];
        FE_Grid.X = [FE_Grid.X; FE_Grid_outerPlatefront.X];
        FE_Grid.Y = [FE_Grid.Y; FE_Grid_outerPlatefront.Y];
        FE_Grid.Z = [FE_Grid.Z; FE_Grid_outerPlatefront.Z];
        
        ID_frontoutEL = [1:length(t2Dplatefront)]' + deltaIDplate;
        ctria_out = FE_OUT_ctria3([], ID_frontoutEL, rigidRibPID, t2Dplatefront(:,1), t2Dplatefront(:,2), t2Dplatefront(:,3));
        meshOut.ctria = [meshOut.ctria; ctria_out];
    end
    
    %%%
    if addMidPlate
        deltaIDplate = deltaIDplate+4430*i;
        t2DplateMID = t2DplateMID+deltaIDplate;
        FE_Grid_outerPlatefront.ID = [1:length(p2DplateMID)]' + deltaIDplate;
        FE_Grid_outerPlatefront.X = p2DplateMID(:,1);
        FE_Grid_outerPlatefront.Y = p2DplateMID(:,2);
        FE_Grid_outerPlatefront.Z = ones(length(p2DplateMID),1)*coordpointsUP(1,3);
        tolPlate = 0.00001;
        tipplateDoubleA = [];
        
        for j = 1:length(FE_Grid_outerPlatefront.X)
            tipplateDouble = ~isempty(FE_Grid.ID(abs(FE_Grid.X - FE_Grid_outerPlatefront.X(j))+abs(FE_Grid.Y - FE_Grid_outerPlatefront.Y(j))+abs(FE_Grid.Z - FE_Grid_outerPlatefront.Z(j)) < tolPlate));
            if tipplateDouble
                t2DplateMID(t2DplateMID==FE_Grid_outerPlatefront.ID(j)) = FE_Grid.ID(abs(FE_Grid.X - FE_Grid_outerPlatefront.X(j))+abs(FE_Grid.Y - FE_Grid_outerPlatefront.Y(j))+abs(FE_Grid.Z - FE_Grid_outerPlatefront.Z(j))<tolPlate);
                tipplateDoubleA = [tipplateDoubleA;j];
            end
        end
        FE_Grid_outerPlatefront.ID(tipplateDoubleA) = [];
        FE_Grid_outerPlatefront.X(tipplateDoubleA) = [];
        FE_Grid_outerPlatefront.Y(tipplateDoubleA) = [];
        FE_Grid_outerPlatefront.Z(tipplateDoubleA) = [];
        
        FE_Grid.ID = [FE_Grid.ID; FE_Grid_outerPlatefront.ID];
        FE_Grid.X = [FE_Grid.X; FE_Grid_outerPlatefront.X];
        FE_Grid.Y = [FE_Grid.Y; FE_Grid_outerPlatefront.Y];
        FE_Grid.Z = [FE_Grid.Z; FE_Grid_outerPlatefront.Z];

        ID_frontoutEL = [1:size(t2DplateMID,1)]' + deltaIDplate;
        ctria_out = FE_OUT_ctria3([], ID_frontoutEL, rigidRibPID, t2DplateMID(:,1), t2DplateMID(:,2), t2DplateMID(:,3));
        meshOut.ctria = [meshOut.ctria; ctria_out];
    end
    
end


%% MIRROR WING
dID = 1000000;
FE_GridM.ID = FE_Grid.ID + dID;
FE_GridM.X = FE_Grid.X;
FE_GridM.Y = FE_Grid.Y;
deltaZ_FE_GridZ = min(FE_Grid.Z);
FE_Grid.Z = FE_Grid.Z - deltaZ_FE_GridZ;
FE_GridM.Z = -FE_Grid.Z;

FE_GridM_double = find(FE_GridM.Z == 0);
FE_GridM_double_ID = FE_GridM.ID(FE_GridM_double);
FE_GridM.ID(FE_GridM_double) = [];
FE_GridM.X(FE_GridM_double) = [];
FE_GridM.Y(FE_GridM_double) = [];
FE_GridM.Z(FE_GridM_double) = [];

FE_Grid.ID = [FE_Grid.ID; FE_GridM.ID];
FE_Grid.X = [FE_Grid.X; FE_GridM.X];
FE_Grid.Y = [FE_Grid.Y; FE_GridM.Y];
FE_Grid.Z = [FE_Grid.Z; FE_GridM.Z];

FE_SkinNodeListUM = FE_SkinNodeListU + dID;
FE_SkinNodeListDM = FE_SkinNodeListD + dID;

for i = 1:length(FE_GridM_double_ID)
    FE_SkinNodeListUM(find(FE_SkinNodeListUM == FE_GridM_double_ID(i))) = FE_SkinNodeListUM(find(FE_SkinNodeListUM == FE_GridM_double_ID(i))) - dID;
    FE_SkinNodeListDM(find(FE_SkinNodeListDM == FE_GridM_double_ID(i))) = FE_SkinNodeListDM(find(FE_SkinNodeListDM == FE_GridM_double_ID(i))) - dID;
end

FE_SkinNodeListU = [FE_SkinNodeListU; FE_SkinNodeListUM];
FE_SkinNodeListD = [FE_SkinNodeListD; FE_SkinNodeListDM];
wingModelStructure.SkinNodeListU = unique(FE_SkinNodeListU);
wingModelStructure.SkinNodeListD = unique(FE_SkinNodeListD);

SPC_node_listM = SPC_node_list + dID;
SPC_node_list = [SPC_node_list, SPC_node_listM];
SPC_node_list = unique(SPC_node_list);

FE_Grid_rev = sparse([FE_Grid.ID], ones(length(FE_Grid.ID),1), (1:length(FE_Grid.ID)), max(FE_Grid.ID), 1, length(FE_Grid.ID));


%% Actuation positions
for iRib = numRibs-2:-2:numRibs-2*wingDesign.nRibsC
    wingModelStructure.actuation.actuationSparAP(iRib,:) = [rib(iRib).actuation.location_spar_Pos(1), rib(iRib).actuation.location_spar_Pos(2), rib(iRib).actuation.location_spar_Pos(3)];
    wingModelStructure.actuation.actuationsdvAP(iRib,:) = [rib(iRib).actuation.location_sdv_Pos(1), rib(iRib).actuation.location_sdv_Pos(2), rib(iRib).actuation.location_sdv_Pos(3)];
    wingModelStructure.actuation.dist_0(iRib) = rib(iRib).actuation.rodlength;
end
if wingDesign.nRibsC > 0
    wingModelStructure.actuation.actuationSparAP(:,3) = wingModelStructure.actuation.actuationSparAP(:,3) + extraRibOffset;
    wingModelStructure.actuation.actuationsdvAP(:,3) = wingModelStructure.actuation.actuationsdvAP(:,3) + extraRibOffset;
end
       

%% Add skin properties ID

SkinProperties = [];
SkinProperties = [SkinProperties;[PID_FS, MID_FS, wingDesign.flexSkinT]];
SkinProperties = [SkinProperties;[rigidRibPID, MID_Skin, wingDesign.rigidRibT]];
SkinProperties = [SkinProperties;[actuationID, MID_Skin, wingDesign.actuationT]];
SkinProperties = [SkinProperties;[PID_Skin, MID_Skin, wingDesign.t_skin]];
SkinProperties = [SkinProperties;[PID_Spar, MID_Spar, wingDesign.t_spar]];


%% run YetAnotherFECode to calculate the stiffness and mass matrix

wingProperties.FE_Grid = [FE_Grid.ID,FE_Grid.X,FE_Grid.Y,FE_Grid.Z];
wingProperties.SkinProperties = SkinProperties;
wingProperties.quadElements = meshOut.cquad;
wingProperties.triaElements = meshOut.ctria;
wingProperties.SPC_nodes = SPC_node_list;
wingProperties.isPIDout = [PID_IS, MID_IS, wingDesign.voronoiT];
wingProperties.sparPIDout = [PID_Spar, MID_Spar, wingDesign.t_spar];
wingProperties.rigidRibPIDout = [rigidRibPID, MID_Skin, wingDesign.rigidRibT];
wingProperties.actuationPIDout = [actuationID, MID_Skin, wingDesign.actuationT];
wingProperties.skinPIDout = [PID_Skin, MID_Skin, wingDesign.t_skin];
wingProperties.fsPIDout = [PID_FS, MID_FS, wingDesign.flexSkinT];

wingProperties.PID.MID_Skin = MID_Skin;
wingProperties.PID.MID_IS = MID_IS;
wingProperties.PID.MID_FS = MID_FS;

[M,K,wP] = runYetAnotherFEcode(wingProperties, wingDesign);

wingModelStructure.wP = wP; % rigid wing properties mass, inertia, cog


%% Setting constrains
nDOF = size(K,1);
SPC_node_list = unique(SPC_node_list);
GID_ASET = [1:(nDOF/6)]';
GID_ASET(FE_Grid_rev(SPC_node_list)) = [];
DOF_ASET = repmat(0:5, size(GID_ASET,1), 1) + repmat((GID_ASET-1)*6+1, 1, 6);

DOF_RSET = full(FE_Grid_rev(SPC_node_list));
DOF_RSET = repmat(0:5, size(DOF_RSET,1), 1) + repmat((DOF_RSET-1)*6+1, 1, 6);
DOF_RSET = DOF_RSET(:); % SPC_mat_ID contains the DOFs in the R-SET

K_ASET = K;
K_ASET(:, DOF_RSET) = [];
K_ASET(DOF_RSET, :) = [];

M_ASET = M;
M_ASET(:, DOF_RSET) = [];
M_ASET(DOF_RSET, :) = [];


% precalculate transpose and (:)
DOF_ASETt = DOF_ASET';
DOF_ASET = DOF_ASETt(:);
wingModelStructure.DOF_ASET                  = DOF_ASET;    
wingModelStructure.DOF_RSET                  = DOF_RSET;               
wingModelStructure.FE_Grid       = FE_Grid;  

wingModelStructure.M_ASET = M_ASET; % Full mass matrix
wingModelStructure.K_ASET = K_ASET; % Full stuffness matrix


%% Get actuation parameters initialization and actuation modes
forceACTinp.L = -simParam.actForce ;
forceACTinp.R = -simParam.actForce ;
[L_ACT, R_ACT, ACTparam] = getActuationModes(wingModelStructure, wingDesign, numRibs, FE_Grid, K_ASET, DOF_RSET, forceACTinp);
wingModelStructure.ACTparam = ACTparam;
if simParam.addMM
    ACT = [L_ACT, R_ACT];
else
    ACT = [];
    simParam.nmodes = simParam.nmodes + 2;
end


%% Get vibration modes
[phi_ASET, omegaSquared] = eigs(K_ASET, -M_ASET, simParam.nmodes, 'sm');

% natural frequencies  
wingModelStructure.omega = imag(diag(omegaSquared).^.5./(2*pi));

% modal stiffness matrix
phi_ASET_T = [phi_ASET, ACT];         %Add actuation mode to deformation mode
K_MODES = phi_ASET_T'*K_ASET*phi_ASET_T;
wingModelStructure.K_MODES = K_MODES;
wingModelStructure.phi_ASET = phi_ASET_T;

% modal mass matrix
M_MODES = phi_ASET_T'*M_ASET*phi_ASET_T;
wingModelStructure.M_MODES = M_MODES;

% modal damping matrix
epsilonD2 = 0.005; % damping ratio https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/usb/default.htm?startat=pt05ch20s01abm43.html
alphaDomega12 = 2*epsilonD2*wingModelStructure.omega(1)*wingModelStructure.omega(2)*2*pi/(wingModelStructure.omega(1)+wingModelStructure.omega(2));
betaDomega12 = 2*epsilonD2/(wingModelStructure.omega(1)*2*pi+wingModelStructure.omega(2)*2*pi);   % Alpha and beta parameters for dumping, see UF aeroactuationelastic ptimization
wingModelStructure.C_ASET = alphaDomega12*M_ASET + betaDomega12*K_ASET;
wingModelStructure.C_MODES = alphaDomega12*M_MODES + betaDomega12*K_MODES;



