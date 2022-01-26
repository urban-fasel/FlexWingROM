function [FE_nodes_3D, FE_iNodes_up_CCW, FE_iNodes_dn_CCW, FE_nodes_spar, FE_delta_iNodes_EdgeToEdge, corner_nodes, FE_Grid, FE_nodes_spar2, meshOut] = createRibStructure(meshOut, wingDesign, targetShape, voronoi_settings, domain_settings, rib_width, nodeOffset, elemOffset, Z_coord0, numElemThickness, domain_segments_PID, createOnlyOuterNodes, PID_FS, PID_IS)

%%
%
%  By Giulio Molinari
%
%  Modified by Urban Fasel
%
%

%% Constants
interp1method = 'spline'; % interpolation method for fitting all points inside the profile

%% Init
k_nodespacing = 0.75; % 0: linear, 1: 1-cosine
voronoi_ext_merge_tol = 0.0025; % merging distance, m

sites_X = voronoi_settings.sites_X;
sites_Y = voronoi_settings.sites_Y;
sites_thickness = voronoi_settings.sites_thickness;

domain_bounds = domain_settings.bounds; % x1, y1; x2, y2; x3, y3
elem_max_size = domain_settings.elem_max_size;
elem_min_number = domain_settings.elem_min_number;
elem_max_size_ext = domain_settings.elem_max_size_ext;
elem_min_number_ext = domain_settings.elem_min_number_ext;
domain_segments_thickness = domain_settings.domain_segments_thickness;
elem_max_size_LE = domain_settings.elem_max_size_LE;

domain_segments_x = [domain_bounds(1,1), domain_bounds(2,1), domain_bounds(3,1), domain_bounds(2,1), domain_bounds(1,1);
                     domain_bounds(2,1), domain_bounds(3,1), domain_bounds(2,1), domain_bounds(1,1), domain_bounds(1,1)];
domain_segments_y = [domain_bounds(1,2), domain_bounds(1,2), domain_bounds(2,2), domain_bounds(3,2), domain_bounds(3,2);
                     domain_bounds(1,2), domain_bounds(2,2), domain_bounds(3,2), domain_bounds(3,2), domain_bounds(1,2)];

%targetShape_nodes = targetShape.nodes;
% note: first point is repeated in upper and lower section
targetShape_nodes = [targetShape.nodes(:, 1:2); targetShape.nodes(2:end, 3:4)];
[xLE, iLE] = min(targetShape_nodes(:, 1));
targetShape_up = targetShape_nodes(1:iLE, :);
targetShape_dn = targetShape_nodes(iLE:end, :);

targetShape_X_start = targetShape.X_start;
targetShape_X_stop = targetShape.X_stop;
targetShape_dX = targetShape_X_stop - targetShape_X_start;



yDomainBottom =@ (x)	(domain_bounds(1,2).*(x>=domain_bounds(1,1)).*(x<=domain_bounds(2,1))) + ...
						((domain_bounds(1,2) + (domain_bounds(2,2)-domain_bounds(1,2)).*(x - domain_bounds(2,1))./(domain_bounds(3,1) - domain_bounds(2,1))).*(x>domain_bounds(2,1)).*(x<=domain_bounds(3,1)));
yDomainTop =@ (x)		(domain_bounds(3,2).*(x>=domain_bounds(1,1)).*(x<=domain_bounds(2,1))) + ...
						((domain_bounds(3,2) + (domain_bounds(2,2)-domain_bounds(3,2)).*(x - domain_bounds(2,1))./(domain_bounds(3,1) - domain_bounds(2,1))).*(x>domain_bounds(2,1)).*(x<=domain_bounds(3,1)));

dt = delaunayTriangulation([sites_X, sites_Y]);
[vx, vy, c_combined, cells] = voronoi_triangulation(dt, '', 100);


%% Limit to boundaries
% create segments array
segments = [];
segments_cells = [];
for i = 1:length(cells)
 	for j = 1:size(cells{i}, 2)
		if numel(segments) > 0
			iSegm = -1;
			for k = 1:size(segments,2)
				if sum(segments(:, k) == cells{i}(:, j))==2 || sum(segments([2,1], k) == cells{i}(:, j))==2
					iSegm = k;
					break
				end
			end
			if iSegm ~= -1
				% this segment is already in the array
			else
				segments = [segments, cells{i}(:, j)];
				iSegm = size(segments, 2);
			end
		else
			segments = [segments, cells{i}(:, j)];
			iSegm = size(segments, 2);
		end
		try
			if segments_cells(iSegm, 1) == 0
				segments_cells(iSegm, 1) = i;
			else
				segments_cells(iSegm, 2) = i;
			end
		catch
			segments_cells(iSegm, 1) = i;
		end
	end
end


% check for intersections
% for all segments
intersections = {};
intersections{size(domain_segments_x,2)} = []; % was: intersections{size(segments,2)} = [];
for i = 1:size(segments,2)
	node1 = segments(1, i);
	node2 = segments(2, i);
	% get the intersection (if exists) with the boundaries of the domain
	p = [c_combined(node1, 1); c_combined(node1, 2)];
	dp = [c_combined(node2, 1); c_combined(node2, 2)] - p;
	for j = 1:size(domain_segments_x,2)
		q = [domain_segments_x(1, j); domain_segments_y(1, j)];
		dq = [domain_segments_x(2, j); domain_segments_y(2, j)] - q;
		[t, u] = segment_intersect_2d(p, dp, q, dq);
		if ((t >= 0 && t <= 1) && (u >= 0 && u <= 1))
			% there's an intersection!
			% identify which of the two points (p or p+dp) is outside the domain
			if twodcross([dq(1), dq(2)], [p(1)-q(1), p(2)-q(2)]) < 0
				% first point is outside the domain
				nodeOut = node1;
			else
				% second point is outside the domain
				nodeOut = node2;
			end
			
			% check all other segments having nodeOut as extreme point and create a copy
			needToCopy = false;
			for k = 1:size(segments,2)
				if k ~= i
					if segments(1,k) == nodeOut
						segments(1,k) = size(c_combined, 1) + 1;
						needToCopy = true;
					end
					if segments(2,k) == nodeOut
						segments(2,k) = size(c_combined, 1) + 1;
						needToCopy = true;
					end
				end
			end
			
			if needToCopy
				c_combined = [c_combined; c_combined(nodeOut,:)];
			end
			
			c_combined(nodeOut, 1) = q(1) + u*dq(1);
			c_combined(nodeOut, 2) = q(2) + u*dq(2);

			try
				intersections{j} = [intersections{j} [nodeOut; u]];
			catch
				intersections{j} = [nodeOut; u];
			end
		end
    end
end


% check for segments completely outside the domain (and remove them)
% for all segments
for i = 1:size(segments,2)
	node1 = segments(1, i);
	node2 = segments(2, i);
	if ((~inpolygon(c_combined(node1,1), c_combined(node1,2), domain_segments_x(1,:), domain_segments_y(1,:))) && (~inpolygon(c_combined(node2,1), c_combined(node2,2), domain_segments_x(1,:), domain_segments_y(1,:))))
		c_combined(node1,1) = NaN;
		c_combined(node1,2) = NaN;
		c_combined(node2,1) = NaN;
		c_combined(node2,2) = NaN;
	end
end


% Add external domain nodes to the list of nodes
c_combined_start_extnodes = size(c_combined, 1);
for i = 1:size(domain_segments_x, 2)
	c_combined(c_combined_start_extnodes+i, :) = [domain_segments_x(1, i), domain_segments_y(1, i)];
end

segments_ext = [];
segments_ext_parent = [];
% Add segments of the external domain
for i = 1:size(domain_segments_x, 2)
	node1 = c_combined_start_extnodes + i;
	intersections_segment = intersections{i};
	if size(intersections_segment,2) > 0
		[~, iSort] = sort(intersections_segment(2,:));
		intersections_segment_sorted = intersections_segment(1,iSort);
		
		% if two intersecting nodes are VERY close, merge them together and keep the last one
		intersections_pos_sorted = intersections_segment(2,iSort);
		for j = 1:size(intersections_pos_sorted,2)-1
			if intersections_pos_sorted(j+1)-intersections_pos_sorted(j) < voronoi_ext_merge_tol/targetShape_dX
				old_nodeID = intersections_segment_sorted(j);
				new_nodeID = intersections_segment_sorted(j+1);
				segments(segments == old_nodeID) = new_nodeID;
				
				c_combined(new_nodeID, :) = (c_combined(old_nodeID, :) + c_combined(new_nodeID, :))/2;
				c_combined(old_nodeID, :) = NaN;
				
				intersections_segment_sorted(j) = NaN;
			end
		end
		intersections_segment_sorted = intersections_segment_sorted(~isnan(intersections_segment_sorted));
		
		for j = 1:size(intersections_segment_sorted,2)
			node2 = intersections_segment_sorted(j);
			segments_ext = [segments_ext [node1; node2]];
			segments_ext_parent = [segments_ext_parent i];
			node1 = node2;
		end
	end
	node2 = c_combined_start_extnodes + i + 1;
	segments_ext = [segments_ext [node1; node2]];
	segments_ext_parent = [segments_ext_parent i];
end
segments_ext(2,end) = c_combined_start_extnodes + 1;


%% Squeeze internal nodes to fit in the profile
if createOnlyOuterNodes
	iOuterNodes = unique(segments_ext);
	iSelectOuterNodes = zeros(size(c_combined,1),1);
	iSelectOuterNodes(iOuterNodes) = 1;
	x_nodes = c_combined(:, 1);
	y_nodes = c_combined(:, 2);
	x_nodes(iSelectOuterNodes == 0) = NaN;
	y_nodes(iSelectOuterNodes == 0) = NaN;
else
	x_nodes = c_combined(:, 1);
	y_nodes = c_combined(:, 2);
end
x_rel_nodes = (x_nodes - domain_segments_x(1,1))/(domain_segments_x(1,3) - domain_segments_x(1,1) + eps); % eps added to the denominator to fix 0/0. Now it's 0/eps and Trailing Edge point works...

y_rel_nodes = (y_nodes - yDomainBottom(x_rel_nodes))./(yDomainTop(x_rel_nodes) - yDomainBottom(x_rel_nodes));

x_nodes_new = targetShape_X_start + x_rel_nodes*targetShape_dX;
y_nodes_new = interp1(targetShape_dn(:,1), targetShape_dn(:,2), x_nodes_new, interp1method) + ...
			  y_rel_nodes.*(interp1(targetShape_up(:,1), targetShape_up(:,2), x_nodes_new, interp1method) - interp1(targetShape_dn(:,1), targetShape_dn(:,2), x_nodes_new, interp1method));


%% Generate a structural mesh
FE_nodes = [x_nodes_new, y_nodes_new]; % X1, X2
FE_beams = []; % GA, GB
FE_beams_thickness_values = [];
FE_beams_iThickness = [];
FE_nodes_ext_up = [];
FE_nodes_ext_dn = [];
FE_nodes_spar = [];
FE_nodes_parent = {};

if ~createOnlyOuterNodes
	% first on the internal segments
	i_notNaN = 0;
	for i = 1:size(segments,2)
		node1 = segments(1, i);
		node2 = segments(2, i);
		p = [x_nodes_new(node1); y_nodes_new(node1)];
		dp = [x_nodes_new(node2); y_nodes_new(node2)] - p;
		if sum(isnan(p))==0 && sum(isnan(dp))==0
			i_notNaN = i_notNaN + 1;
			num_subdivisions = round(max([elem_min_number, sqrt(dp(1)^2+dp(2)^2)/elem_max_size]));
			dp_subd = dp*linspace(0, 1, num_subdivisions + 1);
			iNodeStartBeam = size(FE_nodes, 1) + 1;
			FE_nodes = [FE_nodes; repmat(p', num_subdivisions-1, 1) + dp_subd(:,2:end-1)'];
			FE_beams = [FE_beams; [node1 iNodeStartBeam + (0:num_subdivisions-2); [iNodeStartBeam + (0:num_subdivisions-2) node2]]'];

			cur_thickness = sites_thickness(segments_cells(i, 1)) + sites_thickness(segments_cells(i, 2));
			FE_beams_thickness_values = [FE_beams_thickness_values cur_thickness];
			FE_beams_iThickness = [FE_beams_iThickness repmat(i_notNaN, 1, num_subdivisions)];
		end
	end
	segments_notNaN = i_notNaN;
else
	segments_notNaN = 0;
end


% then on the external segments, and move their nodes so that they match the aerodyn. profile
i_notNaN = 0;
FE_beams_thickness_values = [FE_beams_thickness_values(:); domain_segments_thickness(:)];
for i = 1:size(segments_ext,2)
	node1 = segments_ext(1, i);
	node2 = segments_ext(2, i);
	p = [x_nodes_new(node1); y_nodes_new(node1)];
	dp = [x_nodes_new(node2); y_nodes_new(node2)] - p;
	if sum(isnan(p))==0 && sum(isnan(dp))==0
		i_notNaN = i_notNaN + 1;
		num_subdivisions = round(max([elem_min_number, sqrt(dp(1)^2+dp(2)^2)/elem_max_size]));
		dp_subd = dp*linspace(0, 1, num_subdivisions + 1);
		iNodeStartBeam = size(FE_nodes, 1) + 1;
		switch segments_ext_parent(i) % 1: 
			%%%% TO BE CHANGED WHEN THE EXTERNAL SEGMENTS ARE MODIFIED ####
			case {1, 2}
				% lower
				num_subdivisions = round(max([elem_min_number_ext, sqrt(dp(1)^2+dp(2)^2)/elem_max_size_ext]));
				dp_subd = dp*linspace(0, 1, num_subdivisions + 1);
				FE_seg_nodes_x = p(1) + dp_subd(1, 2:end-1);
				FE_seg_nodes_y = interp1(targetShape_dn(:,1), targetShape_dn(:,2), FE_seg_nodes_x, interp1method);
				FE_nodes = [FE_nodes; [FE_seg_nodes_x' FE_seg_nodes_y']];
				FE_nodes_ext_dn = [FE_nodes_ext_dn [node1 iNodeStartBeam + (0:num_subdivisions-2) node2]];
			case {3, 4}
				% upper
				num_subdivisions = round(max([elem_min_number_ext, sqrt(dp(1)^2+dp(2)^2)/elem_max_size_ext]));
				dp_subd = dp*linspace(0, 1, num_subdivisions + 1);
				FE_seg_nodes_x = p(1) + dp_subd(1, 2:end-1);
				FE_seg_nodes_y = interp1(targetShape_up(:,1), targetShape_up(:,2), FE_seg_nodes_x, interp1method);
				FE_nodes = [FE_nodes; [FE_seg_nodes_x' FE_seg_nodes_y']];
				FE_nodes_ext_up = [FE_nodes_ext_up [node1 iNodeStartBeam + (0:num_subdivisions-2) node2]];
			otherwise
				FE_nodes = [FE_nodes; repmat(p', num_subdivisions-1, 1) + dp_subd(:,2:end-1)'];
				if segments_ext_parent(i) == 5; FE_nodes_spar = [FE_nodes_spar [node1 iNodeStartBeam + (0:num_subdivisions-2) node2]]; end %%% SPAR is segment 5
		end
		if length(FE_nodes_parent) < node1
			FE_nodes_parent{node1} = segments_ext_parent(i);
		else
 			if sum([FE_nodes_parent{node1}] == segments_ext_parent(i)) == 0
				FE_nodes_parent{node1} = [FE_nodes_parent{node1} segments_ext_parent(i)];
			end
		end
		if length(FE_nodes_parent) < node2
			FE_nodes_parent{node2} = segments_ext_parent(i);
		else
			if sum([FE_nodes_parent{node2}] == segments_ext_parent(i)) == 0
				FE_nodes_parent{node2} = [FE_nodes_parent{node2} segments_ext_parent(i)];
			end
		end
		for j = 1:num_subdivisions-1
			if length(FE_nodes_parent) < iNodeStartBeam+(j-1)
				FE_nodes_parent{iNodeStartBeam+(j-1)} = segments_ext_parent(i);
			else
				if sum([FE_nodes_parent{iNodeStartBeam+(j-1)}] == segments_ext_parent(i)) == 0
					FE_nodes_parent{iNodeStartBeam+(j-1)} = [FE_nodes_parent{iNodeStartBeam+(j-1)} segments_ext_parent(i)];
				end
			end
		end
		
		FE_beams = [FE_beams; [node1 iNodeStartBeam + (0:num_subdivisions-2); [iNodeStartBeam + (0:num_subdivisions-2) node2]]'];
	else
		warning('Something went wrong: nodes on the boundaries have NaN values!');
	end
	FE_beams_iThickness = [FE_beams_iThickness repmat(segments_notNaN + segments_ext_parent(i), 1, num_subdivisions)];
end


% Identify corner nodes
%%% Corner nodes will be: CORRUGATION START BOTTOM, CORRUGATION END (= Trailing Edge), CORRUGATION START TOP, SPAR TOP, SPAR BOTTOM, LEADING EDGE (added later)
corner_nodes = zeros(6,1);
for i = 1:length(FE_nodes_parent)
	if length(FE_nodes_parent{i}) > 1
		for j = 1:length(FE_nodes_parent{i})
			if FE_nodes_parent{i}(j) == 5
				k = find(FE_nodes_parent{i} == 1, 1);
			else
				k = find(FE_nodes_parent{i} == FE_nodes_parent{i}(j) + 1, 1);
			end
			if ~isempty(k)
				if k ~= 0 && corner_nodes(FE_nodes_parent{i}(j)) == 0
					corner_nodes(FE_nodes_parent{i}(j)) = i;
				end
			end
		end
	end
end

% then on the leading edge of the profile
% FROM 0 TO targetShape_X_start, STEP elem_max_size_ext, ON targetShape_up AND targetShape_dn
num_FE_nodes = size(FE_nodes, 1);
numElem_LE = ceil(targetShape_X_start/elem_max_size_LE);
numPropertyLE_UP = max(FE_beams_iThickness) + 1;
numPropertyLE_DN = max(FE_beams_iThickness) + 2;

x_subd_LE = xLE + ((1-k_nodespacing)*linspace(0, 1, numElem_LE) + k_nodespacing*(1-cos(linspace(0, 1, numElem_LE)*pi/2)))*(targetShape_X_start - xLE); % MIXED NODE SPACING

y_subd_LE_up = interp1(targetShape_up(:,1), targetShape_up(:,2), x_subd_LE, 'linear', 'extrap');
y_subd_LE_dn = interp1(targetShape_dn(:,1), targetShape_dn(:,2), x_subd_LE, 'linear', 'extrap');
FE_nodes_LE_up = [fliplr(x_subd_LE(1:end-1))', fliplr(y_subd_LE_up(1:end-1))'];
FE_nodes_LE_dn = [(x_subd_LE(2:end-1))', (y_subd_LE_dn(2:end-1))'];
FE_nodes_LE = [FE_nodes_LE_up; FE_nodes_LE_dn]; % arranged counter-clockwise
FE_nodes = [FE_nodes; FE_nodes_LE];

spar2_pos = wingDesign.frontSparX*targetShape_X_start;
[~,iNode] = min(abs(FE_nodes_LE_up(:,1)-spar2_pos));
[~,iNode2] = min(abs(FE_nodes_LE_dn(:,1)-spar2_pos));

for i = 0:length(FE_nodes_LE)
    if i == iNode 
        Node_top = num_FE_nodes+i;
    end
    if i == iNode2+length(FE_nodes_LE_up)
        Node_bot = num_FE_nodes+i;
    end
	if i >= 1
		node1 = num_FE_nodes + i;
	else
		% the first element starts from the top of the main spar, CORNER_NODES(4)
		node1 = corner_nodes(4);
		if length(FE_nodes_parent) < node1
			FE_nodes_parent{node1} = 6; % Parent: LEADING EDGE PORTION = 6
		else
			FE_nodes_parent{node1} = [FE_nodes_parent{node1} 6]; % Parent: also LEADING EDGE PORTION = 6
		end
	end
	if i < length(FE_nodes_LE)
		node2 = num_FE_nodes + i + 1;
	else
		% the last element ends at the bottom of the main spar, CORNER_NODES(3)
		node2 = corner_nodes(5);
	end
	FE_beams = [FE_beams; node1, node2];

	if length(FE_nodes_parent) < node2
		FE_nodes_parent{node2} = 6; % Parent: LEADING EDGE PORTION = 6
	else
		FE_nodes_parent{node2} = [FE_nodes_parent{node2} 6]; % Parent: also LEADING EDGE PORTION = 6
	end
end
FE_beams_iThickness = [FE_beams_iThickness repmat(numPropertyLE_UP, 1, numElem_LE-1) repmat(numPropertyLE_DN, 1, numElem_LE-1)];
FE_nodes_LE = num_FE_nodes + numElem_LE - 1;

corner_nodes(6) = FE_nodes_LE;

boundary_nodes = {};
for i = 1:length(FE_nodes_parent)
	for j = 1:length(FE_nodes_parent{i})
		if length(boundary_nodes) < FE_nodes_parent{i}(j)
			boundary_nodes{FE_nodes_parent{i}(j)} = i;
		else
			boundary_nodes{FE_nodes_parent{i}(j)} = [boundary_nodes{FE_nodes_parent{i}(j)} i];
		end
	end
end

% sort nodes on the bottom, unique and in ascending order
botNodeIDs = unique([boundary_nodes{1} boundary_nodes{2}]);
[~, iBN1s] = sort(FE_nodes(botNodeIDs, 1), 1, 'ascend');
BN1s = botNodeIDs(iBN1s);

% sort nodes on the top, unique and in descending order
topNodeIDs = unique([boundary_nodes{3} boundary_nodes{4}]);
[~, iBN3s] = sort(FE_nodes(topNodeIDs, 1), 1, 'descend');
BN3s = topNodeIDs(iBN3s);

FE_iNodes_up_CCW = [BN3s, ((num_FE_nodes+1):FE_nodes_LE)];
FE_iNodes_dn_CCW = [(FE_nodes_LE):size(FE_nodes, 1), BN1s]; % the LE node is repeated

% Remove duplicates in the spar
FE_nodes_spar = unique(FE_nodes_spar);
% Sort spar nodes from bottom to top
[~, iSortSpar] = sort(FE_nodes(FE_nodes_spar,2), 'ascend');
FE_nodes_spar = FE_nodes_spar(iSortSpar);


%% Add second spar
[~,iNode] = min(abs(FE_nodes_LE_up(:,1)-spar2_pos));
[~,iNode2] = min(abs(FE_nodes_LE_dn(:,1)-spar2_pos));
p = FE_nodes_LE_up(iNode,:)';
dp = FE_nodes_LE_dn(iNode2,:)'-p;
num_subdivisions2 = round(max([elem_min_number, sqrt(dp(1)^2+dp(2)^2)/elem_max_size]));
dp_subd = dp*linspace(0, 1, num_subdivisions2 + 1);
Nodes_Spar2 = [FE_nodes_LE_up(iNode,:); repmat(p', num_subdivisions2-1, 1) + dp_subd(:,2:end-1)'; FE_nodes_LE_dn(iNode2,:)];
ID_nodes_spar2 = [Node_top length(FE_nodes)+(1:length(Nodes_Spar2)-2) Node_bot];
FE_nodes = [FE_nodes; Nodes_Spar2]; 
FE_nodes = unique(FE_nodes,'stable','rows');
for i = 1:length(ID_nodes_spar2)-1
    FE_beams = [FE_beams; ID_nodes_spar2(i), ID_nodes_spar2(i+1)];
end
FE_nodes_spar2.ID = ID_nodes_spar2;
FE_nodes_spar2.pos = Nodes_Spar2;

numProperty_Spar2 = segments_notNaN+5;
FE_beams_iThickness = [FE_beams_iThickness repmat(numProperty_Spar2, 1, num_subdivisions2)];

%% Extrusion of the rib along the third dimension (generation of new points and QUAD elements)
FE_delta_iNodes_EdgeToEdge = size(FE_nodes, 1)*(numElemThickness);

FE_nodes = [FE_nodes(:,1), FE_nodes(:,2), repmat(Z_coord0, size(FE_nodes, 1), 1)];

numNodesPerSlice = length(FE_nodes);
numBeamsPerSlice = length(FE_beams);

ribSlicesRelLocation = linspace(-.5, .5, numElemThickness+1);
FE_nodes_3D = ones(size(FE_nodes,1) * (numElemThickness+1), 3)*NaN;
for i = 1:numElemThickness+1
	FE_nodes_3D((1+(i-1)*numNodesPerSlice):(i*numNodesPerSlice),:) = [...
		FE_nodes(:,1), ...
		FE_nodes(:,2), ...
		zeros(size(FE_nodes, 1), 1) + ribSlicesRelLocation(i)*rib_width];
end

FE_nodes_3D = FE_nodes_3D+repmat([0 0 Z_coord0], size(FE_nodes_3D, 1), 1);

FE_panels = ones(size(FE_beams,1)*numElemThickness, 4) * NaN;
FE_panels_PID = ones(size(FE_beams,1)*numElemThickness, 1) * NaN;

for i = 1:numElemThickness
	for j = 1:length(FE_beams)
		FE_panels(numBeamsPerSlice*(i-1)+j, 4) = FE_beams(j, 1) + numNodesPerSlice*(i-1);
		FE_panels(numBeamsPerSlice*(i-1)+j, 1) = FE_beams(j, 2) + numNodesPerSlice*(i-1);
		FE_panels(numBeamsPerSlice*(i-1)+j, 2) = FE_beams(j, 2) + numNodesPerSlice*(i);
		FE_panels(numBeamsPerSlice*(i-1)+j, 3) = FE_beams(j, 1) + numNodesPerSlice*(i);
		FE_panels_PID(numBeamsPerSlice*(i-1)+j) = PID_IS; 
		if FE_beams_iThickness(j) > length(FE_beams_thickness_values) - length(domain_segments_thickness)
			FE_panels_PID(numBeamsPerSlice*(i-1)+j) = domain_segments_PID(FE_beams_iThickness(j) - (length(FE_beams_thickness_values) - length(domain_segments_thickness)));
		end
	end
end


%% Renumbering of nodes
FE_nodes_3D_notNaN = ((~isnan(FE_nodes_3D(:,1))) & (~isnan(FE_nodes_3D(:,2))) & (~isnan(FE_nodes_3D(:,3))));
oneToNumNodes = 1:size(FE_nodes_3D_notNaN);
FE_panels_notNaN = FE_nodes_3D_notNaN(FE_panels(:,1)) & FE_nodes_3D_notNaN(FE_panels(:,2)) & FE_nodes_3D_notNaN(FE_panels(:,3)) & FE_nodes_3D_notNaN(FE_panels(:,4));
oneToNumPanels = 1:size(FE_panels_notNaN, 1);


%% Output of panels and nodes

% Output the skin node IDs and coordinates and the node IDs of the upper and lower part of the profile
FE_Grid = [nodeOffset + oneToNumNodes(logical(FE_nodes_3D_notNaN))', FE_nodes_3D(logical(FE_nodes_3D_notNaN),1), FE_nodes_3D(logical(FE_nodes_3D_notNaN),2), FE_nodes_3D(logical(FE_nodes_3D_notNaN),3)];


%% No corrugation at TE 
FE_panels_PIDtest = [];
panelsCoordsX = [FE_nodes_3D(FE_panels(:,1), 1) FE_nodes_3D(FE_panels(:,2), 1) FE_nodes_3D(FE_panels(:,3), 1) FE_nodes_3D(FE_panels(:,4), 1)];
for i=1:length(panelsCoordsX)
    if max(panelsCoordsX(i,:)) > targetShape.X_trans*1.02 && min(panelsCoordsX(i,:)) < targetShape.X_trans2*0.98 && FE_panels_PID(i) == domain_segments_PID(1)
        FE_panels_PID(i) = PID_FS; 
        FE_panels_PIDtest = [FE_panels_PIDtest; PID_FS];
    end
    if FE_panels_PID(i) == 10 % set lower panel PID to PID_Skin if it is not part of the flexible skin
        FE_panels_PID(i) = 1;
    end
end    

%% Create output file
cquad_out = FE_OUT_cquad4(elemOffset + oneToNumPanels(logical(FE_panels_notNaN)), FE_panels_PID(oneToNumPanels(logical(FE_panels_notNaN))), FE_panels(logical(FE_panels_notNaN), 1) + nodeOffset, FE_panels(logical(FE_panels_notNaN), 2) + nodeOffset, FE_panels(logical(FE_panels_notNaN), 3) + nodeOffset, FE_panels(logical(FE_panels_notNaN), 4) + nodeOffset);
meshOut.cquad = [meshOut.cquad; cquad_out];

end
