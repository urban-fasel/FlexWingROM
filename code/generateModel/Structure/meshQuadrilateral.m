function [vertices, origVertices, sectionsNodes] = meshQuadrilateral(xVertices, yVertices, edge1, edge2, edge3, edge4)
%
% y   V1---e4----V4
%     |          |
% ^  e1         e3
% |   |          |
% |  V2---e2----V3
% o-------->   x
%
% edge1...4: [pos_node1; pos_node2...]
%
% EXAMPLE:
% xVertices = [0; 0; 1; 1];
% yVertices = [1; 0; 0; 1];
% 
% edge1 = [0:0.1:1]';
% edge2 = [0:0.1:1]';
% edge3 = [0:0.2:1]';
% edge4 = [0:0.1:1]';
% 
% [vertices, origVertices] = meshQuadrilateral(xVertices, yVertices, edge1, edge2, edge3, edge4);
% 
% vertices = [vertices; origVertices];
% 
% tri = delaunay(vertices(:,1), vertices(:,2));
% triplot(tri, vertices(:,1), vertices(:,2));

vertices = [];
sectionsNodes = nan(max([size(edge1,1) size(edge3,1)]), max([size(edge2,1) size(edge4,1)]));
doPlot = false;

if 4 ~= numel(xVertices) || 4 ~= numel(yVertices)
	error('Must have 4 vertices')
end
if size(edge2, 1) ~= size(edge4, 1)
	error('Edge 2 and edge 4 must have the same number of nodes')
end

num_join_patches = size(edge3, 1) - size(edge1, 1);
%{
if num_join_patches < 0
	xVerticesN = [xVertices(3), xVertices(4), xVertices(1), xVertices(2)];
	yVerticesN = [yVertices(3), yVertices(4), yVertices(1), yVertices(2)];
	% [vertices, origVertices] = meshQuadrilateral(xVerticesN, yVerticesN, edge3, edge4, edge1, edge2);
	[vertices, ~] = meshQuadrilateral(xVerticesN, yVerticesN, edge3, edge4, edge1, edge2);
end
%}
k_addremove = sign(num_join_patches);
num_join_patches = abs(num_join_patches);

pos_join_patches = (1:num_join_patches)/(num_join_patches+1);

edge1V = [(xVertices(1) + (xVertices(2) - xVertices(1)) * edge1(:, 1)), (yVertices(1) + (yVertices(2) - yVertices(1)) * edge1(:, 1))];
edge2V = [(xVertices(2) + (xVertices(3) - xVertices(2)) * edge2(:, 1)), (yVertices(2) + (yVertices(3) - yVertices(2)) * edge2(:, 1))];
edge3V = flipud([(xVertices(3) + (xVertices(4) - xVertices(3)) * edge3(:, 1)), (yVertices(3) + (yVertices(4) - yVertices(3)) * edge3(:, 1))]);
edge4V = flipud([(xVertices(4) + (xVertices(1) - xVertices(4)) * edge4(:, 1)), (yVertices(4) + (yVertices(1) - yVertices(4)) * edge4(:, 1))]);

% origVertices = [edge1V(1:end-1,:); edge2V(1:end-1,:); edge3V(1:end-1,:); edge4V(1:end-1,:)];
origVertices = [edge1V(1:end-1,:); edge2V(1:end-1,:); flipud(edge3V(2:end,:)); flipud(edge4V(2:end,:))];
if num_join_patches < 0
	return;
end

if doPlot
	figure
	hold all
	plot(edge1V(:,1), edge1V(:,2), '.', edge2V(:,1), edge2V(:,2), '.', edge3V(:,1), edge3V(:,2), '.', edge4V(:,1), edge4V(:,2), '.')
end

% create vertical lines
vert_A = zeros(size(edge2, 1), 1);
vert_B = ones(size(edge2, 1), 1);
vert_E = zeros(size(edge2, 1), 1);
for i = 1:size(edge2, 1)
	vert_A(i) = (edge2V(i, 2) - edge4V(i, 2));
	vert_B(i) = (edge4V(i, 1) - edge2V(i, 1));
	vert_E(i) = (edge2V(i, 1)*edge4V(i, 2) - edge4V(i, 1)*edge2V(i, 2));
end

% first "vertical" section = e1
prevVertEdge = edge1V;
iPrevPosX = 0;
numAddedPoints = 0;

for i = 1:num_join_patches+1
	% identify horizontal position of the patch
	if i <= num_join_patches
		iPosX = find(edge2(:,1) - pos_join_patches(i) > 0, 1) - 1;
	else
		iPosX = size(edge2, 1);
	end
	xBot = edge2V(iPosX, 1);
	yBot = edge2V(iPosX, 2);
	xTop = edge4V(iPosX, 1);
	yTop = edge4V(iPosX, 2);

	% get number of elements along y
	num_elem_y = size(edge1, 1) + (i-1)*k_addremove;

	if i <= num_join_patches
		% create horizontal lines from previous section to current section
		% curVertEdge = flipud([(xBot + (xTop - xBot) * linspace(0, 1, num_elem_y)'), (yBot + (yTop - yBot) * linspace(0, 1, num_elem_y)')]); %%%% CONSTANT SPACING
		edge1density = interp1(linspace(0,1,length(edge1)), 1-flipud(edge1), linspace(0, 1, num_elem_y));
		edge3density = interp1(linspace(0,1,length(edge3)), edge3, linspace(0, 1, num_elem_y));
		%curRelPos = (xBot-edge1(1))/(edge1(end)-edge1(1));
		%curRelPos = (xBot-xVertices(1))/(xVertices(end)-xVertices(1));
		%curRelPos = (yBot-yVertices(1))/(yVertices(3)-yVertices(1));
		curRelPos = i/(num_join_patches+1);
		curEdgeDensity = edge1density*curRelPos + edge3density*(1-curRelPos);
		curVertEdge = flipud([(xBot + (xTop - xBot) * curEdgeDensity'), (yBot + (yTop - yBot) * curEdgeDensity')]);
	else
		% last patch, until the end of the rectangle
		if i ~= 1 % don't know why it's needed, but it's needed...
			curVertEdge = edge3V;
		else
			curVertEdge = edge3V;
		end
	end
	
	% create horizontal lines
	vert_C = zeros(size(curVertEdge, 1), 1);
	vert_D = ones(size(curVertEdge, 1), 1);
	vert_F = zeros(size(curVertEdge, 1), 1);
	for j = 1:size(curVertEdge, 1)
		vert_C(j) = (prevVertEdge(j, 2) - curVertEdge(j, 2));
		vert_D(j) = (curVertEdge(j, 1) - prevVertEdge(j, 1));
		vert_F(j) = (prevVertEdge(j, 1)*curVertEdge(j, 2) - curVertEdge(j, 1)*prevVertEdge(j, 2));
	end
	
	% find intersections between EACH horizontal line and EACH vertical line
	% for j = 1:size(edge2, 1)
	for j = (iPrevPosX+1):iPosX
		if j ~= 1 && j ~= size(edge2,1)
			%for k = 1:size(curVertEdge, 1)
			for k = 2:size(curVertEdge, 1)-1
				A = [vert_A(j) vert_B(j); vert_C(k), vert_D(k)];
				B = [-vert_E(j); -vert_F(k)];
				intPoint = A\B;
				if doPlot
					plot(intPoint(1), intPoint(2), 'o');
				end
				numAddedPoints = numAddedPoints + 1;
				vertices = [vertices; intPoint'];
				sectionsNodes(k, j) = numAddedPoints; % = size(vertices, 1);
			end
		end
	end
	
	iPrevPosX = iPosX;
	
	% have to increase by one the number of elements
	% get number of elements along y
	num_elem_y = size(edge1, 1) + i*k_addremove;
	
	if i <= num_join_patches
		
		% curVertEdge_p1 = [(xBot + (xTop - xBot) * linspace(0, 1, num_elem_y)'), (yBot + (yTop - yBot) * linspace(0, 1, num_elem_y)')]; %%%% CONSTANT SPACING
		edge1density = interp1(linspace(0,1,length(edge1)), 1-flipud(edge1), linspace(0, 1, num_elem_y));
		edge3density = interp1(linspace(0,1,length(edge3)), edge3, linspace(0, 1, num_elem_y));
		%curRelPos = (xBot-edge1(1))/(edge1(end)-edge1(1));
		%curRelPos = (xBot-xVertices(1))/(xVertices(end)-xVertices(1));
		%curRelPos = (yBot-yVertices(1))/(yVertices(3)-yVertices(1));
		curRelPos = i/(num_join_patches+1);
		curEdgeDensity = edge1density*curRelPos + edge3density*(1-curRelPos);
		curVertEdge_p1 = [(xBot + (xTop - xBot) * curEdgeDensity'), (yBot + (yTop - yBot) * curEdgeDensity')];

		prevVertEdge = flipud(curVertEdge_p1);
	else
		prevVertEdge = flipud(curVertEdge);
	end
end