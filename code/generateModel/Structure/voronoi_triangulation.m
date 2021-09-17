function [vxx,vy,c_combined,cells] = voronoi_ext(varargin)
%VORONOI Voronoi diagram.
%   VORONOI(X,Y) plots the Voronoi diagram for the points X,Y. Lines-to-
%   infinity are approximated with an arbitrarily distant endpoint.
%
%   VORONOI(X,Y,OPTIONS) specifies a cell array of strings OPTIONS 
%   that were previously used by Qhull. Qhull-specific OPTIONS are no longer 
%   required and are currently ignored. Support for these options will be 
%   removed in a future release. 
%
%   VORONOI(X,Y,TRI) uses the Delaunay triangulation TRI instead of 
%   computing it internally.
%
%   VORONOI(DT) uses the Delaunay triangulation DT instead of computing it 
%   internally, where DT is a DelaunayTri.
%
%   VORONOI(AX,...) plots into AX instead of GCA.
%
%   H = VORONOI(...,'LineSpec') plots the diagram with color and linestyle
%   specified and returns handles to the line objects created in H.
%
%   VORONOI(DT,'LineSpec',NM) creates segments "to infinity" with a length of NM
%
%   [VX,VY] = VORONOI(...) returns the vertices of the Voronoi edges in VX 
%   and VY. 
%
%   For the topology of the voronoi diagram, i.e. the vertices for
%   each voronoi cell, use DelaunayTri/voronoiDiagram as follows: 
%
%         dt = DelaunayTri(X(:),Y(:))
%         [V,C] = voronoiDiagram(dt)
%
%   See also DelaunayTri, VORONOIN, DELAUNAY, CONVHULL.

%   Copyright 1984-2009 The MathWorks, Inc.
%   $Revision: 1.15.4.16 $  $Date: 2010/11/22 02:46:48 $

[cax,args,nargs] = axescheck(varargin{:});
if nargs < 1
    error(message('MATLAB:narginchk:notEnoughInputs'));
elseif nargs > 4
    error(message('MATLAB:narginchk:tooManyInputs'));
end

if isa(args{1}, 'delaunayTriangulation')       
    dt = args{1};
    if dt.size(1) == 0
        error('MATLAB:voronoi:EmptyTri',...
          'The triangulation is empty.');
    elseif dt.size(2) ~= 3
        error('MATLAB:voronoi:NonPlanarTri',...
          'The triangulation must be composed of triangles.');
    end
    x = dt.Points(:,1);
    y = dt.Points(:,2);
    tri = dt(:,:); 
    if nargs == 1
        ls = '';
		nm = NaN;
	elseif nargs == 2
        ls = args{2};
		nm = NaN;
	else
		ls = args{2};
		nm = args{3};
	end
else
    x = args{1};
    y = args{2};
    if ~isequal(size(x),size(y))
        error('MATLAB:voronoi:InputSizeMismatch',...
              'X and Y must be the same size.');
    end
    if ndims(x) > 2 || ndims(y) > 2
        error('MATLAB:voronoi:HigherDimArray',...
              'X,Y cannot be arrays of dimension greater than two.');
    end
    
    x = x(:);
    y = y(:);
    if nargs == 2,
        tri = delaunay(x,y);
        ls = '';
    else 
        arg3 = args{3};
        if nargs == 3,
            ls = '';
        else
            arg4 = args{4};
            ls = arg4;
        end 
        if isempty(arg3),
            tri = delaunay(x,y);
        elseif ischar(arg3),
            tri = delaunay(x,y); 
            ls = arg3;
        elseif iscellstr(arg3),
             warning('MATLAB:voronoi:DeprecatedOptions',...
            ['VORONOI will not support Qhull-specific options in a future release.\n',...
             'Please remove these options when calling VORONOI.']);
            tri = delaunay(x,y);
        else
            tri = arg3;
        end
    end    
end

if isempty(tri)
    return;
end


% Compute centers of triangles
tr = TriRep(tri,x,y);
c = tr.circumcenters();


% Create matrix T where i and j are endpoints of edge of triangle T(i,j)
n = numel(x);
t = repmat((1:size(tri,1))',1,3);
T = sparse(tri,tri(:,[3 1 2]),t,n,n); 

% i and j are endpoints of internal edge in triangle E(i,j)
E = (T & T').*T; 
% i and j are endpoints of external edge in triangle F(i,j)
F = xor(T, T').*T;

% v and vv are triangles that share an edge
[~,~,v] = find(triu(E));
[~,~,vv] = find(triu(E'));

% Internal edges
vx = [c(v,1) c(vv,1)]';
vy = [c(v,2) c(vv,2)]';

%%% Compute lines-to-infinity
% i and j are endpoints of the edges of triangles in z
[i,j,z] = find(F);
% Counter-clockwise components of lines between endpoints
dx = x(j) - x(i);
dy = y(j) - y(i);

% Calculate scaling factor for length of line-to-infinity
% Distance across range of data
rx = max(x)-min(x); 
ry = max(y)-min(y);
% Distance from vertex to center of data
cx = (max(x)+min(x))/2 - c(z,1); 
cy = (max(y)+min(y))/2 - c(z,2);
% Sum of these two distances
if isnan(nm) % can be supplied as a parameter
	nm = sqrt(rx.*rx + ry.*ry) + sqrt(cx.*cx + cy.*cy);
end
% Compute scaling factor
scale = nm./sqrt((dx.*dx+dy.*dy));
    
% Lines from voronoi vertex to "infinite" endpoint
% We know it's in correct direction because compononents are CCW
ex = [c(z,1) c(z,1)-dy.*scale]';
ey = [c(z,2) c(z,2)+dx.*scale]';

ubound_v = max([v; vv]);
if isempty(ubound_v); ubound_v = 0; end
num_inf_endpoints = length(z);
%v_combined = [v; z];
%vv_combined = [vv; ubound_v + (1:num_inf_endpoints)'];
c_combined = [c; [c(z,1)-dy.*scale, c(z,2)+dx.*scale]];
[iF,jF,~] = find(F);
numbering_inf_endpoints = zeros(size(F));
for i = 1:num_inf_endpoints
	numbering_inf_endpoints(iF(i), jF(i)) = ubound_v + i;%1:num_inf_endpoints;
end

cells = {};
E_up = E;%triu(E);
E_dn = E;%tril(E);
cur_point_infty = 0;
for i = 1:n
	[~,~,segF_from] = find(E_up(i,:));
	[~,~,segF_to] = find(E_dn(:,i)');
	[~,segI_from_j1,segI_from_z1] = find(F(i,:));
	[segI_from_i2,~,segI_from_z2] = find(F(:,i));
	segI_from = [segI_from_z1(:) segI_from_z2(:)];
	segI_to1 = numbering_inf_endpoints(i, segI_from_j1);
	segI_to2 = numbering_inf_endpoints(segI_from_i2, i);
	% segI_to = ubound_v + cur_point_infty + (1:length(segI_from));
	segI_to = [segI_to1(:) segI_to2(:)];
	cur_point_infty = cur_point_infty + length(segI_from);
	if size(segI_from,1) > 0
		cells{i} = [segF_from segI_from; segF_to segI_to];
	else
		cells{i} = [segF_from; segF_to];
	end
end

%{
% TEST:Plot
figure
plot(x, y, '.b')
hold on
for i = 1:length(cells)
	for j = 1:size(cells{i}, 2)
		line([c_combined(cells{i}(1,j),1) c_combined(cells{i}(2,j),1)], [c_combined(cells{i}(1,j),2) c_combined(cells{i}(2,j),2)])
	end
end
%}

% Combine with internal edges
vx = [vx ex];
vy = [vy ey];

if nargout<2
    % Plot diagram
    if isempty(cax)
        % If no current axes, create one
        cax = gca;
    end
    if isempty(ls)
        % Default linespec
        ls = '-';
    end
    [l,c,mp,msg] = colstyle(ls); error(msg) % Extract from linespec
    if isempty(mp)
        % Default markers at points        
        mp = '.';
    end
     if isempty(l)
        % Default linestyle
        l = get(ancestor(cax,'figure'),'DefaultAxesLineStyleOrder'); 
    end
    if isempty(c), 
        % Default color        
        co = get(ancestor(cax,'figure'),'DefaultAxesColorOrder');
        c = co(1,:);
    end
    % Plot points
    h1 = plot(x,y,'marker',mp,'color',c,'linestyle','none','parent',cax);
    % Plot voronoi lines
    h2 = line(vx,vy,'color',c,'linestyle',l,'parent',cax,...
        'yliminclude','off','xliminclude','off');
    if nargout==1, vxx = [h1; h2]; end % Return handles
else
    vxx = vx; % Don't plot, just return vertices
end