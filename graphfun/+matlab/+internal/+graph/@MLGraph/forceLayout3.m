function [x,y,z,iterations] = forceLayout3(G,x0,y0,z0,iterations)
% FORCELAYOUT   Force-directed node layout
%
%   FOR INTERNAL USE ONLY -- This feature is intentionally undocumented.
%   Its behavior may change, or it may be removed in a future release.
%

% Reference: T. Fruchterman and E. Reingold, "Graph drawing by 
% force-directed placement", Software-Practice & Experience, vol. 21 (11),
% pp. 1129-1164, 1991.

%   Copyright 2016 The MathWorks, Inc.

nn = G.NodeCount;
if ~isnumeric(iterations) || ~isreal(iterations) || ~(iterations >= 0) ...
        || ~isfinite(iterations) || fix(iterations) ~= iterations
    error(message('MATLAB:graphfun:graphbuiltin:InvalidIterations'));
end
if ~isnumeric(x0) || ~isvector(x0) || numel(x0) ~= nn || ~isreal(x0)
    error(message('MATLAB:graphfun:graphbuiltin:InvalidXCoordinates',nn));
end
if any(~isfinite(x0))
    error(message('MATLAB:graphfun:graphbuiltin:NonfiniteXCoordinates'));
end
if ~isnumeric(y0) || ~isvector(y0) || numel(y0) ~= nn || ~isreal(y0)
    error(message('MATLAB:graphfun:graphbuiltin:InvalidYCoordinates',nn));
end
if any(~isfinite(y0))
    error(message('MATLAB:graphfun:graphbuiltin:NonfiniteYCoordinates'));
end
if ~isnumeric(z0) || ~isvector(z0) || numel(z0) ~= nn || ~isreal(z0)
    error(message('MATLAB:graphfun:graphbuiltin:InvalidZCoordinates',nn));
end
if any(~isfinite(z0))
    error(message('MATLAB:graphfun:graphbuiltin:NonfiniteZCoordinates'));
end
x = double(full(x0(:)));
y = double(full(y0(:)));
z = double(full(z0(:)));
if iterations > 0
    comp  = connectedComponents(G);
    ncomp = max(comp);
    if ncomp > 1
        % Apply force-directed layout for each connected component.
        A = adjacency(G);
        compcell = cell(ncomp,1);
        for k = 1:ncomp
            % Extract a component's nodes and edges.
            nk = find(comp == k);
            [sk,tk] = find(tril(A(nk,nk)));
            [x(nk),y(nk),z(nk)] = layoutOneConnComp(x(nk),y(nk),z(nk),sk,tk,iterations);
            compcell{k} = nk;
        end
        % Pack the components together.
        xy = matlab.internal.graph.packLayouts([x y],compcell,comp);
        x = xy(:,1);
        y = xy(:,2);
    else
        edges = G.Edges;
        sources = edges(:,1);
        targets = edges(:,2);
        [x,y,z] = layoutOneConnComp(x,y,z,sources,targets,iterations);
    end
end

%--------------------------------------------------------------------------
function [x,y,z] = layoutOneConnComp(x,y,z,sources,targets,iterations)
% Force-directed layout for a graph with edges defined as (source,target)
% pairs. The node ids in sources and targets range from 1 to length(x).
nn = length(x);
if nn <= 1
    x = zeros(nn,1);
    y = zeros(nn,1);
    z = zeros(nn, 1);
elseif nn == 2
    x = [0; 0];
    y = [0; 1];
    z = [0; 0];
else
    meanside = mean([max(x)-min(x), max(y)-min(y), max(z)-min(z)]);
    meanside(meanside == 0) = 1;
    % Spring constant.
    k = meanside/nn^(1/3);
    % Threshold for nudging nodes that are too close.
    gap = 0.1*meanside/nn^(1/3);
    % Cool down schedule to limit the node movement at each iteration.
    temperature = linspace(0.1*meanside,0,iterations+1);
    % Apply force-directed layout.
    dx = zeros(nn,1);
    dy = zeros(nn,1);
    dz = zeros(nn,1);
    oldstate = rng(0,'twister');
    for i = 1:iterations
        [dx,dy,dz] = applyRepulsiveForce(x,y,z,dx,dy,dz,nn,k,gap);
        [dx,dy,dz] = applyAttractiveForce(x,y,z,dx,dy,dz,nn,sources,targets,k,gap);
        [x,y,z]    = moveNodes(x,y,z,dx,dy,dz,temperature(i));
    end
    rng(oldstate);
    % Normalize coordinates.
    x = x - mean(x);
    y = y - mean(y);
    z = z - mean(z);
    r = max(sqrt(x.^2 + y.^2 + z.^2)); % ~= 0
    r(r == 0) = 1;
    rnew = log(nn)/r;
    x = rnew*x;
    y = rnew*y;
    z = rnew*z;
end

%--------------------------------------------------------------------------
function [dx,dy,dz] = applyRepulsiveForce(x,y,z,dx,dy,dz,nn,k,gap)
% Compute node displacement after applying repulsive forces for each node.
% Re-sets all elements of displacements dx and dy.
gapsq = gap^2;
ksq = k^2;
for i = 1:nn
    deltax = x(i) - x;
    deltay = y(i) - y;
    deltaz = z(i) - z;
    distsq = deltax.^2 + deltay.^2 + deltaz.^2;
    % Nudge nodes that are too close.
    ind = (distsq < gapsq);
    deltaind = randn(sum(ind),3);
    deltax(ind) = deltaind(:,1);
    deltay(ind) = deltaind(:,2);
    deltaz(ind) = deltaind(:,3);
    distsq(ind) = gapsq;
    % Compute repulsive displacement.
    tempdx = ksq*(deltax./distsq);
    tempdy = ksq*(deltay./distsq);
    tempdz = ksq*(deltaz./distsq);
    tempdx(i) = 0;
    tempdy(i) = 0;
    tempdz(i) = 0;
    dx(i) = sum(tempdx);
    dy(i) = sum(tempdy);
    dz(i) = sum(tempdz);
end
%--------------------------------------------------------------------------
function [dx,dy,dz] = applyAttractiveForce(x,y,z,dx,dy,dz,nn,sources,targets,k,gap)
% Compute node displacement after applying attractive forces for each edge.
% Updates displacements dx and dy.
deltaex = x(sources) - x(targets);
deltaey = y(sources) - y(targets);
deltaez = z(sources) - z(targets);
diste = sqrt(deltaex.^2 + deltaey.^2 + deltaez.^2);
tmpdx = deltaex.*(diste/k);
tmpdy = deltaey.*(diste/k);
tmpdz = deltaez.*(diste/k);
% Nudge nodes that are too close.
ind = (diste < gap);
r = 2+rand(sum(ind),1);
tmpind = randn(sum(ind),3);
normtmpind = sqrt(sum(tmpind.^2,2));
tmpind = gap*((r./normtmpind).*tmpind);
tmpdx(ind) = tmpind(:,1);
tmpdy(ind) = tmpind(:,2);
tmpdz(ind) = tmpind(:,3);
% Update displacement.
dx = dx - accumarray(sources,tmpdx,[nn 1]);
dy = dy - accumarray(sources,tmpdy,[nn 1]);
dz = dz - accumarray(sources,tmpdz,[nn 1]);
dx = dx + accumarray(targets,tmpdx,[nn 1]);
dy = dy + accumarray(targets,tmpdy,[nn 1]);
dz = dz + accumarray(targets,tmpdz,[nn 1]);
%--------------------------------------------------------------------------
function [x,y,z] = moveNodes(x,y,z,dx,dy,dz,temperature)
% Move nodes according to displacements and cool down temperature.
d = sqrt(dx.^2 + dy.^2 + dz.^2);
a = min(d,temperature)./d;
a(d == 0) = 0;
x = x + a.*dx;
y = y + a.*dy;
z = z + a.*dz;
