function H = rmnode(G, N)
%RMNODE Remove nodes from a graph
%   H = RMNODE(G, nodeID) returns a graph H equivalent to G with nodes
%   specified by the string (or cell array of strings) or numeric IDs
%   nodeID removed from it.  All edges in G incident upon the nodes to be
%   removed are also removed.
%
%   Example:
%       % Create and plot a graph. Remove node 'C', and then plot the new
%       % graph.
%       s = {'A' 'A' 'B' 'C' 'D' 'B' 'C' 'B'};
%       t = {'B' 'C' 'C' 'D' 'A' 'E' 'E' 'D'};
%       G = graph(s,t)
%       plot(G)
%       G = rmnode(G,'C')
%       figure, plot(G)
%
%   See also GRAPH, NUMNODES, ADDNODE, RMEDGE

%   Copyright 2014-2016 The MathWorks, Inc.

ind = findnode(G, N);
ind(ind == 0) = [];
% Determine new node ids.
nodesToKeep = 1:numnodes(G);
nodesToKeep(ind) = [];

if size(G.EdgeProperties, 2) == 0
    N = adjacency(G.Underlying);
    N = N(nodesToKeep, nodesToKeep);
    H = graph(matlab.internal.graph.MLGraph(N));
else
    N = adjacency(G.Underlying, 1:numedges(G));
    N = N(nodesToKeep, nodesToKeep);
    edgeind = nonzeros(tril(N));
    H = graph(matlab.internal.graph.MLGraph(N), ...
              G.EdgeProperties(edgeind, :));
end

H.NodeProperties = G.NodeProperties(nodesToKeep, :);

if nargout < 1
    warning(message('MATLAB:graphfun:rmnode:NoOutput'));
end
