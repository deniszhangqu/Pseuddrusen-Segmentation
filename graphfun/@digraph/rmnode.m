function H = rmnode(G, N)
%RMNODE Remove nodes from a digraph
%   H = RMNODE(G, NodeID) returns a digraph H equivalent to G with nodes
%   specified by the string (or cell array of strings) or numeric IDs 
%   NodeID removed from it.  All edges in G incident upon the nodes to be 
%   removed are also removed.
%
%   Example:
%       % Create and plot a digraph. Remove node 'C', and then plot the new
%       % digraph.
%       s = {'A' 'A' 'B' 'C' 'D' 'B' 'C' 'B'};
%       t = {'B' 'C' 'C' 'D' 'A' 'E' 'E' 'D'};
%       G = digraph(s,t)
%       plot(G)
%       G = rmnode(G,'C')
%       figure, plot(G)
%
%   See also DIGRAPH, NUMNODES, ADDNODE, RMEDGE

%   Copyright 2014-2016 The MathWorks, Inc.

ind = findnode(G, N);
ind(ind == 0) = [];
% Determine new node ids.
nodesToKeep = 1:numnodes(G);
nodesToKeep(ind) = [];

if size(G.EdgeProperties, 2) == 0
    N = adjacency(G.Underlying, 'transp');
    N = N(nodesToKeep, nodesToKeep);
    H = digraph(matlab.internal.graph.MLDigraph(N, 'transp'));
else
    N = adjacency(G.Underlying, 1:numedges(G), 'transp');
    N = N(nodesToKeep, nodesToKeep);
    edgeind = nonzeros(N);
    H = digraph(matlab.internal.graph.MLDigraph(N, 'transp'), ...
                G.EdgeProperties(edgeind, :));
end

H.NodeProperties = G.NodeProperties(nodesToKeep, :);

if nargout < 1
    warning(message('MATLAB:graphfun:rmnode:NoOutput'));
end
