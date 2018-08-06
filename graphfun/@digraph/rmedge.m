function H = rmedge(G, s, t)
%RMEDGE Remove edges from a digraph
%   H = RMEDGE(G, s, t) returns a new digraph H equivalent to G but with
%   edges specifed by pairs of node IDs s and t removed.  s and t must both
%   be strings or cell arrays of strings specifying names of nodes or
%   numeric node IDs.
%
%   H = RMEDGE(G, ind) returns a new digraph H equivalent to G but with
%   edges specified by the edge index ind removed.
%
%   RMEDGE never removes nodes.
%
%   Example:
%       % Create and plot a digraph. Remove two edges, then plot the new 
%       % digraph.
%       s = {'A' 'A' 'B' 'C' 'D' 'B' 'C' 'B'};
%       t = {'B' 'C' 'C' 'D' 'A' 'E' 'E' 'D'};
%       G = digraph(s,t)
%       plot(G)
%       G = rmedge(G,{'A' 'B'},{'C' 'D'})
%       figure, plot(G)
%
%   Example:
%       % Create a digraph and view the edge list. Remove edge 3, then view
%       % the new edge list.
%       s = {'BOS' 'NYC' 'NYC' 'NYC' 'LAX'};
%       t = {'NYC' 'LAX' 'DEN' 'LAS' 'DCA'};
%       G = digraph(s,t);
%       G.Edges
%       G = rmedge(G,3);
%       G.Edges
%
%   See also DIGRAPH, NUMEDGES, ADDEDGE, RMNODE

%   Copyright 2014-2015 The MathWorks, Inc.

H = G;
% Determine the indices of the edges to be removed.
if nargin == 2
  ind = s;
else
  ind = findedge(G, s, t);
  ind(ind == 0) = [];
end
% Convert the edge indices to pairs of Node IDs.
[s, t] = findedge(G, ind);
% Remove edges from the graph.
nn = numnodes(H);
A = adjacency(H.Underlying, 'transp');
A(sub2ind([nn, nn], t, s)) = 0;
H.Underlying = matlab.internal.graph.MLDigraph(A, 'transp');
% Remove corresponding edge properties.
H.EdgeProperties(ind, :) = [];

if nargout < 1
    warning(message('MATLAB:graphfun:rmedge:NoOutput'));
end
