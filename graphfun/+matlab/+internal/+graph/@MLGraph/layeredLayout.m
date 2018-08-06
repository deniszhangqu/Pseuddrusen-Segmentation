function [nodesCoord, edgesCoord] = layeredLayout(g, sources, sinks, asgnLay)
% LAYEREDLAYOUT   Compute layered layout
%
%   FOR INTERNAL USE ONLY -- This feature is intentionally undocumented.
%   Its behavior may change, or it may be removed in a future release.
%

%   Copyright 2015 The MathWorks, Inc.

% Construct a directed acyclic graph from undirected graph g
% and feed that into MLDigraph/layeredLayout

% Make acyclic (replace every edge by an edge in one direction, 
% such that the resulting graph is acyclic

if g.NodeCount == 0
    gdag = matlab.internal.graph.MLDigraph;
else
    gdag = makeAcyclic(g, sources, sinks);
end

[nodesCoord, edgesCoord] = layeredLayout(gdag, sources, sinks, asgnLay);

p = findedge(g, gdag.Edges(:, 1), gdag.Edges(:, 2));

edgesCoord(p) = edgesCoord;

function gdag = makeAcyclic(g, sources, sinks)
% Orient edges of g to create acyclic graph gdag. 

% Remove self-loops and revert some edges of g to create acyclic graph
% gdag. isreverted has length numedges(g) and is true if the edge is
% reverted in gdag.

nn = g.NodeCount;
M = adjacency(g);
Mnew = M;

% Remove self-loops
Mnew(1:nn+1:end) = 0;

if ~isempty(sources)
    % Remove all incoming edges of source nodes
    Mnew(:, sources) = 0;
end

if ~isempty(sinks)
    % Remove all outgoing edges of sink nodes
    Mnew(sinks, :) = 0;
end
 
if isempty(sources) && isempty(sinks)
% Construct a sinks or sources array (only for the purpose of removing cycles): 
    sources = 1;
end

% If no sources are set, revert all edges and use sinks instead
swapsrcsink = isempty(sources);
if swapsrcsink
    % Swap sources and sinks, revert after the following
    sources = sinks;
    Mnew = Mnew.';
end

% Make gdag acyclic by reverting the back-edges found by dfsearch starting
% at nodes sources
edgeToDiscovered = false(1, 6);
edgeToDiscovered(4) = true;

Mhelper = Mnew;
Mhelper(nn+1, nn+1) = 0;    % add a helper node
Mhelper(nn+1, sources) = 1; % with edges to all sources

h = matlab.internal.graph.MLDigraph(Mhelper);

restart = true; % needed if there are several weak components
formattable = false;
revedges = depthFirstSearch(h, nn+1, edgeToDiscovered, restart, formattable);

Mnew(sub2ind([nn, nn], revedges(:, 1), revedges(:, 2))) = 0;
Mnew(sub2ind([nn, nn], revedges(:, 2), revedges(:, 1))) = 1;

if swapsrcsink
    Mnew = Mnew.';
end

%gdag = matlab.internal.graph.MLDigraph(Mnew);
%assert(dfsTopologicalSort(gdag));

Mnew(1:nn+1:end) = diag(M);
Mnew(sources, sources) = triu(M(sources, sources));
Mnew(sinks, sinks) = triu(M(sinks, sinks));

gdag = matlab.internal.graph.MLDigraph(Mnew);
