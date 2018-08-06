function h = flipedge(g, s, t)
%FLIPEDGE Flip edge directions
%   H = FLIPEDGE(G) returns digraph H, which contains the same edges as
%   digraph G, with reversed directions. H has the same node and edge
%   properties as G.
%
%   H = FLIPEDGE(G,S,T) reverses the direction of edges specified by
%   pairs of node IDs S and T.
%
%   H = FLIPEDGE(G,IND) reverses the direction of edges specified by the
%   edge indices IND.
%
%  See also: DIGRAPH

%   Copyright 2016 The MathWorks, Inc.

if nargin <= 1
    if size(g.EdgeProperties, 2) == 0
        h = g;
        h.Underlying = matlab.internal.graph.MLDigraph(adjacency(g.Underlying, 'transp'));
    else
        Hadjacency = adjacency(g.Underlying, 1:numedges(g));
        h = g;
        h.Underlying = matlab.internal.graph.MLDigraph(Hadjacency, 'transp');
        
        % Permute edge properties
        p = nonzeros(Hadjacency);
        h.EdgeProperties = h.EdgeProperties(p, :);
    end
else
    % Determine the indices of the edges to be removed.
    if nargin == 2
        edgeind = s;
        % Reuse input checking in findedge
        [~, ~] = findedge(g, edgeind);
        useNodeNames = hasNodeNames(g);
    else
        edgeind = findedge(g, s, t);
        
        useNodeNames = hasNodeNames(g) && ~isnumeric(s);
        
        % Error if edge doesn't exist
        if any(edgeind == 0)
            ind = find(edgeind == 0, 1);
            if ~useNodeNames
                error(message('MATLAB:graphfun:flipedge:InvalidEdge', s(ind), t(ind)));
            else
                s = findnode(g, s);
                t = findnode(g, t);
                nodeNames = g.NodeProperties.Name;
                error(message('MATLAB:graphfun:flipedge:InvalidEdge', nodeNames{s(ind)}, nodeNames{t(ind)}));
            end
        end
    end
    
    % Check to make sure that there are no duplicates specified.
    if any(accumarray(full(edgeind(:)), 1) > 1)
        error(message('MATLAB:graphfun:flipedge:DuplicateEdge'));
    end
    
    nn = numnodes(g);
    ed = g.Underlying.Edges;
    ed(edgeind, :) = fliplr(ed(edgeind, :));
    Hadjacency = sparse(ed(:, 2), ed(:, 1), 1:numedges(g), nn, nn);
    
    % Check if only one of the edges in a two-cycle is flipped
    if nnz(Hadjacency) ~= numedges(g)
        Hduplicates = sparse(ed(:, 2), ed(:, 1), 1, nn, nn) > 1;

        [s, t] = find(Hduplicates, 1);
        if useNodeNames
           s = g.NodeProperties.Name{s};
           t = g.NodeProperties.Name{t};
        end
        error(message('MATLAB:graphfun:flipedge:TwoCycle', s, t, t, s));
    end
    
    h = g;
    h.Underlying = matlab.internal.graph.MLDigraph(Hadjacency, 'transp');
    
    % Permute edge properties if they exist
    if size(h.EdgeProperties, 2) > 0
        p = nonzeros(Hadjacency);
        h.EdgeProperties = h.EdgeProperties(p, :);
    end
    
end

if nargout < 1
    warning(message('MATLAB:graphfun:flipedge:NoOutput'));
end