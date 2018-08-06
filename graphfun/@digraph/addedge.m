function H = addedge(G, s, t, weights)
%ADDEDGE Add edges to a digraph
%   H = ADDEDGE(G,s,t) returns digraph H that is equivalent to G but with
%   edges specified by s and t added to it.  s and t must both refer to
%   string node names or numeric node indices.  If a node specified by s or
%   t is not present in the digraph G, that node is added as well.  s and t
%   must specify edges not already present in G.
%
%   H = ADDEDGE(G,s,t,w), where G is a weighted graph, adds edges with
%   corresponding edge weights defined by w.  w must be numeric.
%
%   H = ADDEDGE(G,EdgeTable) adds edges with attributes specified by
%   the table EdgeTable.  EdgeTable must be able to be concatenated with
%   G.Edges.
%
%   H = ADDEDGE(G,s,t,EdgeTable) adds edges with attributes specified by
%   the table EdgeTable.  EdgeTable must not contain a variable EndNodes,
%   and must be able to be concatenated with G.Edges(:, 2:end).
%
%   Example:
%       % Construct a digraph with three edges, then add two new edges.
%       G = digraph([1 2 3],[2 3 4])
%       G.Edges
%       G = addedge(G,[2 1],[4 6])
%       G.Edges
%
%   See also DIGRAPH, NUMEDGES, RMEDGE, ADDNODE

%   Copyright 2014-2017 The MathWorks, Inc.

tableOK = false;
if istable(s)
    if nargin > 2
        error(message('MATLAB:graphfun:addedge:TableMaxRHS'));
    end
    if size(s,2) < 1
        error(message('MATLAB:graphfun:addedge:TableSize'));
    end
    if ~strcmp('EndNodes', s.Properties.VariableNames{1})
        error(message('MATLAB:graphfun:addedge:TableFirstVar'));
    end
    if size(s.EndNodes,2) ~= 2 || ~(isnumeric(s.EndNodes) || ...
        iscellstr(s.EndNodes))
        error(message('MATLAB:graphfun:addedge:BadEndNodes'));
    end
    % Extract into s, t, w.
    t = s.EndNodes(:,2);
    weights = s(:,2:end);
    s = s.EndNodes(:,1);
    tableOK = true;
elseif nargin >= 4 && istable(weights)
    if any(strcmp('EndNodes', weights.Properties.VariableNames))
        error(message('MATLAB:graphfun:addedge:DuplicateEndNodes'));
    end
    tableOK = true;
elseif nargin < 4 && hasEdgeWeights(G)
    error(message('MATLAB:graphfun:addedge:SpecifyWeight'));
end

inputsAreStrings = digraph.isvalidstring(s) && digraph.isvalidstring(t);
if ~inputsAreStrings && ~(isnumeric(s) && isnumeric(t))
    error(message('MATLAB:graphfun:addedge:InconsistentNodeNames'));
end

% Add any nodes that are not present.
if inputsAreStrings
    if ~iscellstr(s), s = {s}; end
    if ~iscellstr(t), t = {t}; end
    if numel(t) == 1 && numel(s) >= 1
        refdNodes = [t; s(:)]; refdNodes([1 2]) = refdNodes([2 1]);
    elseif numel(s) == numel(t)
        refdNodes = [s(:).'; t(:).'];
    else
        refdNodes = [s(:); t(:)];
    end
    refdNodes = unique(refdNodes(:), 'stable');
    if hasNodeNames(G)
        % Lookup node names and add any that we might need.
        inds = ismember(refdNodes, G.NodeProperties.Name);
        newNodes = refdNodes(~inds);
    else
        newNodes = refdNodes;
    end
    H = addnode(G, newNodes);
    % Replace s/t with numeric indices.
    [~, s] = ismember(s(:), H.NodeProperties.Name);
    [~, t] = ismember(t(:), H.NodeProperties.Name);
else
    s = double(s(:));
    t = double(t(:));
    ms = validateNodeIDs(s);
    mt = validateNodeIDs(t);
    N = max(ms, mt);
    if N > numnodes(G)
        H = addnode(G, N-numnodes(G));
    else
        H = G;
    end
end

% Compute adjacency matrix of edges to add
nn = numnodes(H);
A = sparse(t, s, 1, nn, nn);

% Check to make sure that there are no duplicates specified.
if nnz(A > 1) ~= 0
    error(message('MATLAB:graphfun:addedge:DuplicateEdge'));
end

% Check to see if we are attempting to add any edges that are already
% present and error appropriately.
Hadjacency = adjacency(H.Underlying, 'transp');

if nnz(A ~= 0 & Hadjacency ~= 0) ~= 0
    [t, s] = find(A ~= 0 & Hadjacency ~= 0, 1);
    if ~inputsAreStrings
        error(message('MATLAB:graphfun:addedge:EdgeExists', s, t));
    else
        nodeNames = H.NodeProperties.Name;
        error(message('MATLAB:graphfun:addedge:EdgeExists', nodeNames{s}, nodeNames{t}));
    end
end

% Finally, add the edges.
A = A + Hadjacency;
H.Underlying = matlab.internal.graph.MLDigraph(A, 'transp');

p = findedge(H.Underlying, s, t);

q = true(numedges(H), 1);
q(p) = false;

EdgePropTable = H.EdgeProperties;
if nargin == 4 || tableOK
    if isnumeric(weights)
        if numedges(G) > 0
            if ~hasEdgeWeights(G)
                error(message('MATLAB:graphfun:addedge:NoWeights'));
            end
            EdgePropTable = expandTable(G.EdgeProperties, q);
        end
        if ~isscalar(weights) && numel(p) ~= numel(weights)
            error(message('MATLAB:graphfun:addedge:NumWeightsMismatch'));
        end
        EdgePropTable{p,'Weight'} = weights(:);
    elseif tableOK && istable(weights)
        if numedges(G) > 0
            if size(G.EdgeProperties, 2) ~= size(weights, 2)
                error(message('MATLAB:table:VarDimensionMismatch'));
            end
            if ~isequal(G.EdgeProperties.Properties.VariableNames, ...
                    weights.Properties.VariableNames)
                error(message('MATLAB:table:VarDimensionMismatch'));
            end
            EdgePropTable = expandTable(G.EdgeProperties, q);
        end
        EdgePropTable(p,:) = weights;
    else
        error(message('MATLAB:graphfun:addedge:FourthInput'));
    end
else
    if size(G.EdgeProperties, 2) == 0
        EdgePropTable = table.empty(numedges(H), 0);
    else
        EdgePropTable = expandTable(G.EdgeProperties, q);
    end
end
H.EdgeProperties = EdgePropTable;
if nargout < 1
    warning(message('MATLAB:graphfun:addedge:NoOutput'));
end

function m = validateNodeIDs(ids)
if ~isreal(ids) || any(fix(ids)~=ids) || any(ids < 1)
    error(message('MATLAB:graphfun:addedge:InvalidNodeID'));
end
m = max(ids(:));

function tnew = expandTable(t, q)
% t is a table, q a logical array with nnz(q) == size(t, 1).
% Return value is a table tnew with numel(q) rows, where tnew(q) = t
% and all other rows of tnew are the result of table expansion.

tnew = t([], :);
tnew(numel(q)+1, :) = t(1, :);
tnew(numel(q)+1, :) = [];

tnew(q, :) = t;
