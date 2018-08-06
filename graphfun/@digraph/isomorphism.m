function perm = isomorphism(G1, G2, varargin)
%ISOMORPHISM Compute an isomorphism between G and G2
%   P = ISOMORPHISM(G,G2) returns an isomorphism between G and G2, if one
%   exists. If no isomorphism exists, then P is an empty array. If an
%   isomorphism exists, then REORDERNODES(G2,P) has the same structure as G.
%
%   P = ISOMORPHISM(..., 'NodeVariables', PROPLIST) requires the
%   isomorphism P to preserve all node properties listed in PROPLIST.
%   PROPLIST must be a string vector, a character vector or a cell array of
%   character vectors.
%
%   P = ISOMORPHISM(..., 'EdgeVariables', PROPLIST) requires the
%   isomorphism P to preserve all edge properties listed in PROPLIST.
%   PROPLIST must be a string vector, a character vector or a cell array of
%   character vectors. If PROPLIST contains 'EndNodes', it is ignored.
%
%   See also DIGRAPH, ISISOMORPHIC, GRAPH/ISOMORPHISM

%   Copyright 2016 The MathWorks, Inc.

if ~isa(G1, 'digraph') || ~isa(G2, 'digraph')
    error(message('MATLAB:graphfun:isomorphism:ClassMismatch'));
end

[nodeCell, edgeCell] = parseInputs(varargin{:});

% Find indices of node and edge variables in the properties tables. Error
% if variable not in properties table, or if datatypes in G1 and G2 do not
% match.
[nodeVarInd1, nodeVarInd2] = extractVarInd(nodeCell, G1.NodeProperties, G2.NodeProperties, 'Node');
[edgeVarInd1, edgeVarInd2] = extractVarInd(edgeCell, G1.EdgeProperties, G2.EdgeProperties, 'Edge');

% Early return if no match is possible:
if numnodes(G1) ~= numnodes(G2) || numedges(G1) ~= numedges(G2)
    perm = [];
    return;
elseif numnodes(G1) == 0
    perm = zeros(1, 0);
    return;
end

% Compute unique indices for each edge and node in G1 and G2, to be used
% later. Error if unique not supported for input datatype
[nodeCat, nodeMatchPossible] = extractInd(nodeVarInd1, nodeVarInd2, G1.NodeProperties, G2.NodeProperties);
[edgeCat, edgeMatchPossible] = extractInd(edgeVarInd1, edgeVarInd2, G1.EdgeProperties, G2.EdgeProperties);

if ~nodeMatchPossible || ~edgeMatchPossible
    perm = [];
    return;
end

% Add node degree and size of containing connected component to nodeCat
[nodeCat, matchPossible] = refineNodeCategories(G1.Underlying, G2.Underlying, nodeCat, edgeCat);

if matchPossible
    perm = computeIsomorphism(G1, G2, nodeCat, edgeCat);
else
    perm = [];
end


function [nodeCell, edgeCell] = parseInputs(varargin)

nodeCell = {};
edgeCell = {};

if numel(varargin) == 0
    return;
end

for ii=1:2:numel(varargin)
    name = varargin{ii};
    if ~digraph.isvalidoption(name)
        error(message('MATLAB:graphfun:isomorphism:ParseFlags'));
    end
    
    if digraph.partialMatch(name, "NodeVariables")
        edgeMatch = false;
        key = 'NodeVariables';
    elseif digraph.partialMatch(name, "EdgeVariables")
        edgeMatch = true;
        key = 'EdgeVariables';
    else
        error(message('MATLAB:graphfun:isomorphism:ParseFlags'));
    end
    
    if ii+1 > numel(varargin)
        error(message('MATLAB:graphfun:isomorphism:KeyWithoutValue', key));
    end
    value = varargin{ii+1};
    if ischar(value) || iscellstr(value)
        value = string(value);
    elseif ~isstring(value)
        error(message('MATLAB:graphfun:isomorphism:ParseVar', key));
    end
    
    if edgeMatch
        % Remove any instances of 'EndNodes' in the array
        value(value == "EndNodes") = [];
        edgeCell = value;
    else
        nodeCell = value;
    end
    
end

function [varInd1, varInd2] = extractVarInd(varnameString, t1, t2, nodeEdgeStr)

if isempty(varnameString)
    varInd1 = [];
    varInd2 = [];
    return;
end

% Index into node variables - error if variable doesn't exist
[~, varInd1] = ismember(varnameString, t1.Properties.VariableNames);
[~, varInd2] = ismember(varnameString, t2.Properties.VariableNames);

if any(varInd1==0)
    % Convert to char because message does not support string at the moment
    varName = char(varnameString(find(varInd1==0, 1)));
    error(message(['MATLAB:graphfun:isomorphism:' nodeEdgeStr 'VarNotInFirstGraph'], varName));
elseif any(varInd2 == 0)
    varName = char(varnameString(find(varInd2==0, 1)));
    error(message(['MATLAB:graphfun:isomorphism:' nodeEdgeStr 'VarNotInSecondGraph'], varName));
end

% Check that all variables have identical type in t1 and t2 - not automatic
% because cat will silently convert between some types.
for ii=1:numel(varInd1)
    var1 = t1{:, varInd1(ii)};
    var2 = t2{:, varInd2(ii)};
        
    if ~isequal(class(var1), class(var2))
        error(message(['MATLAB:graphfun:isomorphism:' nodeEdgeStr 'VarTypeMismatch'], ...
            char(varnameString(ii)), class(var1), class(var2)));
    end
    
    sz1 = size(var1);
    sz2 = size(var2);
    if ~isequal(sz1(2:end), sz2(2:end))
        error(message(['MATLAB:graphfun:isomorphism:' nodeEdgeStr 'VarSizeMismatch'], char(varnameString(ii))));
    end
    
end


function [ind, matchPossible] = extractInd(varInd1, varInd2, t1, t2)
% Extract an mx2 matrix of indices from the node or edge variables
% specified. ind(:, 1) represents the variables for G1, ind(:, 2)
% represents the variables for G2.

if isempty(varInd1)
    ind = [];
    matchPossible = true;
else
    propTable = [t1(:, varInd1); t2(:, varInd2)];
    [~, ~, ind] = unique(propTable);
    ind = uint64(reshape(ind, [size(t1, 1) 2]));
    
    % Check if node index distributions are equal; otherwise, there is no
    % possible permutation that matches the node properties
    matchPossible = isequal(accumarray(ind(:, 1), 1), accumarray(ind(:, 2), 1));
end


function [nodeCat, nodeMatchPossible] = refineNodeCategories(mlG1, mlG2, nodeCat, edgeCat)
% Refine the matrix nodeCat based on the structures of graphs G1 and G2

% First step: add additional categories:

% 1) in- and out-degree of each node
outdegrees = [outdegree(mlG1), outdegree(mlG2)];
indegrees = [indegree(mlG1), indegree(mlG2)];

% 2) size of strongly connected component containing each node
bins1 = connectedComponents(mlG1);
bins2 = connectedComponents(mlG2);
strongCompSizes = [compBinSize(bins1), compBinSize(bins2)];

% 3) size of strongly connected component containing each node
bins1 = weakConnectedComponents(mlG1);
bins2 = weakConnectedComponents(mlG2);
weakCompSizes = [compBinSize(bins1), compBinSize(bins2)];

% Compute updated matrix nodeCat
[~, ~, nodeCat] = unique([nodeCat(:) outdegrees(:) indegrees(:) strongCompSizes(:) weakCompSizes(:)], 'rows');
nodeCat = reshape(nodeCat, numnodes(mlG1), 2);

% Test if a match between G1 and G2 is possible based on nodeCat:
nodeMatchPossible = isequal(accumarray(nodeCat(:, 1), 1), accumarray(nodeCat(:, 2), 1));

if nodeMatchPossible
    % Use successive comparison of the categories of all neighboring nodes
    % to refine nodeCat as much as possible.
    [nodeMatchPossible, nodeCat] = refineNodeCat(mlG1, mlG2, uint64(nodeCat-1), edgeCat);
end


function binSize = compBinSize(bins)
% given the bin number for each node, return the bin size (this is
% independent of permutations of the graph, the bin number itself is not)

binSize = accumarray(bins(:), 1);
binSize = binSize(bins);
binSize = binSize(:);


function perm = computeIsomorphism(G1, G2, nodeCat, edgeCat)
% Compute isomorphism of G1 and G2, with numeric node and edge categories.

% Reorder the nodes in G1 such that nodes with low multiplicities
% in nodeCat are at the start:
multip = accumarray(nodeCat(:, 1), 1);
nodeMult = multip(nodeCat(:, 1));
[~, ind] = sort(nodeMult);

% Compute re-ordered version of G1 made undirected
N = adjacency(G1.Underlying, 1:numedges(G1), 'transp');
Nundir = N | N';
modG1 = matlab.internal.graph.MLDigraph(Nundir(ind, ind), 'transp');

% Update order ind to be depth-first search discovery order of modG1:
labelVector = eye(1, 6, 'logical');
ind2 = breadthFirstSearch(modG1, 1, labelVector, true, false);
ind = ind(ind2);

% Construct graph modG1 with nodes reordered according to ind:
N = N(ind, ind);
modG1 = matlab.internal.graph.MLDigraph(N, 'transp');
nodeCat(:, 1) = nodeCat(ind, 1);

if size(edgeCat, 2) == 0
    edgeCat = zeros(numedges(G1), 2, 'uint64');
else
    edgeind = nonzeros(N);
    edgeCat(:, 1) = edgeCat(edgeind, 1);
end

% Compute isomorphism. This algorithms matches nodes 1, 2, ..., n of modG1
% against G2, which is why the permutation ind is important for
% performance.
[isIso, perm] = isomorphism(modG1, G2.Underlying, nodeCat, edgeCat);

if isIso
    perm(ind) = perm;
end
