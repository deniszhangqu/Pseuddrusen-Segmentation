function [G, EdgeProps, NodeProps] = ...
    constructFromEdgeList(underlyingCtor, errTag, s, t, varargin)
% CONSTRUCTFROMEDGELIST Construct graph/digraph from edge list representation.

% Copyright 2015-2016 The MathWorks, Inc.

if nargin > 7
    error(message('MATLAB:maxrhs'));
end
% Peel off the last arg, and check it for OmitSelfLoops.
omitLoops = false;
if numel(varargin) > 0
    flag = varargin{end};
    if (ischar(flag) && isrow(flag)) || (isstring(flag) && isscalar(flag))
        omitLoops = startsWith("OmitSelfLoops", flag, 'IgnoreCase', true) && strlength(flag) > 0;
        if ~omitLoops
            error(message(['MATLAB:graphfun:' errTag ':InvalidFlag']));
        end
        varargin(end) = [];
    elseif numel(varargin) == 3
        error(message(['MATLAB:graphfun:' errTag ':InvalidFlag']));
    end
end
% Discover and set NodeProps.
NodeProps = [];
% NodeStuff is the input, if present.  It may be
%  * A collection of Node Names.
%  * A numeric scalar indicating the number of nodes.
%  * A table containing node properties.
NodeStuff = {};
totalNodes = [];
userSetNumNodes = false;
explicitNodeNames = false;
if numel(varargin) > 1
    NodeStuff = varargin{2};
    if isvalidstring(NodeStuff)
        if ~iscellstr(NodeStuff), NodeStuff = {NodeStuff}; end
        NodeStuff = validateName(NodeStuff(:), errTag);
        totalNodes = numel(NodeStuff);
        explicitNodeNames = true;
        NodeProps = table(NodeStuff, 'VariableNames', {'Name'});
    elseif isnumeric(NodeStuff) && isscalar(NodeStuff)
        totalNodes = NodeStuff;
        if ~isreal(totalNodes) || ~isfinite(totalNodes) ...
                || fix(totalNodes) ~= totalNodes
            error(message(['MATLAB:graphfun:' errTag ':InvalidNumNodesProps']));
        end
        userSetNumNodes = true;
        NodeStuff = {};
    elseif istable(NodeStuff)
        totalNodes = size(NodeStuff,1);
        % Validate Nodes Table.
        if any(strcmp('Name', NodeStuff.Properties.VariableNames))
            if ~iscolumn(NodeStuff.Name)
                error(message(['MATLAB:graphfun:' errTag ':InvalidNames']));
            end
            validateName(NodeStuff.Name, errTag);
            explicitNodeNames = true;
        end
        NodeProps = NodeStuff;
    else
        error(message(['MATLAB:graphfun:' errTag ':EdgeListFourthArg']));
    end
end
if isvalidstring(s) && isvalidstring(t)
    if userSetNumNodes
        error(message(['MATLAB:graphfun:' errTag ':EdgeListNumNodes']));
    end
    if ~iscellstr(s), s = {s}; end
    if ~iscellstr(t), t = {t}; end
    if explicitNodeNames
        [present, s] = ismember(s(:), NodeProps.Name);
        if ~all(present)
            error(message(['MATLAB:graphfun:' errTag ':InvalidNodeRefd']));
        end
        [present, t] = ismember(t(:), NodeProps.Name);
        if ~all(present)
            error(message(['MATLAB:graphfun:' errTag ':InvalidNodeRefd']));
        end
    else
        if istable(NodeStuff)
            error(message(['MATLAB:graphfun:' errTag ':NodesTableNeedsName']));
        end
        if numel(t) == 1 && numel(s) >= 1
            Name = [t; s(:)]; Name(1:2) = Name([2 1]);
        elseif numel(s) == numel(t)
            Name = [s(:).'; t(:).'];
        else
            Name = [s(:); t(:)];
        end 
        Name = unique(Name(:), 'stable');
        [~, s] = ismember(s(:), Name);
        [~, t] = ismember(t(:), Name);
        NodeProps = table(Name, 'VariableNames', {'Name'});
        totalNodes = numel(Name);
    end
elseif ~(isnumeric(s) && isnumeric(t))
    error(message(['MATLAB:graphfun:' errTag ':InvalidEdges']));
end
implicitTotal = max([s(:);t(:)]);
if isempty(totalNodes)
    totalNodes = implicitTotal;
elseif totalNodes < implicitTotal
    if explicitNodeNames
        error(message(['MATLAB:graphfun:' errTag ':InvalidNumNodeNames'], implicitTotal));
    elseif istable(NodeStuff)
        error(message(['MATLAB:graphfun:' errTag ':InvalidNumNodesTable'], implicitTotal));
    else
        error(message(['MATLAB:graphfun:' errTag ':InvalidNumNodes'], implicitTotal));
    end
end
% Need the following for when we cope with weights below...
specifiedEdges = max(numel(s), numel(t));
if omitLoops
    omittedRows = (s == t);
    if isscalar(s) && ~isscalar(t)
        t(omittedRows) = [];
    elseif isscalar(t) && ~isscalar(s)
        s(omittedRows) = [];
    else
        s(omittedRows) = [];
        t(omittedRows) = [];
    end
end
G = underlyingCtor(double(s), double(t), totalNodes);
% Define the table if we haven't already.
if ~istable(NodeProps)
    NodeProps = table.empty(G.NodeCount, 0);
end
% Set Edge properties.
EdgeProps = table.empty(G.EdgeCount,0);
if numel(varargin) > 0
    w = varargin{1};
    ignoreWeights = false;
    if ~isfloat(w)
        if ~istable(w)
            error(message(['MATLAB:graphfun:' errTag ':InvalidWeights']));
        end
    else
        % Look for [] meaning ignore weights.
        ignoreWeights = isequal(w,[]);
        w = w(:);
    end
    if ~ignoreWeights
        if (size(w,1) ~= specifiedEdges) && ~(isscalar(w) && isnumeric(w))
            error(message(['MATLAB:graphfun:' errTag ':InvalidSizeWeight']));
        end
        if omitLoops && ~isscalar(w)
            w(omittedRows,:) = [];
        end
        ind = findedge(G, s, t);
        if isfloat(w)
            EdgeProps.Weight(ind,1) = w;
        else
            EdgeProps = table;
            EdgeProps(ind,:) = w;
        end
    end
end

function v = isvalidstring(s)
%v = matlab.internal.datatypes.isCharStrings(s, false, false);
v = iscellstr(s) || (ischar(s) && isrow(s));

function Name = validateName(Name, errTag)
if ~matlab.internal.datatypes.isCharStrings(Name, false, false)
    error(message(['MATLAB:graphfun:' errTag ':InvalidNames']));
end
if numel(unique(Name)) ~= numel(Name)
    error(message(['MATLAB:graphfun:' errTag ':NonUniqueNames']));
end
