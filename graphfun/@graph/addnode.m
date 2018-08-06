function H = addnode(G, N)
%ADDNODE Add nodes to a graph
%   H = ADDNODE(G,nodeIDs) returns graph H that is equivalent to G but
%   with nodes specified by the string (or cell array of strings) nodeIDs
%   added to it. nodeIDs must specify nodes not already present in G.
%
%   H = ADDNODE(G,N) returns graph H that is equivalent to G but with N
%   nodes added.  N must be a nonnegative numeric scalar.
%
%   H = ADDNODE(G,NodeProps) returns graph H that is equivalent to G but
%   with as many nodes as there are rows in the table NodeProps added to
%   it.  NodeProps must be able to be concatenated with G.Nodes.
%
%   Example:
%       % Construct a graph with named nodes, then add two named nodes.
%       G = graph({'A' 'B' 'C'},{'D' 'C' 'D'})
%       G.Nodes
%       G = addnode(G,{'E' 'F'})
%       G.Nodes
%
%   Example:
%       % Construct a graph with four nodes, then add two nodes.
%       G = graph([1 2 3],[2 3 4])
%       G = addnode(G,2)
%
%   See also GRAPH, NUMNODES, RMNODE, ADDEDGE

%   Copyright 2014-2015 The MathWorks, Inc.

H = G;
if graph.isvalidstring(N)
    if ~iscellstr(N), N = {N}; end
    newnodes = numel(N);
    if ~hasNodeNames(H)
        H.NodeProperties.Name = makeNodeNames(numnodes(H), 0);
    end
    N = graph.validateName(N(:));
    if any(ismember(N, H.NodeProperties.Name))
        error(message('MATLAB:graphfun:graph:NonUniqueNames'));
    end
    H.NodeProperties.Name(end+1:end+newnodes,1) = N;
elseif isnumeric(N) && isscalar(N)
    if ~isreal(N) || fix(N) ~= N || N < 0
        error(message('MATLAB:graphfun:addnode:InvalidNrNodes'));
    end
    if hasNodeNames(H)
        H.NodeProperties.Name(end+1:end+N,1) = makeNodeNames(N, numnodes(H));
    else
        if isempty(H.NodeProperties)
            H.NodeProperties = table.empty(numnodes(G)+N,0);
        else
            % Copy in the first row of Node Props to get expansion behavior
            % then delete it!
            H.NodeProperties(end+N+1,:) = H.NodeProperties(1,:);
            H.NodeProperties(end,:) = [];
        end
    end
    newnodes = N;
elseif istable(N)
    N = graph.validateNodeProperties(N);
    if hasNodeNames(H) && any(strcmp(N.Properties.VariableNames, 'Name'))
        if any(ismember(N.Name, H.NodeProperties.Name))
            error(message('MATLAB:graphfun:graph:NonUniqueNames'));
        end
    end
    H.NodeProperties = [H.NodeProperties; N];
    newnodes = size(N,1);
else
    error(message('MATLAB:graphfun:addnode:SecondInput'));
end
if newnodes > 0
    A = adjacency(H.Underlying);
    A(end+newnodes, end+newnodes) = 0;
    H.Underlying = matlab.internal.graph.MLGraph(A);
end
if nargout < 1
    warning(message('MATLAB:graphfun:addnode:NoOutput'));
end

function C = makeNodeNames(numNodes, firstVal)
if numNodes > 0
    C = cellstr([repmat('Node', numNodes, 1) ...
            num2str(firstVal+(1:numNodes)', '%-d')]);    
else
    C = cell(0,1);
end
