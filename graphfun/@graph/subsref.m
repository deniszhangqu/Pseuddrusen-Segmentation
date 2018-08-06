function [out, varargout] = subsref(G, S)

%   Copyright 2014-2016 The MathWorks, Inc.

%  Only dot subscripting is allowed on graph.
if ~strcmp(S(1).type, '.')
  error(message('MATLAB:graphfun:graph:ScalarObject'));
end
% Only allow access of public properties/methods.  Otherwise flow through...
if ~any(strcmp(S(1).subs, {'Nodes', 'Edges'}))
    if ~ismethod(G, S(1).subs)
        if ~isprop(G, S(1).subs)
            error(message('MATLAB:noSuchMethodOrField', S(1).subs, class(G)));
        else
            error(message('MATLAB:class:GetProhibited', S(1).subs, class(G)));
        end
    end
end
% Short-circuit to EdgeProperties if appropriate.
if numel(S) > 1 && strcmp(S(1).subs, 'Edges')
    if strcmp(S(2).type, '.') && ~any(strcmp(S(2).subs, {'EndNodes', 'Properties'}))
        S(1).subs = 'EdgeProperties';
    end
end
% Fall through to the builtin.
if nargout <= 1
    out = builtin('subsref', G, S);
else
    [out, varargout{1:nargout-1}] = builtin('subsref', G, S);
end
