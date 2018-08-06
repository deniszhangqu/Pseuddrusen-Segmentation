function shortestdistalgo(Weight_matrix,startterminal_matrix,endterminal_matrix,startid,finishid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%      THIS WORK IS SUBMITTED BY:
%%
%%      ALOK PATEL
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This function finds shortest distance between nodes using Dijkstra
% Algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function contains five input algorithms.
% startid represents the starting node.

% finishid represents the end node.

%Weight_matrix denotes the weight matrix for different paths. 

%startterminal_matrix and endterminal_matrix denotes the matrixes that are

%starting and end terminals arranged in line for example:

%Weight_matrix=[2 2 6 1 2 4 7 3 2 3 2];

%startterminal_matrix=[1 2 1 5 5 7 2 3 6 3 4];

%endterminal_matrix=[2 5 7 7 6 8 3 6 8 4 8]; 

%this will give following output
%  (1,2)        2 (i.e path is from 1 to 2 with weight 2.)
%  (2,3)        7 
%  (3,4)        3
%  (2,5)        2
%  (3,6)        3
%  (5,6)        2
%  (1,7)        6
%  (5,7)        1
%  (4,8)        2
%  (6,8)        2
%  (7,8)        4


%% Implementation Example:

%            a=[2 2 6 1 2 4 7 3 2 3 2];
%            b=[1 2 1 5 5 7 2 3 6 3 4];
%            c=[2 5 7 7 6 8 3 6 8 4 8];
%            shortestdistalgo(a,b,c,1,8)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This file contains algorithm for both directed and non-directed graphs. 
% comment out Directed graph part of the function if you dont need it.
% If you need only non-directed graph comment out the lines below it.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
W = Weight_matrix;
k=max(startterminal_matrix);
l=max(endterminal_matrix);
m=max([k l])
%% Shortest Distance algorithm for Directed Graph
DG = sparse(startterminal_matrix,endterminal_matrix,W,m,m)
h = view(biograph(DG,[],'ShowWeights','on'));
[dist,path,pred] = graphshortestpath(DG,startid,finishid);
set(h.Nodes(path),'Color',[1 0.4 0.4]);
edges = getedgesbynodeid(h,get(h.Nodes(path),'ID'));
set(edges,'LineColor',[1 0 0]);
set(edges,'LineWidth',1.5);
%% Shortest Distance algorithm for Non-Directed Graph
UG = tril(DG + DG');
h = view(biograph(UG,[],'ShowArrows','off','ShowWeights','on'));
[dist,path,pred] = graphshortestpath(UG,startid,finishid,'directed',false);
set(h.Nodes(path),'Color',[1 0.4 0.4]);
 fowEdges = getedgesbynodeid(h,get(h.Nodes(path),'ID'));
 revEdges = getedgesbynodeid(h,get(h.Nodes(fliplr(path)),'ID'));
 edges = [fowEdges;revEdges];
 set(edges,'LineColor',[1 0 0]);
set(edges,'LineWidth',1.5);
end