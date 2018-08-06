function xy = packLayouts(xy,compcell,comp)
%packLayouts Pack different layouts into a rectangle
%
%   FOR INTERNAL USE ONLY -- This feature is intentionally undocumented.
%   Its behavior may change, or it may be removed in a future release.
%
%   XY = packLayouts(XY,COMPCELL,COMP) Packs a collection of different
%   layouts into a rectangle, such that the layouts don't overlap. XY has
%   size NUMNODES-by-2 and contains the node coordinates. COMPCELL is a
%   cell array where cell k contains the nodes forming layout k. COMP is a
%   vector of length NUMNODES specifying the layout number for each node.

%   Copyright 2015 The MathWorks, Inc.

% Stack components left-justified on a vertical strip of fixed width.
[xy,complevel,wmax,levelh] = packIntoStripFFDH(xy,compcell);
% We now have a TALL strip of width wmax. Make it as square as possible.
xy = foldStripIntoSquare(xy,complevel,wmax,levelh,comp);

%--------------------------------------------------------------------------
function [xy,complevel,wmax,levelh] = packIntoStripFFDH(xy,compcell)
% Stack components left-justified on a vertical strip of fixed width, i.e.,
% First-Fit Decreasing-Height (FFDH) strip packing. Reference:
%   E. G. Coffman, Jr., M. R. Garey, D. S. Johnson, and R. E. Tarjan,
%   "Performance Bounds for Level-Oriented Two-Dimensional Packing
%   Algorithms", SIAM J. Comput., 9(4), pp. 808-826, 1980.
ncomp = numel(compcell);
% Bounding boxes for each connected component.
wh = zeros(ncomp,2); % [width, height]
for i = 1:ncomp
    xyi = xy(compcell{i},:);
    minxyi = min(xyi,[],1);
    wh(i,:) = max(xyi,[],1) - minxyi;
    % Prefer tall boxes.
    if wh(i,1) > wh(i,2)
        wh(i,:) = wh(i,[2 1]);
        minxyi = minxyi(:,[2 1]);
        xyi = xyi(:,[2 1]);
    end
    % Bounding boxes with lower left corners at (0,0).
    xy(compcell{i},:) = xyi - ones(size(xyi,1),1)*minxyi;
end
% Pad boxes so boundary nodes don't overlap with neighboring boxes.
padxy = 0.15*max(wh,[],1);
padxy( padxy == 0 ) = 0.5;
xy = xy + ones(size(xy,1),1)*padxy;
wh = wh + ones(size(wh,1),1)*2*padxy;
% Sort boxes according to decreasing height.
[~,hind] = sort(wh(:,2),'descend');
wmax = max(wh(:,1));
% Place highest bounding box left-justified on the first level.
complevel = zeros(ncomp,1);
complevel(hind(1)) = 1;    % levels where the components sit
levelw(1) = wh(hind(1),1); % occupied width of each level
levelh(1) = wh(hind(1),2); % height of each level
% Greedy placement of remaining boxes: First-Fit Decreasing-Height packing.
for i = 2:ncomp
    compind = hind(i);
    [lvl,dxy,levelw,levelh] = findLevel(wh(compind,:),wmax,levelw,levelh);
    tmpindi = compcell{compind};
    xy(tmpindi,:) = xy(tmpindi,:) + ones(length(tmpindi),1)*dxy;
    complevel(compind) = lvl;
end
%--------------------------------------------------------------------------
function [level,dxy,levelw,levelh] = findLevel(wh,wmax,levelw,levelh)
% For a partially filled strip of width wmax, find a level to place a box
% that has width wh(1), height wh(2) and lower left corner at (0,0).
level = find((wmax - levelw) >= wh(1),1);
if isempty(level)
    dxy = [0, sum(levelh)]; % x and y offset
    levelw(end+1) = wh(1);
    % Use the fact that boxes are ordered in descending height.
    levelh(end+1) = wh(2);
    level = length(levelw);
else
    dxy = [levelw(level), sum(levelh(1:(level-1)))]; % x and y offset
    levelw(level) = levelw(level) + wh(1);
end
%--------------------------------------------------------------------------
function xy = foldStripIntoSquare(xy,complevel,stripw,levelheight,comp)
% Re-arrange a tall strip of rectangles of same width into an almost square
% shape. The rectangles are sorted according to height: highest at bottom.
levely = [0 cumsum(levelheight(1:end-1))];
striph = sum(levelheight);
c = ceil(sqrt(striph/stripw)); % stripw is never 0, because we padded
if c >= 2 % if new width is a multiple of old width
    newh = striph/c; % neww = c*stripw;
    % Place top of strip at y = 0 and shift it towards right by stripw.
    % Repeat until nothing is above the new height.
    levelind = find(levely > newh); %levels to shift down,towards right
    c = c-1; % safeguard
    while ~isempty(levelind) && c > 0
        dxy = [stripw, -levely(levelind(1))]; % shift
        ind = find(complevel >= levelind(1)); % components above newh
        ind = ismember(comp,ind);             % nodes above newh
        xy(ind,:) = xy(ind,:) + ones(nnz(ind),1)*dxy;
        levely(levelind) = levely(levelind) + dxy(2);
        levelind = find(levely > newh); %levels to shift down,towards right
        c = c-1;
    end
end
