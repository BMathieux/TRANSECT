function D2 = modifyD(D,DEM,FD,varargin)
%MODIFYD Modify a DIVIDEobj by order, distance, or interactive selection
%
% Syntax
%   D2 = modifyD(D,DEM,FD,'distance',500)
%   D2 = modifyD(D,DEM,FD,'distance',[500 2000])
%   D2 = modifyD(D,DEM,FD,'order','>=3')
%   D2 = modifyD(D,DEM,FD,'interactive','reachselect')
%   D2 = modifyD(D,DEM,FD,'interactive','polyselect')
%
% Description
%   modifyD filters or interactively edits a DIVIDEobj while preserving the
%   NaN-separated path structure used by TopoToolbox.
%
% Options
%   'order'        scalar or expression, e.g. 3, '>=3', '==2'
%   'distance'     scalar threshold or [min max] interval
%   'interactive'  'reachselect', 'polyselect', 'rectselect', 'ellipseselect'
%   'mode'         'single' or 'multi', used by reachselect
%   'verbose'      true or false

%% Parse inputs

p = inputParser;
p.FunctionName = 'modifyD';

addRequired(p,'D',@(x) isa(x,'DIVIDEobj'));
addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
addRequired(p,'FD',@(x) isa(x,'FLOWobj'));

addParameter(p,'verbose',true,@islogical);
addParameter(p,'mode','single',@(x) ischar(x) && any(strcmpi(x,{'single','multi'})));
addParameter(p,'interactive',[],@(x) isempty(x) || any(strcmpi(x,{'reachselect','polyselect','rectselect','ellipseselect'})));
addParameter(p,'order',[]);
addParameter(p,'distance',[],@(x) isnumeric(x) && numel(x)<=2);

parse(p,D,DEM,FD,varargin{:});

vrb = p.Results.verbose;
mod = p.Results.mode;
imt = p.Results.interactive;

if ~isempty(imt)
    ax0 = gca;
    hold(ax0,'on')
else
    ax0 = [];
end

%% Prepare divide object

if ~D.issorted
    D = sort(D);
end

if isempty(D.ordertype) || ~any(D.ordertype)
    D = divorder(D,'topo');
end

%% Interactive selection

if ~isempty(imt)

    imt = validatestring(imt,{'reachselect','polyselect','rectselect','ellipseselect'});

    switch imt

        case 'reachselect'
            IX = selectReach(D,mod,vrb,ax0);
            D2 = rebuildReach(D,IX);
            plotResult(D2,ax0)

        otherwise
            IX = selectShape(D,DEM,imt);
            kp = ismember(D.IX,IX) | isnan(D.IX);
            D2 = rebuildNet(D,FD,D.IX(kp));
            plotResult(D2,ax0)
    end

    return
end

%% Automatic filtering by distance

if ~isempty(p.Results.distance)

    if isempty(D.distance)
        D = divdist(D);
    end

    dr = p.Results.distance(:);

    if isscalar(dr)
        kp = D.distance > dr | isnan(D.IX);
    else
        kp = (D.distance >= min(dr) & D.distance <= max(dr)) | isnan(D.IX);
    end

    D2 = rebuildNet(D,FD,D.IX(kp));
    return
end

%% Automatic filtering by divide order

if ~isempty(p.Results.order)

    if isempty(D.order)
        D = divorder(D,'topo');
    end

    req = p.Results.order;

    if isnumeric(req)

        kp = D.order > req | isnan(D.IX);

    else

        tok = regexp(strtrim(char(req)),'^(>=|<=|==|~=|>|<)\s*([-+]?\d+\.?\d*)$','tokens','once');

        if isempty(tok)
            error('modifyD:WrongOrderExpression','Use expressions such as ''>=3'', ''==2'', or ''<4''.')
        end

        op = tok{1};
        th = str2double(tok{2});

        switch op
            case '>'
                kp = D.order > th | isnan(D.IX);
            case '>='
                kp = D.order >= th | isnan(D.IX);
            case '<'
                kp = D.order < th | isnan(D.IX);
            case '<='
                kp = D.order <= th | isnan(D.IX);
            case '=='
                kp = D.order == th | isnan(D.IX);
            case '~='
                kp = D.order ~= th | isnan(D.IX);
        end
    end

    D2 = rebuildNet(D,FD,D.IX(kp));
    return
end

error('modifyD:NoOption','Provide one option: ''interactive'', ''order'', or ''distance''.')
end

function D2 = rebuildNet(D,FD,IX)

IX = cleanIX(D.size,IX,2,true);

D2 = D;
D2.IX = IX;

if isempty(D2.IX)
    D2 = clearD(D2);
    return
end

D2 = divnet(D2,FD);
D2 = sort(D2);

if ~isempty(D.distance)
    D2 = divdist(D2);
end

if ~isempty(D.order)
    D2 = divorder(D2,D.ordertype);
end
end

function D2 = rebuildReach(D,IX)

IX = cleanIX(D.size,IX,1,false);

D2 = D;
D2.IX = IX;

if isempty(D2.IX)
    D2 = clearD(D2);
    return
end

[G,uix,~,~] = graphFromIX(D.size,D2.IX);
dg = degree(G);

try
    D2.ep = uix(dg == 1);
catch
end

try
    D2.jct = uix(dg > 2);
catch
end

try
    D2.jctedg = dg(dg > 2);
catch
end

try
    D2.distance = [];
catch
end

try
    D2.order = [];
catch
end

try
    D2.issorted = false;
catch
end
end

function IX = selectReach(D,mod,vrb,ax)

% Select one or several shortest paths between two snapped divide nodes

if nargin < 4 || isempty(ax) || ~isgraphics(ax)
    ax = gca;
end

fig = ancestor(ax,'figure');
axes(ax)
hold(ax,'on')

ixn = D.IX(:);
ixn = ixn(~isnan(ixn));
[xn,yn] = ind2coord(D,ixn);

[G,uix,map,cmp] = graphFromIX(D.size,D.IX);

if vrb
    disp('Starting interactive divide reach selection')
end

plot(D,'color',[0 0 0])
scatter(xn,yn,10,'.','MarkerEdgeColor',[0.5 0.5 0.5])
axis(ax,'image')
box(ax,'on')

lm = axis(ax);
xlim(ax,[lm(1)-(lm(2)-lm(1))/20,lm(2)+(lm(2)-lm(1))/20])
ylim(ax,[lm(3)-(lm(4)-lm(3))/20,lm(4)+(lm(4)-lm(3))/20])

IX = [];
done = false;

while ~done

    pth = [];
    ida = 0;
    idb = 0;
    hpath = [];

    title(ax,'Set first divide node location')
    hp1 = impoint(ax,'PositionConstraintFcn',@snapNode);
    setColor(hp1,[0 1 0])
    addNewPositionCallback(hp1,@(~)drawPath);
    setPosition(hp1,getPosition(hp1))

    title(ax,'Set last divide node location')
    hp2 = impoint(ax,'PositionConstraintFcn',@snapNode);
    setColor(hp2,[1 0 0])
    addNewPositionCallback(hp2,@(~)drawPath);
    setPosition(hp2,getPosition(hp2))

    title(ax,'Move divide nodes and press any key to extract reach')

    set(fig,'WindowKeyPressFcn',@(~,~)uiresume(fig))
    uiwait(fig)

    if isgraphics(hpath)
        set(hpath,'Color',[0.10 0.35 0.85],'LineWidth',2.2)
    elseif ~isempty(pth)
        [xp,yp] = ind2coord(D,pth);
        hpath = plot(ax,xp,yp,'Color',[0.10 0.35 0.85],'LineWidth',2.2);
    end

    delete(hp1)
    delete(hp2)

    if ~isempty(pth)
        IX = [IX; pth(:); NaN];
    end

    if strcmpi(mod,'single')
        done = true;
    else
        q = questdlg('Select another reach?','','Yes','No','No');
        done = ~strcmp(q,'Yes');
    end
end

IX = cleanIX(D.size,IX,1,false);

    function pos2 = snapNode(pos)

        [~,id] = min((xn-pos(1)).^2 + (yn-pos(2)).^2);
        pos2 = [xn(id) yn(id)];

    end

    function drawPath(~)

        p1 = getPosition(hp1);
        p2 = getPosition(hp2);

        [~,id1] = min((xn-p1(1)).^2 + (yn-p1(2)).^2);
        [~,id2] = min((xn-p2(1)).^2 + (yn-p2(2)).^2);

        ida = double(map(ixn(id1)));
        idb = double(map(ixn(id2)));

        if ida == 0 || idb == 0 || cmp(ida) ~= cmp(idb)
            pth = [];
            if isgraphics(hpath)
                set(hpath,'XData',nan,'YData',nan)
            end
            return
        end

        [P,d] = shortestpath(G,ida,idb);

        if isempty(P) || isinf(d)
            pth = [];
            if isgraphics(hpath)
                set(hpath,'XData',nan,'YData',nan)
            end
            return
        end

        pth = uix(P);
        [xp,yp] = ind2coord(D,pth);

        if isgraphics(hpath)
            set(hpath,'XData',xp,'YData',yp)
        else
            hpath = plot(ax,xp,yp,'r','LineWidth',1.5);
        end

        drawnow limitrate

    end
end

IX = cleanIX(D.size,IX,1,false);

    function pos2 = snapNode(pos)
        [~,id] = min((xn-pos(1)).^2 + (yn-pos(2)).^2);
        pos2 = [xn(id) yn(id)];
    end

    function updatePath(pos,typ)

        [~,id] = min((xn-pos(1)).^2 + (yn-pos(2)).^2);
        gid = double(map(ixn(id)));

        if typ == 1
            ida = gid;
        else
            idb = gid;
        end

        if ida == 0 || idb == 0 || cmp(ida) ~= cmp(idb)
            pth = [];
            clearPreview()
            return
        end

        [P,d] = shortestpath(G,ida,idb);

        if isempty(P) || isinf(d)
            pth = [];
            clearPreview()
            return
        end

        pth = uix(P);
        [xp,yp] = ind2coord(D,pth);

        if isgraphics(hp)
            set(hp,'XData',xp,'YData',yp)
        else
            hp = plot(ax,xp,yp,'r','LineWidth',1.5);
        end

        drawnow limitrate
    end

    function clearPreview()
        if isgraphics(hp)
            set(hp,'XData',nan,'YData',nan)
        end
    end
end

function IX = selectShape(D,DEM,meth)

fig = figure('Color','w');
imageschs(DEM,[],'colormap',[0.7 0.7 0.7],'colorbar',false)
hold on
plot(D,'color',[0 0 0])
axis image
box on

switch meth
    case 'polyselect'
        title('Create polygon and double-click to finalize')
        h = drawpolygon;
    case 'rectselect'
        title('Create rectangle and double-click to finalize')
        h = drawrectangle;
    case 'ellipseselect'
        title('Create ellipse and double-click to finalize')
        h = drawellipse;
end

wait(h)

if isprop(h,'Vertices')
    pos = h.Vertices;
else
    pos = h.Position;
end

if ~isequal(pos(1,:),pos(end,:))
    pos(end+1,:) = pos(1,:);
end

ix = D.IX(:);
ok = ~isnan(ix);

[x,y] = ind2coord(D,ix(ok));
in = inpolygon(x,y,pos(:,1),pos(:,2));

IX = ix(ok);
IX = unique(IX(in),'stable');

delete(h)
close(fig)
end

function [G,uix,map,cmp] = graphFromIX(sz,IX)

IX = IX(:);
uix = unique(IX(~isnan(IX)),'stable');

map = zeros(sz,'uint32');
map(uix) = uint32(1:numel(uix));

a = IX(1:end-1);
b = IX(2:end);

ok = ~isnan(a) & ~isnan(b) & a ~= b;
id = find(ok);

if ~isempty(id)
    [ra,ca] = ind2sub(sz,a(id));
    [rb,cb] = ind2sub(sz,b(id));
    bad = abs(ra-rb) > 1 | abs(ca-cb) > 1;
    ok(id(bad)) = false;
end

s = double(map(a(ok)));
t = double(map(b(ok)));

if isempty(s)
    G = graph([],[],[],numel(uix));
else
    E = unique(sort([s t],2),'rows','stable');
    G = graph(E(:,1),E(:,2),[],numel(uix));
end

cmp = conncomp(G);
end

function IX = cleanIX(sz,IX,nmin,breakgap)

IX = IX(:);

if isempty(IX)
    return
end

bad = isnan(IX(1:end-1)) & isnan(IX(2:end));
IX([false; bad]) = [];

while ~isempty(IX) && isnan(IX(1))
    IX(1) = [];
end

while ~isempty(IX) && isnan(IX(end))
    IX(end) = [];
end

if isempty(IX)
    return
end

if breakgap && numel(IX) > 1
    a = IX(1:end-1);
    b = IX(2:end);
    ok = ~isnan(a) & ~isnan(b);
    cut = false(numel(a),1);

    if any(ok)
        [ra,ca] = ind2sub(sz,a(ok));
        [rb,cb] = ind2sub(sz,b(ok));
        cut(ok) = abs(ra-rb) > 1 | abs(ca-cb) > 1;
    end

    if any(cut)
        IX2 = nan(numel(IX)+nnz(cut),1);
        j = 1;

        for i = 1:numel(a)
            IX2(j) = IX(i);
            j = j+1;

            if cut(i)
                IX2(j) = NaN;
                j = j+1;
            end
        end

        IX2(j) = IX(end);
        IX = IX2;
    end
end

br = [0; find(isnan(IX)); numel(IX)+1];
C = {};

for k = 1:numel(br)-1
    ix = IX(br(k)+1:br(k+1)-1);
    ix = ix(~isnan(ix));

    if numel(ix) >= nmin
        C{end+1,1} = ix(:);
    end
end

if isempty(C)
    IX = [];
    return
end

IX = nan(sum(cellfun(@numel,C)) + numel(C) - 1,1);
j = 1;

for k = 1:numel(C)
    ix = C{k};
    IX(j:j+numel(ix)-1) = ix;
    j = j+numel(ix);

    if k < numel(C)
        IX(j) = NaN;
        j = j+1;
    end
end
end

function plotResult(D,ax)

if isempty(D.IX) || ~isgraphics(ax)
    return
end

axes(ax)
hold(ax,'on')

IX = D.IX(:);
br = [0; find(isnan(IX)); numel(IX)+1];

for k = 1:numel(br)-1
    ix = IX(br(k)+1:br(k+1)-1);
    ix = ix(~isnan(ix));

    if numel(ix) > 1
        [x,y] = ind2coord(D,ix);
        plot(ax,x,y,'Color',[0.10 0.35 0.85],'LineWidth',2.2)
    end
end

drawnow
end

function D = clearD(D)

try
    D.IX = [];
catch
end

try
    D.distance = [];
catch
end

try
    D.order = [];
catch
end

try
    D.ep = [];
catch
end

try
    D.jct = [];
catch
end

try
    D.jctedg = [];
catch
end
end