function [Hmin,FP] = SelectHmin(DEM,HW,hr,hs)
%% ========================================================================
%  SELECTHMIN — INTERACTIVE FLOODPLAIN THRESHOLD SELECTION (GRAPHFLOOD)
%  Author: Bastien Mathieux (2026)
%
%  DESCRIPTION:
%  Interactive selection of Hmin to convert graphflood water height (HW)
%  into a binary floodplain mask (FP).
%
%  SPEED STRATEGY:
%  Panel (B) is precomputed as RGB images (uint8) using imageschs so that
%  moving the slider only swaps:
%     - the displayed RGB image
%     - the red outline coordinates
%     - the title text
%
%  SLIDER CONTROL:
%    hr = [Hmin_min Hmin_max]   (default [0 20])
%    hs = slider step           (default 1)        
%
%  DISPLAY (2 PANELS):
%    (A) Fixed reference with Hmin = 0 m (colorbar shown)
%    (B) Current Hmin (slider) + FP outline in red (bwboundaries, 8-conn)
%
%  ========================================================================

if nargin<2, error('SelectHmin(DEM,HW) needs at least two inputs.'); end
if nargin<3 || isempty(hr), hr=[0 20]; end
if nargin<4 || isempty(hs), hs=1; end

if isa(HW,'GRIDobj')
    HWZ = HW.Z;
elseif isstruct(HW) && isfield(HW,'Z')
    HWZ = HW.Z;
else
    HWZ = HW;
end

if ~isnumeric(HWZ), error('HW must be numeric, GRIDobj, or struct with field .Z.'); end
if ~isequal(size(HWZ),size(DEM.Z)), error('HW grid must match DEM.Z size.'); end

hr = sort(hr(:))'; hr = hr(1:2);
hs = abs(hs); if hs==0, hs=1; end

hv = (hr(1):hs:hr(2))';
if isempty(hv) || hv(end)~=hr(2), hv = unique([hv; hr(2)]); end
nh = numel(hv);

Hmin = hv(1);
FP   = [];

mA = (HWZ>=0) & ~isnan(HWZ);
ZA = HWZ; 
ZA(~mA) = nan;

cmap = flip(ttscm('batlowW'));
nanC = [80 85 80]./255;

[~,R] = GRIDobj2im(DEM);

%% ======================= PRECOMPUTE RGB + OUTLINES ======================
fprintf('SelectHmin: precomputing %d images...\n',nh);

rgb = cell(nh,1);
are = zeros(nh,1);
bdx = cell(nh,1);
bdy = cell(nh,1);

ft = figure('Visible','off','Color','w');
at = axes(ft);

for i=1:nh
    m = (HWZ>=hv(i)) & ~isnan(HWZ);
    are(i) = sum(m(:))*DEM.cellsize^2/1e6;

    Zc = HWZ; 
    Zc(~m) = nan;

    G = DEM; 
    G.Z = Zc;

    cla(at)
    RGB = imageschs(G,[],'colormap',cmap,'nancolor',nanC,'colorbar',false,'caxis',[0 prctile(ZA(:),99)]);
    rgb{i} = uint8(double(RGB));

    B = bwboundaries(m,8,'noholes');
    xx=[]; yy=[];
    for j=1:numel(B)
        id = sub2ind(size(m),B{j}(:,1),B{j}(:,2));
        [x,y] = ind2coord(DEM,id);
        xx=[xx; x; nan];
        yy=[yy; y; nan];
    end
    bdx{i}=xx; 
    bdy{i}=yy;
end

close(ft)

%% ======================= PANEL (A): FIXED Hmin = 0 =======================

ft = figure('Visible','off','Color','w');
at = axes(ft);

G = DEM; 
G.Z = ZA;

RGB0 = imageschs(G,[],'colormap',cmap,'nancolor',nanC,'colorbar',false,'caxis',[0 prctile(ZA(:),99)]);
rgb0 = uint8(double(RGB0));

close(ft)

%% ======================= FIGURE LAYOUT ==================================
fig = figure('Name','Select Hmin','Color','w','Units','normalized','Position',[.08 .08 .84 .82]);

tl  = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
axA = nexttile(tl,1);
axB = nexttile(tl,2);

if hr(2)==hr(1)
    stp = 1;
else
    stp = min(1,max(hs/(hr(2)-hr(1)),1/200));
end

sld = uicontrol(fig,'Style','slider','Min',hr(1),'Max',hr(2),'Value',Hmin,...
    'Units','normalized','Position',[.25 .015 .5 .03],'SliderStep',[stp stp],'Interruptible','off');

txt = uicontrol(fig,'Style','text','Units','normalized','Position',[.01 .012 .1 .035],...
    'String',sprintf('Hmin = %.3g m',Hmin),'FontWeight','bold','BackgroundColor','w','HorizontalAlignment','left');

btn = uicontrol(fig,'Style','pushbutton','Units','normalized','Position',[.78 .012 .2 .035],...
    'String','OK','Interruptible','off');

uistack(sld,'top'); uistack(txt,'top'); uistack(btn,'top');

%% ======================= PANEL (A): SHOW FIXED RGB + COLORBAR ============
axes(axA)
imshow(flipud(rgb0),R); set(axA,'YDir','normal')
colormap(axA,cmap); caxis(axA,[0 prctile(ZA(:),99)])
cbA = colorbar(axA); ylabel(cbA,'Water height (m)')
title(axA,'A: Reference (Hmin = 0 m)','FontWeight','bold')

%% ======================= PANEL (B): INIT OBJECTS =========================
axes(axB)
imB = imshow(flipud(rgb{1}),R); hold on
set(axB,'YDir','normal')
colormap(axB,cmap); caxis(axB,[0 prctile(ZA(:),99)])
lnB = plot(axB,nan,nan,'r','LineWidth',2);
title(axB,sprintf('B: Hmin = %.3g m | FP area = %.2f km^2',hv(1),are(1)),'FontWeight','bold')

%% ======================= STATE ==========================================
st.hv  = hv;
st.hs  = hs;
st.hr  = hr;

st.rgb = rgb;
st.are = are;
st.bdx = bdx;
st.bdy = bdy;

st.axB = axB;
st.imB = imB;
st.lnB = lnB;

st.sld = sld;
st.txt = txt;
st.btn = btn;

st.hm  = Hmin;

guidata(fig,st);

set(sld,'Callback',@(s,e)upd(fig));
set(btn,'Callback',@(s,e)done(fig));
set(fig,'CloseRequestFcn',@(s,e)done(fig));

%% ======================= INITIAL DRAW ===================================
upd(fig);
uiwait(fig);

if isvalid(fig)
    st   = guidata(fig);
    Hmin = st.hm;
    FP   = (HWZ>=Hmin) & ~isnan(HWZ);
    delete(fig);
end

%% ======================= NESTED CALLBACKS ===============================
    function upd(fig)
        st = guidata(fig);

        hm = get(st.sld,'Value');
        hm = st.hr(1) + round((hm-st.hr(1))/st.hs)*st.hs;
        hm = max(st.hr(1),min(st.hr(2),hm));
        set(st.sld,'Value',hm);
        set(st.txt,'String',sprintf('Hmin = %.3g m',hm));

        k = round((hm-st.hv(1))/st.hs)+1;
        k = max(1,min(numel(st.hv),k));
        st.hm = st.hv(k);

        set(st.imB,'CData',flipud(st.rgb{k}));
        set(st.lnB,'XData',st.bdx{k},'YData',st.bdy{k});
        title(st.axB,sprintf('B: Hmin = %.3g m | FP area = %.2f km^2',st.hm,st.are(k)),'FontWeight','bold')

        uistack(st.sld,'top'); uistack(st.txt,'top'); uistack(st.btn,'top');
        guidata(fig,st);
    end

    function done(fig)
        if isvalid(fig), uiresume(fig); end
    end
end
