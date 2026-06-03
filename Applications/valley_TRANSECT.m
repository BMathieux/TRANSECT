%% ========================================================================
%  SINGLE VALLEY ANALYSIS WITHOUT FLOODPLAIN EXTRACTION
%
%  DESCRIPTION:
%  This script performs valley geometry extraction with TRANSECT along a
%  user-selected river reach. 
%
%  WORKFLOW OVERVIEW:
%   (1) Load and preprocess the DEM.
%   (2) Build a flow network (FLOWobj), extract streams (STREAMobj),
%       and compute drainage divides (DIVIDEobj).
%   (3) Select interactively a stream reach and a divide segment.
%   (4) Run TRANSECT along the selected reach and pair transects with:
%        - the selected divide segment (Td)
%   (5) Compute and plot valley geometry metrics along the selected reach.
%
%  NOTE ON DEM RESAMPLING:
%  The script allows tracing TRANSECT paths on a coarser DEM (speed-up).
%  Only the transect geometry (node connectivity / tracing) is computed on
%  the coarser DEM. The TRANSECT object is then resampled onto the original
%  DEM so that elevations and cross-section profiles are interpolated from
%  the native-resolution topography (used for plotting and metrics).
%
%  REQUIREMENTS:
%     - TopoToolbox v3
%     - TRANSECT
%
%  Author: Bastien Mathieux (bastien.mathieux[at]gmail.com)
%  Date: February, 2026
%  ========================================================================

close all; clear; clc
fprintf('\n=== Single valley analysis ===\n');

%% ======================= USER PARAMETERS ================================
DEMfile = 'dem_lauch_10m.tif';

A0      = 1e5;                 % Critical drainage area (m^2)
ite     = 5;                   % TRANSECT iterations

doResT  = true;                % true: resample DEM for TRANSECT only
csResT  = 30;                  % target cellsize (m) for TRANSECT if doResT=true
fprintf('Parameters loaded.\n');

%% ======================= LOAD AND PREPARE DEM ===========================
fprintf('Loading and preprocessing DEM...\n');
DEM0 = GRIDobj(DEMfile);
DEM0 = fillsinks(DEM0);
DEM0.Z(DEM0.Z<280) = nan;      % remove low-elevation artefacts (Lauch-specific)

%% ======================= PREPARE TRANSECT INPUTS ========================

DEM = DEM0;

if doResT
    fprintf('Resampling DEM for TRANSECT...\n');
    DEM = resample(DEM0,csResT);
    DEM = fillsinks(DEM);
end

%% ===================== FLOW NETWORK, STREAMS, DIVIDES ===================
fprintf('Extracting flow network and divides...\n');
FD = FLOWobj(DEM,'preprocess','carve');
A  = flowacc(FD);
S  = STREAMobj(FD,A.*DEM.cellsize^2>A0);
D  = DIVIDEobj(FD,S);
D  = shrink(D,FD,3e2);

fprintf('Stream and divide networks created.\n');

%% =================== INTERACTIVE SELECTION (REACH + DIVIDE) ==============
fprintf('Opening interactive selection window...\n');
f = figure('Name','Select channel','Color','w');
imageschs(DEM,[],'colormap',[.7 .7 .7],'colorbar',false); hold on

Sall = S;
S    = modify(Sall,'interactive','reachselect');

D = modifyD(D,DEM,FD,'interactive','reachselect');
ixD = D.IX;
ixD(isnan(ixD)) = [];
[xd,yd] = ind2coord(D,ixD);

close(f)
fprintf('Stream and divide segments selected.\n');

%% ========================== RUN TRANSECT ================================
fprintf('Computing TRANSECT...\n');

DB   = drainagebasins(FD,S.IXgrid(end));
mask = logical(DB.Z);

[w,~] = maskWidth(DEM,S.x,S.y,mask,ite);

T = TRANSECT(DEM,S.x,S.y,'w',w,'ite',ite,'verbose',true);

% If TRANSECT was built on a resampled DEM, resample it back to DEM0 for plotting
if doResT
    fprintf('Resampling TRANSECT back to original DEM resolution...\n');
    T = resample(T,DEM0);
end

fprintf('Pairing and computing stats...\n');

Td = pairing(T,xd,yd);    % pair transects with divide segment

% Compute stats
Td = stats(Td);

fprintf('Computation completed.\n');

%% ============================= PLOTTING ================================
fprintf('Generating figures...\n');

% Fig1: map view of the selected reach and divide pairing
figure('Name','Transect–Divide pairing','Color','w')
imageschs(DEM0,[],'colormap',[.7 .7 .7],'colorbar',false); hold on
plot(Td)
title('TRANSECT paired with selected divide','FontSize',14,'FontWeight','bold')
axis equal


% Fig2: 3D transects + ksn along the same baseline
figure('Name','3D transects + ksn','Color','w')

du = cumsum([0;hypot(diff(Td.x),diff(Td.y))])./1e3;

for i1 = 1:numel(Td.x)
    zi  = smooth(Td.stats.z{i1}',10);
    dii = Td.stats.d{i1}'./1e3;
    plot3(dii,ones(1,numel(dii)).*du(i1),zi,'color',[0 0 0 .5]); hold on
end

k  = ksn(Sall,DEM,A,.45,1e2);
[~,j] = pdist2([Sall.x(:) Sall.y(:)],[Td.x(:) Td.y(:)],'euclidean','Smallest',1);
j  = j';
kf = k(j);

id = coord2ind(DEM0,Td.x(:),Td.y(:));
zt = DEM0.Z(id);

i1 = (1:numel(du)-1)'; 
i2 = i1+1;

cmap = ttscm('imola');

surface(zeros(2,numel(i1)),[du(i1) du(i2)]'-0.1,[zt(i1) zt(i2)]'+15,[kf(i1) kf(i2)]', ...
    'FaceColor','none','EdgeColor','interp','LineWidth',3);

colormap(cmap)
cb = colorbar;
cb.Label.String = 'ksn (m^{0.9})';

xlim([-5.25 4.5])
zlim([200 1400])
grid on; grid minor
xlabel('Distance along transect (km)');
ylabel('Distance upstream (km)')
zlabel('Elevation (m)')

set(gca,'FontWeight','bold','FontSize',22,'FontName','Arial')
set(gca,'CameraPosition',[6.1623,-136.9727,6249.3])


% Fig3: baseline long profile (ksn-colored) + transect-derived metrics
figure('Name','Long profile + valley metrics','Color','w')
tl = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

ws   = 5; % window size for metrics smoothing
cmap = flip(ttscm('roma'));

dr = du(:);
zr = zt(:);

i1 = (1:numel(dr)-1)'; 
i2 = i1+1;

z1s = movmean(Td.stats.Z1(:),ws,'omitnan');
z2s = movmean(Td.stats.Z2(:),ws,'omitnan');

m1 = isfinite(dr) & isfinite(zr) & isfinite(z1s);
m2 = isfinite(dr) & isfinite(zr) & isfinite(z2s);

ax = nexttile(tl,1);
hold(ax,'on')
cc = bwconncomp(m1);
for k=1:cc.NumObjects
    id = cc.PixelIdxList{k};
    fill(ax,[dr(id); flipud(dr(id))],[z1s(id); flipud(zr(id))],[0.5 0.5 0.5],'FaceAlpha',0.25,'EdgeColor','none');
end
cc = bwconncomp(m2);
for k=1:cc.NumObjects
    id = cc.PixelIdxList{k};
    fill(ax,[dr(id); flipud(dr(id))],[z2s(id); flipud(zr(id))],[0.5 0.5 0.5],'FaceAlpha',0.25,'EdgeColor','none');
end
surface(ax,[dr(i1) dr(i2)]',[zr(i1) zr(i2)]',zeros(2,numel(i1)),[kf(i1) kf(i2)]', ...
    'FaceColor','none','EdgeColor','interp','LineWidth',2);
plot(ax,dr,z1s,'k','LineWidth',1.2)
plot(ax,dr,z2s,'k','LineWidth',1.2)
axis(ax,'tight'); grid(ax,'on'); grid(ax,'minor')
colormap(ax,cmap)
cb = colorbar(ax); cb.Label.String = 'ksn (m^{0.9})';
xlabel(ax,'Distance upstream (km)')
ylabel(ax,'Elevation (m)')

ax = nexttile(tl,2);
yyaxis(ax,'left')
surface(ax,[dr(i1) dr(i2)]',[zr(i1) zr(i2)]',zeros(2,numel(i1)),[kf(i1) kf(i2)]', ...
    'FaceColor','none','EdgeColor','interp','LineWidth',2);
axis(ax,'tight'); grid(ax,'on'); grid(ax,'minor')
colormap(ax,cmap)
ylabel(ax,'Elevation (m)')
yyaxis(ax,'right')
y = atan(Td.stats.H(:)./Td.stats.W(:)).*180/pi;
plot(ax,dr,movmean(y,ws,'omitnan'),'k','LineWidth',1.5)
ylabel(ax,'Hillslope degree (°)')
xlabel(ax,'Distance upstream (km)')
set(ax,'YColor',[0 0 0])

ax = nexttile(tl,3);
yyaxis(ax,'left')
surface(ax,[dr(i1) dr(i2)]',[zr(i1) zr(i2)]',zeros(2,numel(i1)),[kf(i1) kf(i2)]', ...
    'FaceColor','none','EdgeColor','interp','LineWidth',2);
axis(ax,'tight'); grid(ax,'on'); grid(ax,'minor')
colormap(ax,cmap)
ylabel(ax,'Elevation (m)')
yyaxis(ax,'right')
y = Td.stats.W(:)./Td.stats.H(:);
plot(ax,dr,movmean(movmean(y,ws,'omitnan'),ws,'omitnan'),'k','LineWidth',1.5)
ylabel(ax,'Valley width/height ratio')
xlabel(ax,'Distance upstream (km)')
set(ax,'YColor',[0 0 0])

ax = nexttile(tl,4);
yyaxis(ax,'left')
surface(ax,[dr(i1) dr(i2)]',[zr(i1) zr(i2)]',zeros(2,numel(i1)),[kf(i1) kf(i2)]', ...
    'FaceColor','none','EdgeColor','interp','LineWidth',2);
axis(ax,'tight'); grid(ax,'on'); grid(ax,'minor')
colormap(ax,cmap)
ylabel(ax,'Elevation (m)')
yyaxis(ax,'right')
y = Td.stats.W(:);
plot(ax,dr,movmean(movmean(y,ws,'omitnan'),ws,'omitnan'),'k','LineWidth',1.5)
ylabel(ax,'Valley width (m)')
xlabel(ax,'Distance upstream (km)')
set(ax,'YColor',[0 0 0])

set(findall(gcf,'Type','axes'),'FontWeight','bold','FontSize',18,'FontName','Arial')