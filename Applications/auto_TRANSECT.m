%% ========================================================================
%  AUTOMATIC TRANSECT EXTRACTION WITH GRAPHFLOOD-BASED FLOODPLAIN DETECTION
%  Author: Bastien Mathieux (2025)
%
%  DESCRIPTION:
%  This script automatically extracts valley and floodplain transects 
%  across an entire mountain range using:
%      (1) TopoToolbox flow network analysis
%      (2) Graphflood-based floodplain detection
%      (3) The TRANSECT toolbox for valley geometry extraction
%
%  The script produces two figures:
%      • Map of floodplain width (m)
%      • Map of valley height/width ratio
%
%  REQUIREMENTS:
%      - TopoToolbox v3
%      - TRANSECT toolbox
%
%  ========================================================================

close all; clear; clc
fprintf('\n=== AUTOMATIC TRANSECT EXTRACTION ===\n');

%% ======================= USER PARAMETERS =================================
DEMfile = 'dem_lauch_10m.tif';   % Input DEM file
A0      = 1e5;                   % Minimum drainage area threshold (m²)
P       = (10-4)/3600*100;       % Precipitation proxy for graphflood
iter    = 1e2;                   % Number of graphflood iterations
Hmin    = 5;                     % Minimum water depth for floodplain mask (m)
ite     = 1;                     % Number of iterations in TRANSECT
npool   = 4;                     % Parallel pool size
doPar   = true;                  % Enable/disable parallel processing
fprintf('Parameters loaded.\n');

%% ========================= DEM PREPARATION ===============================
fprintf('Loading and preprocessing DEM...\n');
DEM = GRIDobj(DEMfile);
DEM = fillsinks(DEM);
DEM.Z(DEM.Z<280) = nan; % remove low-elevation artefacts (specific to Lauch DEM)

%% ====================== FLOODPLAIN EXTRACTION ============================
fprintf('Extracting floodplain mask using graphflood...\n');
HW = graphflood(DEM,P,'N_iterations',iter);
FP = HW.Z; 
FP(FP<Hmin | isnan(FP)) = 0; 
FP = logical(FP);
fprintf('Floodplain mask extracted (Hmin = %.1f m)\n',Hmin);

%% ===================== STREAM NETWORK PROCESSING ========================
fprintf('Building stream network and branches...\n');
DEMi = DEM; 
DEM = resample(DEM,30); % downsample for efficiency
FD = FLOWobj(DEM,'preprocess','carve');
A  = flowacc(FD);
S  = STREAMobj(FD,A.*DEM.cellsize^2>A0);
S  = klargestconncomps(S);

MS = sidebranching(S,'mapstruct');

% Remove first-order tributaries and trim minor branches
for i = numel(MS):-1:1
    if MS(i).streamorder==1
        MS(i) = [];
    elseif MS(i).streamorder==2
        ix = coord2ind(DEM,MS(i).X(1:end-1),MS(i).Y(1:end-1));
        S2 = trunk(modify(S,'upstreamto',(ix(end-1))));
        MS(i).X = S2.x(S2.ixc);
        MS(i).Y = S2.y(S2.ixc);
    end
end
fprintf('Stream branches prepared: %d segments\n',numel(MS));

%% ===================== PARALLEL OR SEQUENTIAL MODE =======================
if doPar
    if isempty(gcp('nocreate')), parpool(npool); end
    mode = 'Parallel';
else
    delete(gcp('nocreate'));
    mode = 'Sequential';
end
fprintf('Running in %s mode.\n',mode);

%% ========================= TRANSECT EXTRACTION ==========================
fprintf('Extracting TRANSECT objects (valley and floodplain)...\n');
T_valley = cell(numel(MS),1);
T_floodplain = cell(numel(MS),1);
PB = ProgressBar(numel(MS)*2,'taskname','Extracting TRANSECT objects...','ui','cli');

parfor (i = 1:numel(MS), doPar*npool + ~doPar)
    xa = MS(i).X(1:end-1); 
    ya = MS(i).Y(1:end-1);
    if numel(xa) < 5
        count(PB); count(PB); 
        continue
    end

    % Delineate drainage basin for current reach
    ix = coord2ind(DEM,xa(end),ya(end));
    DB = drainagebasins(FD,ix);
    mask = logical(DB.Z);

    % Extract basin boundary (TopoToolbox v3)
    b = bwboundaries(mask,8);
    [xb,yb] = sub2coord(DEM,b{1}(:,1),b{1}(:,2));

    % Compute valley width and create TRANSECT
    [w,~] = maskWidth(DEM,xa,ya,mask,ite);
    T_valley{i} = TRANSECT(DEM,xa,ya,'w',w,'ite',ite,'verbose',false);
    count(PB)

    % Pair valley transect with floodplain and outline
    T_floodplain{i} = pairing(T_valley{i}.resample(DEMi),FP);
    T_valley{i}     = pairing(T_valley{i},xb,yb);
    count(PB)
end
fprintf('\nTRANSECT extraction complete.\n');

%% ========================== VISUALIZATION ===============================
fprintf('Generating visualization...\n');

% --- Floodplain width map ---
figure('Name','Floodplain Width','Color','w')
imageschs(DEM,[],'colormap',[.4 .4 .4],'colorbar',false); hold on
for i = 1:numel(T_floodplain)
    if ~isempty(T_floodplain{i})
        x = T_floodplain{i}.x; 
        y = T_floodplain{i}.y; 
        w = log10(T_floodplain{i}.stats.W);
        surf([x';x'],[y';y'],zeros(2,numel(x)),[w';w'],...
            'edgecolor','interp','facecolor','none','linew',2)
    end
end
cb = colorbar; cb.Label.String = 'log(Floodplain width)';
colormap(flip(ttscm('batlowW'))); axis equal
title('Floodplain Width','FontWeight','bold')

% --- Valley aspect ratio map ---
figure('Name','Valley Aspect Ratio','Color','w')
imageschs(DEM,[],'colormap',[.4 .4 .4],'colorbar',false); hold on
for i = 1:numel(T_valley)
    if ~isempty(T_valley{i})
        x = T_valley{i}.x; 
        y = T_valley{i}.y; 
        r = T_valley{i}.stats.H ./ T_valley{i}.stats.W;
        surf([x';x'],[y';y'],zeros(2,numel(x)),[r';r'],...
            'edgecolor','interp','facecolor','none','linew',2)
    end
end
cb = colorbar; cb.Label.String = 'Valley height / width ratio';
colormap(ttscm('imola')); axis equal
title('Valley Height/Width Ratio','FontWeight','bold')

fprintf('Visualization complete.\n');
fprintf('=== Analysis finished ===\n\n');
