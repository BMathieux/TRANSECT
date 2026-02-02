%% ========================================================================
%  SINGLE VALLEY ANALYSIS WITH GRAPHFLOOD-BASED FLOODPLAIN EXTRACTION
%  Author: Bastien Mathieux (2025)
%
%  DESCRIPTION:
%  This script performs a morphometric analysis by running a TRANSECT 
%  along a user-selected stream (STREAMobj) and pairing it with:
%     (1) a user-selected drainage divide (DIVIDEobj)
%     (2) a floodplain mask automatically extracted using graphflood
%
%  The analysis quantifies valley and floodplain geometry for comparison 
%  along the transect and visualizes results in two figures:
%     • Map view: Transect paired with selected divide
%     • 3D view: Transect elevation and floodplain geometry
%
%  REQUIREMENTS:
%     - TopoToolbox v3
%     - TRANSECT toolbox
%
%  ========================================================================

close all; clear; clc
fprintf('\n=== Single valley analysis ===\n');

%% ======================= USER PARAMETERS ================================
DEMfile = 'dem_lauch_10m.tif';   % Input DEM file
A0      = 1e5;                   % Minimum drainage area threshold (m²)
ite     = 5;                     % Iterations in TRANSECT
P       = (10-4)/3600*100;       % Precipitation for graphflood
iter    = 1e2;                   % Graphflood iterations
Hmin    = 5;                     % Minimum water depth for floodplain mask (m)
fprintf('Parameters loaded.\n');

%% ======================= LOAD AND PREPARE DEM ===========================
fprintf('Loading and preprocessing DEM...\n');
DEM = GRIDobj(DEMfile);
DEM = fillsinks(DEM);
DEM.Z(DEM.Z<280) = nan; % remove low-elevation artefacts (specific to Lauch DEM)

%% ====================== FLOODPLAIN EXTRACTION ===========================
fprintf('Extracting floodplain mask using graphflood...\n');
HW = graphflood(DEM,P,'N_iterations',iter);
FP = HW.Z; 
FP(FP<Hmin | isnan(FP)) = 0; 
FP = logical(FP);
fprintf('Floodplain mask extracted (Hmin = %.1f m)\n',Hmin);

%% ===================== EXTRACT STREAM AND DIVIDES ======================
fprintf('Extracting flow network and divides...\n');
FD = FLOWobj(DEM,'preprocess','carve');
A  = flowacc(FD);
S  = STREAMobj(FD,A.*DEM.cellsize^2>A0);
D  = DIVIDEobj(FD,S);
D  = shrink(D,FD,'distance',1e3); % optional: smooth divide geometry
fprintf('Stream and divide networks created.\n');

%% =================== SELECT STREAM AND DIVIDE SEGMENTS ==================
fprintf('Opening interactive selection window...\n');
f = figure('Name','Select Stream and Divide','Color','w');
imageschs(DEM,[],'colormap',[.7 .7 .7],'colorbar',false); hold on

% Select stream reach interactively
S = modify(S,'interactive','reachselect');

% Select divide segment interactively
ixD = SelectDivide(D,DEM,'verbose',true,'mode','single');
ixD = ixD(:,1);
ixD(isnan(ixD)) = [];
[xd,yd] = ind2coord(D,ixD);
close(f)
fprintf('Stream and divide segments selected.\n');

%% ========================== RUN TRANSECT ================================
fprintf('Computing TRANSECT and pairing...\n');

% Delineate drainage basin mask for the selected stream
DB = drainagebasins(FD,S.IXgrid(end));
mask = logical(DB.Z);

% Compute valley width from basin outline
[w,~] = maskWidth(DEM,S.x,S.y,mask,ite);

% Create TRANSECT object
T = TRANSECT(DEM,S.x,S.y,'w',w,'ite',ite,'verbose',true);

% Pair with selected divide
Td = pairing(T,xd,yd);

% Pair with floodplain mask
Tfp = pairing(T,FP);

fprintf('TRANSECT computation and pairing complete.\n');

%% ============================= PLOTTING ================================
fprintf('Generating figures...\n');

% --- Figure 1: TRANSECT paired with divide (map view)
figure('Name','Transect–Divide Pairing','Color','w')
imageschs(DEM,[],'colormap',[.7 .7 .7],'colorbar',false); hold on
plot(Td)
title('TRANSECT paired with selected divide','FontSize',14,'FontWeight','bold')
axis equal

% --- Figure 2: 3D transect with floodplain overlay
figure('Name','Transect and Floodplain 3D','Color','w')
di = cumsum([0;sqrt(diff(Td.x).^2+diff(Td.y).^2)])./1e3;
di = flip(di);

for i1 = 1:numel(Td.x)
    zi = Td.stats.z{i1}'; zi = smooth(zi,10);
    dii = Td.stats.d{i1}'./1e3;
    plot3(dii,ones(1,numel(Td.stats.d{i1})).*di(i1),zi,'color',[0 0 0 .5]); hold on
    plot3(Tfp.stats.d{i1}'./1e3,ones(1,numel(Tfp.stats.d{i1})).*di(i1),...
          Tfp.stats.z{i1}'-5,'color',[1 .3 .3],'LineWidth',3);
end

grid on; grid minor; box on
xlabel('Distance along transect (km)');
ylabel('Distance along trunk (km)');
zlabel('Elevation (m)');
set(gca,'FontWeight','bold','FontSize',22,'FontName','Arial')

fprintf('Figures generated successfully.\n');
fprintf('=== Analysis complete ===\n\n');

