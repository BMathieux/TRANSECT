%% ========================================================================
%  TOPOGRAPHIC TRANSIENT EVOLUTION USING TTLEM AND TRANSECT
%  Author: Bastien Mathieux (2025)
%  ========================================================================

clear; close all; clc
fprintf('\n=== TTLEM–TRANSECT SIMULATION ===\n');

%% ========================== CONFIGURATION ===============================
% --- DEM input ---
DEMfile     = 'dem_lauch_10m.tif';   % Input DEM file
resDEM      = 30;                    % DEM resampling resolution (m)
A0          = 5e5;                   % Channel area threshold (m²)
ite         = 1;                     % TRANSECT smoothing iterations

% --- TTLEM model parameters ---
p.TimeSpan   = 25e6;                 % total time (yr)
p.TimeStep   = 2e5;                  % timestep (yr)
p.D          = 1e-3;                 % diffusivity (m²/yr)
p.Kw         = 0.8151e-6;            % stream power coefficient
p.m          = 0.4168;               % area exponent
p.n          = 1.0;                  % slope exponent
p.diffScheme = 'imp_lin_sc';         % diffusion scheme
p.Sc         = 1.2;                  % critical slope
p.Sc_unit    = 'tangent';
p.DrainDir   = 'fixed';
p.DiffToRiv  = false;
p.riverInc   = 'TVD_FVM';
p.AreaThresh = A0;
p.BC_Type    = 'Dirichlet';
p.BC_dir_value = 0;
p.BC_dir_DistSites = 'tblr';
p.BC_dir_Dist_Value = 1;
p.BC_nbGhost = 2;
p.ploteach   = inf;
p.saveeach   = inf;
p = ttlemset(p);

% --- Uplift history ---
%  This uplift history is retrieved from linear inversion of river profiles of 
%  the Lauch catchment in the Vosges Mountains calibrated with TCN-derived 
%  denudation rates (Mathieux 2025, PhD thesis).
%  Users can easily replace these data with their own uplift histories.

uplift_mode = "custom";    % "custom" or "synthetic"

% (A) Example: user-defined uplift history (from inversion or data)
Udata.time_Myr = [1.21 2.45 3.77 5.07 6.43 7.77 8.78 9.72 10.53 11.29 12.06 ...
                  12.83 13.48 14.11 14.78 15.47 16.21 17.1 18.25 19.79 26.25];
Udata.rate_mMyr = [24.41 43.25 37.16 43.31 84.08 69.12 73.09 62.93 42.44 43.87 43.34...
                   38.57 54.03 48.54 42.22 28.68 28.65 25.69 24.99 22.77 19.51];

% (B) Example: simple synthetic uplift scenario
Ucfg.bg_mMyr    = 40;                % background [m/Myr]
Ucfg.pulse_mMyr = [160 180 150];     % pulse rates [m/Myr]
Ucfg.start_Myr  = [1 9 17];          % pulse start times [Myr]
Ucfg.dur_Myr    = 1;                 % duration [Myr]

fprintf('Configuration loaded.\n');

%% ===================== PREPROCESS DEM AND FLOW ==========================
fprintf('Loading and preprocessing DEM...\n');
DEM = GRIDobj(DEMfile);
DEM = resample(DEM,resDEM);
DEM = fillsinks(DEM); DEMi = DEM;

FD = FLOWobj(DEM,'preprocess','carve');
A  = flowacc(FD);
S  = STREAMobj(FD,A.*DEM.cellsize^2 > A0);
D  = DIVIDEobj(FD,S);
D  = shrink(D,FD,'distance',1e3);

baseZ = min(DEM.Z(:));
DEM.Z(:) = baseZ;
fprintf('DEM flattened to %.1f m base level.\n', baseZ);

%% ======================== BUILD UPLIFT HISTORY ==========================
tvec = 0:p.TimeStep:p.TimeSpan;  % time vector [yr]

switch uplift_mode
    case "custom"
        fprintf('Using user-provided uplift history (e.g. from inversion)...\n');
        tMyr = Udata.time_Myr;
        U_mMyr = Udata.rate_mMyr;
        t_yr = tMyr * 1e6; 
        U_yr = U_mMyr / 1e6;
        Uinterp = interp1(linspace(0,max(t_yr),numel(U_yr)),U_yr,tvec,'linear','extrap');

    case "synthetic"
        fprintf('Using synthetic uplift scenario (background + pulses)...\n');
        Uinterp = ones(size(tvec)) * (Ucfg.bg_mMyr/1e6);
        for k = 1:numel(Ucfg.pulse_mMyr)
            tk0 = Ucfg.start_Myr(k)*1e6;
            tk1 = tk0 + Ucfg.dur_Myr*1e6;
            Uinterp(tvec>=tk0 & tvec<=tk1) = Ucfg.pulse_mMyr(k)/1e6;
        end
end

fprintf('Uplift history built: %d steps (%.1f Myr total)\n',numel(Uinterp),p.TimeSpan/1e6);

%% ===================== INTERACTIVE TRANSECT SETUP =======================
fprintf('Select stream and divide...\n');
f=figure('Name','Select Stream and Divide','Color','w');
imageschs(DEMi,[],'colormap',[.7 .7 .7],'colorbar',false);
hold on
S = modify(S,'interactive','reachselect');
ixD = SelectDivide(D,DEMi,'verbose',true,'mode','single');
[xd,yd] = ind2coord(D,ixD(:,1));
close(f);

DB = drainagebasins(FD,S.IXgrid(end));
mask = logical(DB.Z);
[w,~] = maskWidth(DEM,S.x,S.y,mask,ite);
Td = TRANSECT(DEM,S.x,S.y,'w',w,'ite',ite);
Td = pairing(Td,xd,yd);

%% ======================= RUN TTLEM SIMULATION ===========================
fprintf('Running TTLEM evolution...\n');
f=figure('Color','w','Position',[100 100 1300 550]);
cmap = turbo(numel(Uinterp));

DEM_current = DEM;
DEM_series = cell(numel(Uinterp),1);
Td_series = cell(numel(Uinterp),1);

dist = [0; cumsum(sqrt(diff(Td.x).^2 + diff(Td.y).^2))]./1e3;
ix = coord2ind(DEM,Td.x,Td.y);

for i = 1:numel(Uinterp)
    Ui = GRIDobj(DEM_current);
    Ui.Z(:) = Uinterp(i);

    p.TimeSpan = p.TimeStep;
    out = ttlem(DEM_current,Ui,p,FD);
    DEM_current = fillsinks(out.H1);
    DEM_series{i} = DEM_current;

    Td = update(Td,DEM_current);
    Td_series{i} = Td;

    subplot(2,2,1)
    z = DEM_current.Z(ix);
    plot(dist,z,'Color',cmap(i,:),'LineWidth',1.2);
    xlabel('Distance (km)'), ylabel('Elevation (m)'), title('Longitudinal profile')
    grid on; box on; hold on

    subplot(2,2,2)
    plot(tvec./1e6,Uinterp*1e6,'k'); hold on
    scatter(tvec(i)./1e6,Uinterp(i)*1e6,'ro');
    xlabel('Time (Myr)'), ylabel('Uplift rate (m/Myr)')
    title('Uplift history'); grid on; box on; hold on

    subplot(2,2,3)
    if i>1
        G_now  = median([Td.stats.G1 Td.stats.G2],2,'omitnan');
        G_prev = median([Td_series{i-1}.stats.G1 Td_series{i-1}.stats.G2],2,'omitnan');
        plot(dist,G_now-G_prev,'Color',cmap(i,:),'LineWidth',1.2);
        ylabel('Δ Hillslope degree (°)'); title('Hillslope evolution');
    else
        G_now = median([Td.stats.G1 Td.stats.G2],2,'omitnan');
        plot(dist,G_now,'Color',cmap(i,:),'LineWidth',1.2);
        ylabel('Hillslope degree (°)'); title('Initial hillslope');
    end
    xlabel('Distance (km)'); grid on; box on; hold on

    subplot(2,2,4)
    G = gradient8(DEM_current,'degree');
    imageschs(DEM_current,G,'caxis',[0 prctile(G.Z(:),99)],...
            'colormap',ttscm('imola'),'colorbarlabel','Hillslope degree (°)');
    hold on; plot(Td.x,Td.y,'r','LineWidth',1);
    title(sprintf('t = %.2f Myr | U = %.2f m/Myr',tvec(i)/1e6,Uinterp(i)*1e6));
    grid on; box on; hold on

    drawnow
end

fprintf('=== TTLEM–TRANSECT simulation complete ===\n');
