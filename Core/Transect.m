classdef Transect
    %Transect Transect object from baseline coordinates on DEM
    %
    % Syntax
    %
    %     obj = Transect(DEM,x,y)
    %     obj = Transect(DEM,x,y,'pn1',pv1,...)
    %
    % Description
    %
    %     Transect constructs a transect object along the baseline defined by 
    %     coordinates x and y on the digital elevation model DEM (GRIDobj).
    %     The transect is built perpendicular to the baseline with a specified width 
    %     w and number of iterations ite, allowing for detailed topographic profiling.
    %
    %     The object pads the DEM to avoid boundary issues and constructs digraphs 
    %     for path calculations. Data is stored separately for each side of the baseline.
    %
    %     An additional option 'method' allows switching between 'geometric' (default, perpendicular buffer-based) 
    %     and 'flow' (flow routing along topographic divides).
    %
    % Input arguments
    %
    %     DEM        GRIDobj containing the digital elevation model
    %     x          vector of x-coordinates of the baseline
    %     y          vector of y-coordinates of the baseline
    %
    % Parameter name/value pairs
    %
    %     'w'        width of the transect (cells). Default = 5. Must be >3.
    %     'ite'      number of iterations. Default = 1. Must be >=1.
    %     'verbose'  logical scalar. Display progress and warnings. Default = true.
    %     'method'   string: 'geometric' (default) or 'flow'.
    %
    % Properties
    %
    %     type        string: Method used ('geometric' or 'flow')
    %     DEM         GRIDobj: Digital elevation model
    %     x           double: x-coordinates of the baseline
    %     y           double: y-coordinates of the baseline
    %     conn        cell{1x2}: Connection data for each side; each cell contains a struct with:
    %                   ix: cell array of path indices (baseline nodes -> paths -> vector)
    %                   x: cell array of path x coordinates
    %                   y: cell array of path y coordinates
    %                   z: cell array of path elevations
    %                   d: cell array of path distances
    %     int         cell{1x2}: Interpolated data for each side; structured similarly to conn
    %     stats       struct: Statistics with fields x,y,z,d {num_nodes,1}, Z1,G1,Z2,G2,W,H (num_nodes,1)
    %
    % Example
    %
    %     % Assume DEM is loaded
    %     x = [1 2 3]; y = [4 5 6];
    %     obj = Transect(DEM,x,y,'w',5,'ite',2);
    %     % Access data: obj.conn{1}.z{5}{3} for elevations on side 1, node 5, path 3
    %
    % See also: GRIDobj
    %
    % Author: Bastien Mathieux (b.mathieux@unistra.fr)
    % Date: May, 2024
    % Updated: July, 2025
    
    properties
        type        % string: Method used ('geometric' or 'flow')
        DEM         % GRIDobj: Digital elevation model
        x           % double: x-coordinates of the baseline
        y           % double: y-coordinates of the baseline
        conn        % cell{1x2}: Connection data for each side; each cell contains a struct with ix, x, y, z, d
        int         % cell{1x2}: Interpolated data for each side; structured similarly to conn
        stats       % struct: Statistics with fields x,y,z,d {num_nodes,1}, Z1,G1,Z2,G2,W,H (num_nodes,1)
    end
    
    methods
        function obj = Transect(DEM,x,y,varargin)
            % Constructor: Builds the Transect object
            
            %% Input parsing and validation
            
            p = inputParser;
            p.FunctionName = 'Transect';
            p.addRequired('DEM',@(x) isa(x,'GRIDobj'));
            p.addRequired('x',@isvector);
            p.addRequired('y',@isvector);
            p.addParameter('w',5,@(x) isscalar(x) && x > 0);
            p.addParameter('ite',1,@(x) isscalar(x) && x >= 1 && mod(x,1) == 0);
            p.addParameter('verbose',true,@islogical);
            p.addParameter('method','geometric',@(x) any(validatestring(x,{'geometric','flow'})));
            p.parse(DEM,x,y,varargin{:});
            
            DEM = p.Results.DEM;
            x = p.Results.x(:);
            y = p.Results.y(:);
            w = max(p.Results.w, 3);
            ite = max(p.Results.ite, 1);
            verbose = p.Results.verbose;
            method = p.Results.method;
            
            if numel(x) ~= numel(y)
                error('Transect:xAndYMismatch', 'x and y must have the same number of elements.');
            end
            
            if verbose
                if w < p.Results.w
                    warning('Transect:AdjustedWidth', 'Width too short (<3). Set to 3.');
                end
                if ite < p.Results.ite
                    warning('Transect:AdjustedIterations', 'Number of iterations must be at least 1. Set to 1.');
                end
            end
            
            if strcmp(method,'flow')
                if verbose
                    warning('Transect:FlowParams', 'Parameters ''w'' and ''ite'' are not used in ''flow'' method.');
                end
            end
            
            if numel(x) <= 3
                if verbose
                    warning('Transect:ShortBaseline', 'Baseline too short (<=3). Creating empty object.');
                end
                obj.x = [];
                obj.y = [];
                obj.conn = {};
                obj.int = {};
                return;
            end
            
            %% DEM padding (only for geometric method)
            
            if strcmp(method,'geometric')
                wi = (w + 1) * ite;
                [r, c] = size(DEM.Z);
                Z = nan(r + 2*wi, c + 2*wi);
                Z(wi+1:end-wi, wi+1:end-wi) = DEM.Z;
                DEM.Z = Z;
                DEM.refmat(3,1) = DEM.refmat(3,1) - wi * DEM.cellsize;
                DEM.refmat(3,2) = DEM.refmat(3,2) + wi * DEM.cellsize;
                DEM.size = size(Z);
            end
            
            obj.DEM = DEM;
            
            %% Baseline shortening
            
            [baseline_x, baseline_y] = shortpath(x, y, DEM);
            obj.x = baseline_x;
            obj.y = baseline_y;
            
            obj.type = method;
            
            obj.conn = cell(1,2);
            obj.int = cell(1,2);
            
            if strcmp(method,'geometric')
                
                %% Outline computation for sides
                
                outline_x = cell(1,2);
                outline_y = cell(1,2);
                [outline_x{1}, outline_y{1}, outline_x{2}, outline_y{2}] = obj.outline(DEM,1,baseline_x,baseline_y,baseline_x,baseline_y,w);
                
                %% Progress setup
                
                if verbose
                    PB = ProgressBar(200*ite, 'taskname','Extracting Transect object...', 'ui','cli');
                end
                
                %% Main construction loop
                
                for side = 1:2
                    
                    inner_x = cell(ite+1,1);
                    inner_y = cell(ite+1,1);
                    inner_ix = cell(ite+1,1);
                    end_x_cell = cell(ite+1,1);
                    end_y_cell = cell(ite+1,1);
                    end_ix = cell(ite+1,1);
                    source_ix = cell(ite,1);
                    target_ix = cell(ite,1);
                    
                    for iter = 1:ite
                        
                        if verbose
                            PB.count();
                        end
                        
                        if iter == 1
                            current_x = baseline_x;
                            current_y = baseline_y;
                            current_ix = coord2ind(DEM, current_x, current_y);
                            inner_x{iter} = current_x;
                            inner_y{iter} = current_y;
                            inner_ix{iter} = current_ix;
                            end_ix{iter} = current_ix;
                            
                            out_x = outline_x{side};
                            out_y = outline_y{side};
                        else
                            [current_x, current_y] = ind2coord(DEM, unique(target_ix{iter-1}, 'stable'));
                            [current_x, current_y] = shortpath(current_x, current_y, DEM);
                            
                            d1 = sqrt((baseline_x(1) - current_x(1))^2 + (baseline_y(1) - current_y(1))^2);
                            d2 = sqrt((baseline_x(1) - current_x(end))^2 + (baseline_y(1) - current_y(end))^2);
                            
                            if d2 < d1
                                current_x = flip(current_x);
                                current_y = flip(current_y);
                            end
                            
                            current_ix = coord2ind(DEM, current_x, current_y);
                            
                            inner_x{iter} = current_x;
                            inner_y{iter} = current_y;
                            inner_ix{iter} = current_ix;
                            
                            [out_x, out_y] = obj.outline(DEM, iter, current_x, current_y, baseline_x, baseline_y, w);
                        end
                        
                        target_ix{iter} = coord2ind(DEM, out_x, out_y);
                        
                        d1 = sqrt((current_x(1) - out_x(1))^2 + (current_y(1) - out_y(1))^2);
                        d2 = sqrt((current_x(1) - out_x(end))^2 + (current_y(1) - out_y(end))^2);
                        
                        if d2 <= d1
                            out_x = flip(out_x);
                            out_y = flip(out_y);
                        end
                        
                        %% End node extraction
                        
                        if verbose
                            PB.count();
                        end
                        
                        d1 = pdist2([out_x, out_y], [current_x(1), current_y(1)], 'euclidean');
                        id = find(d1 < w * DEM.cellsize + DEM.cellsize);
                        id = min(id):max(id);
                        
                        d1 = pdist2([out_x, out_y], [current_x(end), current_y(end)], 'euclidean');
                        idi = find(d1 < w * DEM.cellsize + DEM.cellsize);
                        id = sort([id min(idi):max(idi)]);
                        
                        if ~isempty(id)
                            if max(diff(id)) > 1
                                id = unique([1:id(1) id id(end):numel(out_x)], 'stable');
                            elseif abs(max(id) - 1) < abs(max(id) - numel(out_x))
                                id = unique([1:id(1) id], 'stable');
                            else
                                id = unique([id id(end):numel(out_x)], 'stable');
                            end
                        end
                        
                        end_x = out_x;
                        end_y = out_y;
                        out_x(id) = [];
                        out_y(id) = [];
                        
                        if numel(out_x) < numel(end_x) - numel(out_x)
                            dif = floor((numel(end_x) - numel(out_x)) / 2);
                            out_x = end_x(dif:end-dif);
                            out_y = end_y(dif:end-dif);
                        end
                        
                        inner_x{iter+1} = out_x;
                        inner_y{iter+1} = out_y;
                        end_x_cell{iter+1} = end_x;
                        end_y_cell{iter+1} = end_y;
                        inner_ix{iter+1} = coord2ind(DEM, inner_x{iter+1}, inner_y{iter+1});
                        end_ix{iter+1} = coord2ind(DEM, end_x_cell{iter+1}, end_y_cell{iter+1});
                        
                        %% Node pairing
                        
                        if verbose
                            PB.count();
                        end
                        
                        Iout = cell(length(current_x), 1);
                        for k = 1:numel(out_x)
                            [~, id] = pdist2([current_x, current_y], [out_x(k), out_y(k)], 'euclidean', 'Smallest', 1);
                            Iout{id} = [Iout{id}; k];
                        end
                        
                        missing_count = 0;
                        for k = 1:numel(current_x)
                            if k == numel(current_x) && isempty(Iout{k})
                                prev = k - missing_count - 1;
                                pi = max(Iout{prev});
                                si = max(Iout{prev});
                                dp = sqrt((out_x(pi) - current_x(prev+1:end-1)).^2 + (out_y(pi) - current_y(prev+1:end-1)).^2);
                                ds = sqrt((out_x(si) - current_x(prev+1:end-1)).^2 + (out_y(si) - current_y(prev+1:end-1)).^2);
                                for m = 1:missing_count
                                    if dp(m) >= ds(m)
                                        Iout{k-m} = [Iout{k-m}; pi];
                                    else
                                        Iout{k-m} = [Iout{k-m}; si];
                                    end
                                end
                                Iout{k} = Iout{k-1}(end);
                                missing_count = 0;
                            elseif isempty(Iout{k})
                                missing_count = missing_count + 1;
                            else
                                if k - missing_count == 1
                                    prev = 0;
                                    succ = k;
                                    pi = min(Iout{succ});
                                else
                                    prev = k - missing_count - 1;
                                    succ = k;
                                    pi = max(Iout{prev});
                                end
                                si = min(Iout{succ});
                                dp = sqrt((out_x(pi) - current_x(prev+1:succ-1)).^2 + (out_y(pi) - current_y(prev+1:succ-1)).^2);
                                ds = sqrt((out_x(si) - current_x(prev+1:succ-1)).^2 + (out_y(si) - current_y(prev+1:succ-1)).^2);
                                for m = 1:missing_count
                                    if dp(m) >= ds(m)
                                        Iout{k-m} = [Iout{k-m}; pi];
                                    else
                                        Iout{k-m} = [Iout{k-m}; si];
                                    end
                                end
                                missing_count = 0;
                            end
                        end
                        
                        s = [];
                        t = [];
                        for k = 1:numel(Iout)
                            s = [s repmat(k, 1, numel(Iout{k}))];
                            t = [t Iout{k}'];
                        end
                        
                        source_ix{iter} = inner_ix{iter}(s);
                        target_ix{iter} = inner_ix{iter+1}(t);
                        
                        %% Fan adjustment
                        
                        if verbose
                            PB.count();
                        end
                        
                        mod1 = 0; mod2 = 0; maxI = numel(unique(target_ix{iter}, 'stable'));
                        iii = 0; prev_maxf = Inf; unchanged_count = 0;
                        
                        while mod1 == 0 && mod2 == 0 && iii <= maxI
                            iii = iii + 1;
                            
                            [~, ~, ii] = unique(target_ix{iter}, 'stable');
                            f = diff([find([true, diff(ii') ~= 0, true])]);
                            maxf = max(f);
                            
                            if maxf <= prev_maxf
                                unchanged_count = unchanged_count + 1;
                            else
                                unchanged_count = 0;
                            end
                            
                            if unchanged_count >= 10
                                break;
                            end
                            
                            prev_maxf = maxf;
                            
                            k = 0;
                            ii = target_ix{iter};
                            
                            while maxf > 3 && k <= maxI
                                k = k + 1;
                                [~, ia, ib] = unique(ii, 'stable');
                                f = diff([find([true, diff(ii') ~= 0, true])]);
                                [maxf, m] = max(f);
                                ind = find(ib == ib(ia(m)));
                                if numel(ind) > 2
                                    ii1 = ii; ii2 = end_ix{iter+1};
                                    [~, ~, idi] = intersect(ii1(ind), ii2);
                                    if idi(1) == 1
                                        ii1(ind(1):ind(1) + floor(numel(ind)/3) - 1) = ii2(idi(1));
                                    else
                                        ii1(ind(1):ind(1) + floor(numel(ind)/3) - 1) = ii2(idi(1) - 1);
                                    end
                                    if numel(ii2) == idi(end)
                                        ii1(ind(end) - floor(numel(ind)/3) + 1:ind(end)) = ii2(idi(end));
                                    else
                                        ii1(ind(end) - floor(numel(ind)/3) + 1:ind(end)) = ii2(idi(end) + 1);
                                    end
                                    ii = ii1;
                                end
                            end
                            target_ix{iter} = ii;
                            
                            if k == 0
                                mod1 = 1;
                            else
                                mod1 = 0;
                            end
                            
                            [~, ~, ii] = unique(source_ix{iter}, 'stable');
                            f = diff([find([true, diff(ii') ~= 0, true])]);
                            maxf = max(f);
                            k = 0;
                            ii = source_ix{iter};
                            
                            while maxf > 3 && k <= maxI
                                k = k + 1;
                                [~, ia, ib] = unique(ii, 'stable');
                                f = diff([find([true, diff(ii') ~= 0, true])]);
                                [maxf, m] = max(f);
                                ind = find(ib == ib(ia(m)));
                                if numel(ind) > 2
                                    ii1 = ii; ii2 = end_ix{iter};
                                    [~, ~, idi] = intersect(ii1(ind), ii2);
                                    if idi(1) == 1
                                        ii1(ind(1):ind(1) + floor(numel(ind)/3) -1) = ii2(idi(1));
                                    else
                                        ii1(ind(1):ind(1) + floor(numel(ind)/3) -1) = ii2(idi(1)-1);
                                    end
                                    if numel(ii2) == idi(end)
                                        ii1(ind(end) - floor(numel(ind)/3) +1:ind(end)) = ii2(idi(end));
                                    else
                                        ii1(ind(end) - floor(numel(ind)/3) +1:ind(end)) = ii2(idi(end) +1);
                                    end
                                    ii = ii1;
                                end
                            end
                            source_ix{iter} = ii;
                            
                            if k == 0
                                mod2 = 1;
                            else
                            mod2 = 0;
                            end
                            
                            if iter > 1
                                [prev_target_x, prev_target_y] = ind2coord(DEM, target_ix{iter - 1});
                                [target_x, target_y] = ind2coord(DEM, target_ix{iter});
                                di = pdist2([prev_target_x(1), prev_target_y(1)], [target_x(1),target_y(1)]);
                                de = pdist2([prev_target_x(end), prev_target_y(end)], [target_x(1), target_y(1)]);
                                if de < di
                                    target_ix{iter} = flip(target_ix{iter});
                                    source_ix{iter} = flip(source_ix{iter});
                                end
                            end
                        end
                    end
                    
                    %% Graph assembly
                    
                    s = [];
                    t = [];
                    d = [];
                    
                    for iter = 1:ite
                        s = [s; source_ix{iter}];
                        t = [t; target_ix{iter}];
                        
                        [xs, ys] = ind2coord(DEM, source_ix{iter});
                        [xt, yt] = ind2coord(DEM, target_ix{iter});
                        
                        d = [d; sqrt((xs - xt).^2 + (ys - yt).^2)];
                    end
                    
                    G = digraph(s, t, d);
                    
                    %% Path computation
                    
                    [last_x, last_y] = ind2coord(DEM, unique(target_ix{end}, 'stable'));
                    
                    d1 = sqrt((baseline_x(1) - last_x(1))^2 + (baseline_y(1) - last_y(1))^2);
                    d2 = sqrt((baseline_x(1) - last_x(end))^2 + (baseline_y(1) - last_y(end))^2);
                    
                    if d2 < d1
                        last_x = flip(last_x);
                        last_y = flip(last_y);
                    end
                    
                    last_ix = coord2ind(DEM, last_x, last_y);
                    
                    inner_x{end} = last_x;
                    inner_y{end} = last_y;
                    inner_ix{end} = last_ix;
                    
                    baseline_nodes = inner_ix{1};
                    outline_nodes = inner_ix{end};
                    num_baseline = numel(baseline_nodes);
                    connected_nodes = cell(num_baseline, 1);
                    all_paths = cell(num_baseline, 1);
                    
                    for i = 1:num_baseline
                        src = baseline_nodes(i);
                        
                        connected_nodes{i} = intersect(outline_nodes, bfsearch(G, src));
                        
                        pS = arrayfun(@(tgt) allpaths(G, src, tgt), connected_nodes{i}, 'UniformOutput', false);
                        
                        cP = vertcat(pS{:});
                        
                        if ~isempty(cP)
                            pSi = cellfun(@mat2str, cP, 'UniformOutput', false);
                            upSi = unique(pSi);
                            
                            uP = cellfun(@str2num, upSi, 'UniformOutput', false);
                            
                            all_paths{i} = uP;
                        else
                            all_paths{i} = {};
                        end
                        
                        if verbose
                            if ismember(i, round(linspace(1, num_baseline, 96*ite)))
                                PB.count();
                            end
                        end
                    end
                    
                    %% Path data extraction
                    
                    path_d = cell(num_baseline, 1);
                    path_e = cell(num_baseline, 1);
                    path_x = cell(num_baseline, 1);
                    path_y = cell(num_baseline, 1);
                    path_id = cell(num_baseline, 1);
                    
                    path_di = cell(num_baseline, 1);
                    path_ei = cell(num_baseline, 1);
                    path_xi = cell(num_baseline, 1);
                    path_yi = cell(num_baseline, 1);
                    path_idi = cell(num_baseline, 1);
                    
                    for i = 1:num_baseline
                        num_paths = numel(all_paths{i});
                        dFS = cell(num_paths, 1);
                        eFS = cell(num_paths, 1);
                        xFS = cell(num_paths, 1);
                        yFS = cell(num_paths, 1);
                        ixFS = cell(num_paths, 1);
                        
                        for j = 1:num_paths
%                             p = all_paths{i}{j};
                            p = unique(all_paths{i}{j},'stable');
                            
                            [xP, yP] = ind2coord(DEM, p);
                            dd = sqrt((xP(2:end) - xP(1:end-1)).^2 + (yP(2:end) - yP(1:end-1)).^2);
                            
                            dFS{j} = [0; cumsum(dd)];
                            
                            eFS{j} = DEM.Z(p);
                            
                            xFS{j} = xP;
                            yFS{j} = yP;
                            ixFS{j} = p;
                        end
                        
                        path_d{i} = dFS;
                        path_e{i} = eFS;
                        path_x{i} = xFS;
                        path_y{i} = yFS;
                        path_id{i} = ixFS;
                    end
                    
                    %% Path interpolation
                    
                    for i = 1:num_baseline
                        for j = 1:numel(path_x{i})
                            xii = path_x{i}{j};
                            yii = path_y{i}{j};
                            dii = path_d{i}{j};
                            
                            intD = 0:DEM.cellsize/2:dii(end);
                            nx = interp1(dii, xii, intD, 'linear');
                            ny = interp1(dii, yii, intD, 'linear');
                            nd = [0 sqrt(diff(nx).^2 + diff(ny).^2)];
                            
                            nid = coord2ind(DEM, nx, ny);
                            [nid, idx] = unique(nid, 'stable');
                            [nx, ny] = ind2coord(DEM, nid);
                            nd = nd(idx);
                            ne = DEM.Z(nid);
                            
                            path_di{i}{j} = cumsum(nd);
                            path_ei{i}{j} = ne;
                            path_xi{i}{j} = nx;
                            path_yi{i}{j} = ny;
                            path_idi{i}{j} = nid;
                        end
                    end
                    
                    %% Data storage
                    
                    conn_struct = struct('ix', path_id, 'x', path_x, 'y', path_y, 'z', path_e, 'd', path_d);
                    int_struct = struct('ix', path_idi, 'x', path_xi, 'y', path_yi, 'z', path_ei, 'd', path_di);
                    
                    obj.conn{side} = conn_struct;
                    obj.int{side} = int_struct;
                end
                
                else % 'flow' method

                    if verbose
                        PB = ProgressBar(numel(baseline_x), 'taskname','Extracting Transect object (flow mode)...', 'ui','cli');
                    end
                
                    FD = FLOWobj(DEM,'preprocess','carve');
                    x = baseline_x; y = baseline_y; sz = size(DEM.Z);
                    lout = cell(numel(x),1); rout = cell(numel(x),1);
                
                    for i = 1:numel(x)
                        if verbose
                            PB.count();
                        end
                        ix = coord2ind(DEM, x(i), y(i));
                        BM = drainagebasins(FD, ix);
                        B = bwboundaries(BM.Z, 8, 'noholes');
                        ixo = sub2ind(size(BM.Z), B{1}(:,1), B{1}(:,2));
                        [xo, yo] = ind2coord(DEM, ixo);
                
                        [~, si] = min(hypot(xo - x(i), yo - y(i)));
                        ixo = circshift(ixo, -(si-1), 1);
                        mid = round(size(ixo,1)/2);
                        lout{i} = ixo(1:mid);
                        rout{i} = ixo(mid+1:end);
                    end
                
                    % Save results to the object
                    conn_struct1 = repmat(struct('ix', {}, 'x', {}, 'y', {}, 'z', {}, 'd', {}), 1, numel(x));
                    conn_struct2 = repmat(struct('ix', {}, 'x', {}, 'y', {}, 'z', {}, 'd', {}), 1, numel(x));
                    for i = 1:numel(x)
                        % Left
                        id = unique(lout{i},'stable');
                        [xl,yl] = ind2coord(DEM,id);
                        zl = DEM.Z(id);
                        dl = [0; cumsum(sqrt(diff(xl).^2 + diff(yl).^2))];
                        conn_struct1(i).ix={id}; conn_struct1(i).x={xl}; conn_struct1(i).y={yl};
                        conn_struct1(i).z={zl};  conn_struct1(i).d={dl};
                        
                        % Right
                        id2 = unique(flip(rout{i}),'stable');
                        [xr,yr] = ind2coord(DEM,id2);
                        zr = DEM.Z(id2);
                        dr = [0; cumsum(sqrt(diff(xr).^2 + diff(yr).^2))];
                        conn_struct2(i).ix={id2}; conn_struct2(i).x={xr}; conn_struct2(i).y={yr};
                        conn_struct2(i).z={zr};  conn_struct2(i).d={dr};

                    end
                
                    obj.conn{1} = conn_struct1;
                    obj.conn{2} = conn_struct2;
                    obj.int{1} = conn_struct1;
                    obj.int{2} = conn_struct2;
            end
            
            %% Compute statistics
            
            xii = cell(2,1);
            yii = cell(2,1);
            zii = cell(2,1);
            dii = cell(2,1);
            
            for i1 = 1:2
                n2 = numel(obj.x);
                xii{i1} = cell(n2, 1);
                yii{i1} = cell(n2, 1);
                zii{i1} = cell(n2, 1);
                dii{i1} = cell(n2, 1);
            
                for i2 = 1:n2
                    if ~isempty(obj.int{i1}(i2).x)
                        midx = ceil(numel(obj.int{i1}(i2).x)/2);
                        xii{i1}{i2} = obj.int{i1}(i2).x{midx}(:);
                        yii{i1}{i2} = obj.int{i1}(i2).y{midx}(:);
                        zii{i1}{i2} = obj.int{i1}(i2).z{midx}(:);
                        dii{i1}{i2} = obj.int{i1}(i2).d{midx}(:);
                    end
                end
            end

            obj.stats.x = cell(numel(obj.x),1);
            obj.stats.y = cell(numel(obj.x),1);
            obj.stats.z = cell(numel(obj.x),1);
            obj.stats.d = cell(numel(obj.x),1);

            obj.stats.Z1 = nan(numel(obj.x),1);
            obj.stats.Z2 = nan(numel(obj.x),1);
            obj.stats.G1 = nan(numel(obj.x),1);
            obj.stats.G2 = nan(numel(obj.x),1);
            obj.stats.W = nan(numel(obj.x),1);
            obj.stats.H = nan(numel(obj.x),1);

            
            for i2 = 1:numel(obj.x)
            
                try                    
                    obj.stats.x{i2} = [flip(xii{1}{i2}); xii{2}{i2}(2:end)];
                    obj.stats.y{i2} = [flip(yii{1}{i2}); yii{2}{i2}(2:end)];
                    obj.stats.z{i2} = [flip(zii{1}{i2}); zii{2}{i2}(2:end)];
                    obj.stats.d{i2} = [flip(dii{1}{i2} * (-1)); dii{2}{i2}(2:end)];
                catch
                    obj.stats.x{i2} = [];
                    obj.stats.y{i2} = [];
                    obj.stats.z{i2} = [];
                    obj.stats.d{i2} = [];
                end
            end
            
            for i1 = 1:numel(obj.x)
            
                try
                    id1 = find(obj.stats.d{i1} < 0);
                    id2 = find(obj.stats.d{i1} > 0);
                    Z1 = obj.stats.z{i1}(id1);
                    Z2 = obj.stats.z{i1}(id2);
                    D1 = obj.stats.d{i1}(id1);
                    D2 = obj.stats.d{i1}(id2);
                catch
                end
            
                try
                    obj.stats.Z1(i1) = max(Z1);
                    G1 = abs(atan(diff(Z1) ./ diff(abs(D1))) * (180/pi)); 
                    obj.stats.G1(i1) = mean(G1);
                catch
                end
                     
                try
                    obj.stats.Z2(i1) = max(Z2);
                    G2 = abs(atan(diff(Z2) ./ diff(abs(D2))) * (180/pi)); 
                    obj.stats.G2(i1) = mean(G2);
                catch
                end
            
                try
                    obj.stats.W(i1) = obj.stats.d{i1}(1) * (-1) + obj.stats.d{i1}(end);
                    obj.stats.H(i1) = mean([obj.stats.Z1(i1) obj.stats.Z2(i1)]) - min(obj.stats.z{i1});
                catch
                end
            end
        end
    end
    
    methods (Access = private, Static)
        function varargout = outline(DEM, i2, x, y, xa, ya, w)
            %outline Extract outline nodes for a layer
            
            A = false(DEM.size);
            ix = coord2ind(DEM, x, y);
            A(ix) = true;
            
            [B, ~] = bwdist(A, 'euclidean');
            mask = B <= w;
            
            bw = bwboundaries(mask, 8);
            [xo, yo] = sub2coord(DEM, bw{1}(:, 1), bw{1}(:, 2));
            D = pdist2([xo, yo], [xo, yo], 'euclidean');
            [~, max_idx] = max(D(:));
            [ii1, ii2] = ind2sub(size(D), max_idx);
            d1 = pdist2([x(1), y(1)], [xo(ii1), yo(ii1)], 'euclidean');
            d2 = pdist2([x(1), y(1)], [xo(ii2), yo(ii2)], 'euclidean');
            if d2 > d1
                ii1 = ii2;
            end
            
            xo = [xo(ii1:end); xo(1:ii1-1)];
            yo = [yo(ii1:end); yo(1:ii1-1)];
            
            n2 = ceil(numel(xo) / 2);
            xo1 = xo(1:n2);
            yo1 = yo(1:n2);
            xo2 = flip(xo(n2+1:end));
            yo2 = flip(yo(n2+1:end));
            
            d1 = median(min(pdist2([xo1, yo1], [xa, ya], 'euclidean'), [], 2));
            d2 = median(min(pdist2([xo2, yo2], [xa, ya], 'euclidean'), [], 2));
            
            if i2 == 1
                varargout = {xo1, yo1, xo2, yo2};
            else
                if d1 > d2
                    varargout = {xo1, yo1};
                else
                    varargout = {xo2, yo2};
                end
            end
        end
    end
end