classdef TRANSECT
    %TRANSECT object 
    %
    % Syntax
    %
    %     obj = TRANSECT(DEM,x,y)
    %     obj = TRANSECT(DEM,S)
    %
    % Description
    %
    %     TRANSECT creates a transect object along a user-defined baseline defined by
    %     coordinates x and y, or by a single STREAMobj S, on a digital elevation
    %     model DEM (GRIDobj). The algorithm generates transects on both sides of
    %     the baseline according to the selected connection mode, enabling consistent
    %     extraction and alignment of geomorphic metrics across multiple landscape
    %     domains.
    %
    %     Two connection modes are available:
    %         (1) 'geometric' – constructs sub-perpendicular transects to the baseline at
    %             regular intervals using Euclidean distance, preserving the native geometry
    %             and node spacing of curved baselines.
    %         (2) 'flow' – traces pathways derived from DEM-based flow directions,
    %             connecting nodes according to natural topographic gradients and following
    %             hydrological or sediment transport routes.
    %
    %     The DEM is padded only in 'geometric' mode and only when the w*ite buffer
    %     around baseline nodes exceeds the DEM limits. In that mode, the parameters
    %     'w' and 'ite' control the geometric construction. In 'flow' mode, 'w' and
    %     'ite' are ignored.
    %
    %     Path data (x, y, z, d, ix) is organized hierarchically as:
    %         T.int{i1}(i2).x{i3}, T.conn{i1}(i2).z{i3}, etc.
    %
    %     Indexing convention:
    %         i1 = side of baseline (1,2)
    %         i2 = baseline node index
    %         i3 = connected path index linking back to baseline node i2
    %
    %     Development note
    %     ----------------
    %     TRANSECT is currently designed for a single continuous river line or
    %     a single continuous x/y baseline. Full STREAMobj networks must be split
    %     into individual branches before creating TRANSECT objects.
    %
    % Input arguments
    %
    %     DEM        GRIDobj containing the digital elevation model
    %     x,y        coordinate vectors defining the baseline
    %     S          single STREAMobj defining the baseline
    %
    % Parameter name/value pairs
    %
    %     'method'   string: 'geometric' (default) or 'flow'
    %     'w'        width parameter (cells) used only in 'geometric' mode. Default = 5. Must be >3.
    %     'ite'      number of iterations used only in 'geometric' mode. Default = 1. Must be >=1.
    %     'verbose'  logical scalar. Display progress and warnings. Default = true.
    %
    % Properties
    %
    %     type        string: Method used ('geometric' or 'flow')
    %     DEM         GRIDobj: Digital elevation model
    %     base        struct: baseline information with fields:
    %                   typ: 'xy' or 'STREAMobj'
    %                   x,y: final internal baseline coordinates used by TRANSECT
    %                   ix: DEM indices of the final internal baseline
    %                   dist: cumulative distance along the final internal baseline
    %                   src: positional index along the original input baseline
    %                   obj: original STREAMobj when typ is STREAMobj
    %     conn        cell{1x2}: Connection data for each side
    %     int         cell{1x2}: Interpolated data for each side
    %     stats       [] or struct: Statistics created by stats. Default = []
    %     pair        [] or struct: Paired object created by pairing. Default = []
    %
    % See also: GRIDobj, FLOWobj, STREAMobj

    properties
        type
        DEM
        base
        conn
        int
        stats
        pair
    end

    methods
        function obj = TRANSECT(DEM,varargin)

            if ~isa(DEM,'GRIDobj')
                error('TRANSECT:InvalidDEM','DEM must be a GRIDobj.');
            end

            [x,y,src,base,varargin] = TRANSECT.parsebase(DEM,varargin{:});

            p = inputParser;
            p.FunctionName = 'TRANSECT';
            p.addParameter('method','geometric',@(x) any(validatestring(x,{'geometric','flow'})));
            p.addParameter('w',5,@(x) isscalar(x) && x > 0);
            p.addParameter('ite',1,@(x) isscalar(x) && x >= 1 && mod(x,1) == 0);
            p.addParameter('verbose',true,@islogical);
            p.parse(varargin{:});

            method = p.Results.method;
            w = max(p.Results.w,3);
            ite = max(p.Results.ite,1);
            verbose = p.Results.verbose;

            if numel(x) ~= numel(y)
                error('TRANSECT:xAndYMismatch','x and y must have the same number of elements.');
            end

            if verbose
                if w < p.Results.w
                    warning('TRANSECT:AdjustedWidth','Width too short (<3). Set to 3 (geometric mode only).');
                end
                if ite < p.Results.ite
                    warning('TRANSECT:AdjustedIterations','Number of iterations must be at least 1. Set to 1 (geometric mode only).');
                end
            end

            if strcmp(method,'flow') && verbose && ...
                    (~ismember('w',p.UsingDefaults) || ~ismember('ite',p.UsingDefaults))
                warning('TRANSECT:FlowParams','Parameters ''w'' and ''ite'' are only used in ''geometric'' mode and are ignored in ''flow'' mode.');
            end

            if numel(x) <= 3
                if verbose
                    warning('TRANSECT:ShortBaseline','Baseline too short (<=3). Creating empty object.');
                end
                ix = coord2ind(DEM,x,y);
                base.x = x(:);
                base.y = y(:);
                base.ix = ix(:);
                base.dist = [0; cumsum(hypot(diff(x(:)),diff(y(:))))];
                base.src = src(:);

                obj.type = method;
                obj.DEM = DEM;
                obj.base = base;
                obj.conn = {};
                obj.int = {};
                obj.stats = [];
                obj.pair = [];
                return
            end

            DEMi = DEM;
            padon = false;

            if strcmp(method,'geometric')
                wi = ceil(w*ite);
                ixp = coord2ind(DEM,x,y);
                ixp = ixp(~isnan(ixp));

                if ~isempty(ixp)
                    [r,c] = ind2sub(DEM.size,ixp);
                    padon = any(r-wi < 1 | r+wi > DEM.size(1) | c-wi < 1 | c+wi > DEM.size(2));
                end

                if padon
                    DEM = pad(DEM,wi,nan);
                end
            end

            obj.DEM = DEM;

            [baseline_x,baseline_y] = shortpath(x,y,DEM);

            d1 = hypot(baseline_x(1)-x(1),baseline_y(1)-y(1)) + ...
                 hypot(baseline_x(end)-x(end),baseline_y(end)-y(end));
            d2 = hypot(baseline_x(end)-x(1),baseline_y(end)-y(1)) + ...
                 hypot(baseline_x(1)-x(end),baseline_y(1)-y(end));

            if d2 < d1
                baseline_x = flipud(baseline_x);
                baseline_y = flipud(baseline_y);
            end

            ix = coord2ind(DEM,baseline_x,baseline_y);
            dst = [0; cumsum(hypot(diff(baseline_x(:)),diff(baseline_y(:))))];
            srcb = TRANSECT.mapsrc(x,y,src,baseline_x,baseline_y);

            base.x = baseline_x(:);
            base.y = baseline_y(:);
            base.ix = ix(:);
            base.dist = dst(:);
            base.src = srcb(:);

            obj.base = base;
            obj.type = method;
            obj.conn = cell(1,2);
            obj.int = cell(1,2);
            obj.stats = [];
            obj.pair = [];

            if strcmp(method,'geometric')

                outline_x = cell(1,2);
                outline_y = cell(1,2);
                [outline_x{1},outline_y{1},outline_x{2},outline_y{2}] = obj.outline(DEM,1,baseline_x,baseline_y,baseline_x,baseline_y,w);

                if verbose
                    np = 2*(4*ite + numel(baseline_x));
                    PB = ProgressBar(np,'taskname','Extracting TRANSECT object...','ui','cli');
                end

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
                            current_ix = coord2ind(DEM,current_x,current_y);
                            inner_x{iter} = current_x;
                            inner_y{iter} = current_y;
                            inner_ix{iter} = current_ix;
                            end_ix{iter} = current_ix;

                            out_x = outline_x{side};
                            out_y = outline_y{side};
                        else
                            [current_x,current_y] = ind2coord(DEM,unique(target_ix{iter-1},'stable'));
                            [current_x,current_y] = shortpath(current_x,current_y,DEM);

                            d1 = hypot(baseline_x(1)-current_x(1),baseline_y(1)-current_y(1));
                            d2 = hypot(baseline_x(1)-current_x(end),baseline_y(1)-current_y(end));

                            if d2 < d1
                                current_x = flip(current_x);
                                current_y = flip(current_y);
                            end

                            current_ix = coord2ind(DEM,current_x,current_y);

                            inner_x{iter} = current_x;
                            inner_y{iter} = current_y;
                            inner_ix{iter} = current_ix;

                            [out_x,out_y] = obj.outline(DEM,iter,current_x,current_y,baseline_x,baseline_y,w);
                        end

                        target_ix{iter} = coord2ind(DEM,out_x,out_y);

                        d1 = hypot(current_x(1)-out_x(1),current_y(1)-out_y(1));
                        d2 = hypot(current_x(1)-out_x(end),current_y(1)-out_y(end));

                        if d2 <= d1
                            out_x = flip(out_x);
                            out_y = flip(out_y);
                        end

                        if verbose
                            PB.count();
                        end

                        d1 = pdist2([out_x,out_y],[current_x(1),current_y(1)],'euclidean');
                        id = find(d1 < w*DEM.cellsize + DEM.cellsize);
                        id = min(id):max(id);

                        d1 = pdist2([out_x,out_y],[current_x(end),current_y(end)],'euclidean');
                        idi = find(d1 < w*DEM.cellsize + DEM.cellsize);
                        id = sort([id min(idi):max(idi)]);

                        if ~isempty(id)
                            if max(diff(id)) > 1
                                id = unique([1:id(1) id id(end):numel(out_x)],'stable');
                            elseif abs(max(id)-1) < abs(max(id)-numel(out_x))
                                id = unique([1:id(1) id],'stable');
                            else
                                id = unique([id id(end):numel(out_x)],'stable');
                            end
                        end

                        end_x = out_x;
                        end_y = out_y;
                        out_x(id) = [];
                        out_y(id) = [];

                        if numel(out_x) < numel(end_x)-numel(out_x)
                            dif = floor((numel(end_x)-numel(out_x))/2);
                            out_x = end_x(dif:end-dif);
                            out_y = end_y(dif:end-dif);
                        end

                        inner_x{iter+1} = out_x;
                        inner_y{iter+1} = out_y;
                        end_x_cell{iter+1} = end_x;
                        end_y_cell{iter+1} = end_y;
                        inner_ix{iter+1} = coord2ind(DEM,inner_x{iter+1},inner_y{iter+1});
                        end_ix{iter+1} = coord2ind(DEM,end_x_cell{iter+1},end_y_cell{iter+1});

                        if verbose
                            PB.count();
                        end

                        Iout = cell(length(current_x),1);
                        for k = 1:numel(out_x)
                            [~,id] = pdist2([current_x,current_y],[out_x(k),out_y(k)],'euclidean','Smallest',1);
                            Iout{id} = [Iout{id}; k];
                        end

                        missing_count = 0;
                        for k = 1:numel(current_x)
                            if k == numel(current_x) && isempty(Iout{k})
                                prev = k - missing_count - 1;
                                pi = max(Iout{prev});
                                si = max(Iout{prev});
                                dp = hypot(out_x(pi)-current_x(prev+1:end-1),out_y(pi)-current_y(prev+1:end-1));
                                ds = hypot(out_x(si)-current_x(prev+1:end-1),out_y(si)-current_y(prev+1:end-1));
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
                                dp = hypot(out_x(pi)-current_x(prev+1:succ-1),out_y(pi)-current_y(prev+1:succ-1));
                                ds = hypot(out_x(si)-current_x(prev+1:succ-1),out_y(si)-current_y(prev+1:succ-1));
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
                            s = [s repmat(k,1,numel(Iout{k}))];
                            t = [t Iout{k}'];
                        end

                        source_ix{iter} = inner_ix{iter}(s);
                        target_ix{iter} = inner_ix{iter+1}(t);

                        if verbose
                            PB.count();
                        end

                        mod1 = 0; mod2 = 0; maxI = numel(unique(target_ix{iter},'stable'));
                        iii = 0; prev_maxf = Inf; unchanged_count = 0;

                        while mod1 == 0 && mod2 == 0 && iii <= maxI
                            iii = iii + 1;

                            [~,~,ii] = unique(target_ix{iter},'stable');
                            f = diff([find([true,diff(ii') ~= 0,true])]);
                            maxf = max(f);

                            if maxf <= prev_maxf
                                unchanged_count = unchanged_count + 1;
                            else
                                unchanged_count = 0;
                            end

                            if unchanged_count >= 10
                                break
                            end

                            prev_maxf = maxf;

                            k = 0;
                            ii = target_ix{iter};

                            while maxf > 3 && k <= maxI
                                k = k + 1;
                                [~,ia,ib] = unique(ii,'stable');
                                f = diff([find([true,diff(ii(:)') ~= 0,true])]);
                                [maxf,m] = max(f);
                                ind = find(ib == ib(ia(m)));
                                if numel(ind) > 2
                                    ii1 = ii; ii2 = end_ix{iter+1};
                                    [~,~,idi] = intersect(ii1(ind),ii2,'stable');
                                    if idi(1) == 1
                                        ii1(ind(1):ind(1)+floor(numel(ind)/3)-1) = ii2(idi(1));
                                    else
                                        ii1(ind(1):ind(1)+floor(numel(ind)/3)-1) = ii2(idi(1)-1);
                                    end
                                    if numel(ii2) == idi(end)
                                        ii1(ind(end)-floor(numel(ind)/3)+1:ind(end)) = ii2(idi(end));
                                    else
                                        ii1(ind(end)-floor(numel(ind)/3)+1:ind(end)) = ii2(idi(end)+1);
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

                            [~,~,ii] = unique(source_ix{iter},'stable');
                            f = diff([find([true,diff(ii') ~= 0,true])]);
                            maxf = max(f);
                            k = 0;
                            ii = source_ix{iter};

                            while maxf > 3 && k <= maxI
                                k = k + 1;
                                [~,ia,ib] = unique(ii,'stable');
                                f = diff([find([true,diff(ii') ~= 0,true])]);
                                [maxf,m] = max(f);
                                ind = find(ib == ib(ia(m)));
                                if numel(ind) > 2
                                    ii1 = ii; ii2 = end_ix{iter};
                                    [~,~,idi] = intersect(ii1(ind),ii2);
                                    if idi(1) == 1
                                        ii1(ind(1):ind(1)+floor(numel(ind)/3)-1) = ii2(idi(1));
                                    else
                                        ii1(ind(1):ind(1)+floor(numel(ind)/3)-1) = ii2(idi(1)-1);
                                    end
                                    if numel(ii2) == idi(end)
                                        ii1(ind(end)-floor(numel(ind)/3)+1:ind(end)) = ii2(idi(end));
                                    else
                                        ii1(ind(end)-floor(numel(ind)/3)+1:ind(end)) = ii2(idi(end)+1);
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
                                [prev_target_x,prev_target_y] = ind2coord(DEM,target_ix{iter-1});
                                [target_x,target_y] = ind2coord(DEM,target_ix{iter});
                                di = pdist2([prev_target_x(1),prev_target_y(1)],[target_x(1),target_y(1)]);
                                de = pdist2([prev_target_x(end),prev_target_y(end)],[target_x(1),target_y(1)]);
                                if de < di
                                    target_ix{iter} = flip(target_ix{iter});
                                    source_ix{iter} = flip(source_ix{iter});
                                end
                            end
                        end
                    end

                    s = [];
                    t = [];
                    d = [];

                    for iter = 1:ite
                        s = [s; source_ix{iter}];
                        t = [t; target_ix{iter}];

                        [xs,ys] = ind2coord(DEM,source_ix{iter});
                        [xt,yt] = ind2coord(DEM,target_ix{iter});

                        d = [d; hypot(xs-xt,ys-yt)];
                    end

                    G = digraph(s,t,d);

                    [last_x,last_y] = ind2coord(DEM,unique(target_ix{end},'stable'));

                    d1 = hypot(baseline_x(1)-last_x(1),baseline_y(1)-last_y(1));
                    d2 = hypot(baseline_x(1)-last_x(end),baseline_y(1)-last_y(end));

                    if d2 < d1
                        last_x = flip(last_x);
                        last_y = flip(last_y);
                    end

                    last_ix = coord2ind(DEM,last_x,last_y);

                    inner_x{end} = last_x;
                    inner_y{end} = last_y;
                    inner_ix{end} = last_ix;

                    baseline_nodes = inner_ix{1};
                    outline_nodes = inner_ix{end};
                    num_baseline = numel(baseline_nodes);
                    nb = max(5,ceil(w));
                    bm = false(DEM.size);
                    bm(baseline_nodes) = true;
                    connected_nodes = cell(num_baseline,1);
                    all_paths = cell(num_baseline,1);

                    for i = 1:num_baseline
                        srcn = baseline_nodes(i);

                        connected_nodes{i} = intersect(outline_nodes,bfsearch(G,srcn));
                        pS = arrayfun(@(tgt) allpaths(G,srcn,tgt),connected_nodes{i},'UniformOutput',false);
                        cP = vertcat(pS{:});

                        if ~isempty(cP)
                            pSi = cellfun(@mat2str,cP,'UniformOutput',false);
                            upSi = unique(pSi);
                            all_paths{i} = cellfun(@str2num,upSi,'UniformOutput',false);
                        else
                            all_paths{i} = {};
                        end

                        if verbose
                            PB.count();
                        end
                    end

                    path_d = cell(num_baseline,1);
                    path_e = cell(num_baseline,1);
                    path_x = cell(num_baseline,1);
                    path_y = cell(num_baseline,1);
                    path_id = cell(num_baseline,1);

                    path_di = cell(num_baseline,1);
                    path_ei = cell(num_baseline,1);
                    path_xi = cell(num_baseline,1);
                    path_yi = cell(num_baseline,1);
                    path_idi = cell(num_baseline,1);

                    for i = 1:num_baseline
                        num_paths = numel(all_paths{i});
                        dFS = cell(num_paths,1);
                        eFS = cell(num_paths,1);
                        xFS = cell(num_paths,1);
                        yFS = cell(num_paths,1);
                        ixFS = cell(num_paths,1);

                        for j = 1:num_paths
                            p = unique(all_paths{i}{j},'stable');
                            ii = max(1,i-nb):min(num_baseline,i+nb);
                            msk = bm;
                            msk(baseline_nodes(ii)) = false;
                            q = find(msk(p),1);

                            if ~isempty(q)
                                p = p(1:max(1,q-1));
                            end

                            if numel(p) < 2
                                continue
                            end

                            [xP,yP] = ind2coord(DEM,p);
                            dd = hypot(diff(xP),diff(yP));

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

                    for i = 1:num_baseline
                        for j = 1:numel(path_x{i})
                            xii = path_x{i}{j};
                            yii = path_y{i}{j};
                            dii = path_d{i}{j};

                            if numel(xii) < 2
                                continue
                            end

                            intD = 0:DEM.cellsize/2:dii(end);
                            nx = interp1(dii,xii,intD,'linear');
                            ny = interp1(dii,yii,intD,'linear');
                            nd = [0 sqrt(diff(nx).^2 + diff(ny).^2)];

                            nid = coord2ind(DEM,nx,ny);
                            [nid,idx] = unique(nid,'stable');
                            [nx,ny] = ind2coord(DEM,nid);
                            nd = nd(idx);
                            ne = DEM.Z(nid);

                            path_di{i}{j} = cumsum(nd);
                            path_ei{i}{j} = ne;
                            path_xi{i}{j} = nx;
                            path_yi{i}{j} = ny;
                            path_idi{i}{j} = nid;
                        end
                    end

                    obj.conn{side} = struct('ix',path_id,'x',path_x,'y',path_y,'z',path_e,'d',path_d);
                    obj.int{side} = struct('ix',path_idi,'x',path_xi,'y',path_yi,'z',path_ei,'d',path_di);
                end

            else

                if verbose
                    PB = ProgressBar(numel(baseline_x),'taskname','Extracting Transect object (flow mode)...','ui','cli');
                end

                FD = FLOWobj(DEM,'preprocess','carve');
                x = baseline_x;
                y = baseline_y;
                lout = cell(numel(x),1);
                rout = cell(numel(x),1);

                for i = 1:numel(x)
                    if verbose
                        PB.count();
                    end

                    ix = coord2ind(DEM,x(i),y(i));
                    BM = drainagebasins(FD,ix);
                    B = bwboundaries(BM.Z,8,'noholes');
                    ixo = sub2ind(size(BM.Z),B{1}(:,1),B{1}(:,2));
                    [xo,yo] = ind2coord(DEM,ixo);

                    [~,si] = min(hypot(xo-x(i),yo-y(i)));
                    ixo = circshift(ixo,-(si-1),1);
                    mid = round(size(ixo,1)/2);
                    lout{i} = ixo(1:mid);
                    rout{i} = ixo(mid+1:end);
                end

                conn_struct1 = repmat(struct('ix',{},'x',{},'y',{},'z',{},'d',{}),1,numel(x));
                conn_struct2 = repmat(struct('ix',{},'x',{},'y',{},'z',{},'d',{}),1,numel(x));

                for i = 1:numel(x)
                    id = unique(lout{i},'stable');
                    [xl,yl] = ind2coord(DEM,id);
                    zl = DEM.Z(id);
                    dl = [0; cumsum(hypot(diff(xl),diff(yl)))];
                    conn_struct1(i).ix = {id};
                    conn_struct1(i).x = {xl};
                    conn_struct1(i).y = {yl};
                    conn_struct1(i).z = {zl};
                    conn_struct1(i).d = {dl};

                    id2 = unique(flip(rout{i}),'stable');
                    [xr,yr] = ind2coord(DEM,id2);
                    zr = DEM.Z(id2);
                    dr = [0; cumsum(hypot(diff(xr),diff(yr)))];
                    conn_struct2(i).ix = {id2};
                    conn_struct2(i).x = {xr};
                    conn_struct2(i).y = {yr};
                    conn_struct2(i).z = {zr};
                    conn_struct2(i).d = {dr};
                end

                obj.conn{1} = conn_struct1;
                obj.conn{2} = conn_struct2;
                obj.int{1} = conn_struct1;
                obj.int{2} = conn_struct2;
            end

            if strcmp(method,'geometric') && padon
                % if verbose
                %     disp('Unpadding...')
                % end

                obj = TRANSECT.unpad(obj,DEMi);

                % if verbose
                %     disp('Finished...')
                % end
            end
        end
    end

    methods (Access = private, Static)
        function [x,y,src,base,varargin] = parsebase(DEM,varargin)

            if isempty(varargin)
                error('TRANSECT:MissingBaseline','Baseline must be x/y or STREAMobj.');
            end

            b = varargin{1};

            if isa(b,'STREAMobj')
                x = b.x(:);
                y = b.y(:);
                src = (1:numel(x))';
                base = struct('typ','STREAMobj','obj',b);
                varargin(1) = [];

            elseif numel(varargin) >= 2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
                x = varargin{1}(:);
                y = varargin{2}(:);
                src = (1:numel(x))';
                base = struct('typ','xy');
                varargin(1:2) = [];

            else
                error('TRANSECT:InvalidBaseline','Baseline must be x/y or STREAMobj.');
            end

            ok = isfinite(x) & isfinite(y);
            x = x(ok);
            y = y(ok);
            src = src(ok);

            ix = coord2ind(DEM,x,y);
            ok = ~isnan(ix);
            x = x(ok);
            y = y(ok);
            src = src(ok);
        end

        function srcb = mapsrc(x,y,src,xb,yb)

            x = x(:); y = y(:); src = src(:);
            xb = xb(:); yb = yb(:);

            dx = diff(x);
            dy = diff(y);
            ds = hypot(dx,dy);
            s = [0; cumsum(ds)];

            ok = ds > 0;

            if ~any(ok)
                srcb = src(ones(numel(xb),1));
                return
            end

            x1 = x(1:end-1); x1 = x1(ok);
            y1 = y(1:end-1); y1 = y1(ok);
            dx = dx(ok);
            dy = dy(ok);
            s0 = s(1:end-1); s0 = s0(ok);
            ds = ds(ok);

            u = ((xb.'-x1).*dx + (yb.'-y1).*dy)./(ds.^2);
            u = max(0,min(1,u));

            xp = x1 + u.*dx;
            yp = y1 + u.*dy;

            [~,j] = min((xb.'-xp).^2 + (yb.'-yp).^2,[],1);
            j = j(:);
            k = (1:numel(xb)).';

            q = s0(j) + u(sub2ind(size(u),j,k)).*ds(j);
            q = max(0,min(s(end),cummax(q)));

            id = interp1(s,(1:numel(src))',q,'nearest','extrap');
            id = max(1,min(numel(src),id));

            srcb = src(id(:));
        end

        function varargout = outline(DEM,i2,x,y,xa,ya,w)

            A = false(DEM.size);
            ix = coord2ind(DEM,x,y);
            A(ix) = true;

            [B,~] = bwdist(A,'euclidean');
            mask = B <= w;

            bw = bwboundaries(mask,8);
            [xo,yo] = sub2coord(DEM,bw{1}(:,1),bw{1}(:,2));

            D = pdist2([xo,yo],[xo,yo],'euclidean');
            [~,max_idx] = max(D(:));
            [ii1,ii2] = ind2sub(size(D),max_idx);

            d1 = pdist2([x(1),y(1)],[xo(ii1),yo(ii1)],'euclidean');
            d2 = pdist2([x(1),y(1)],[xo(ii2),yo(ii2)],'euclidean');

            if d2 > d1
                ii1 = ii2;
            end

            xo = [xo(ii1:end); xo(1:ii1-1)];
            yo = [yo(ii1:end); yo(1:ii1-1)];

            n2 = ceil(numel(xo)/2);
            xo1 = xo(1:n2);
            yo1 = yo(1:n2);
            xo2 = flip(xo(n2+1:end));
            yo2 = flip(yo(n2+1:end));

            d1 = median(min(pdist2([xo1,yo1],[xa,ya],'euclidean'),[],2));
            d2 = median(min(pdist2([xo2,yo2],[xa,ya],'euclidean'),[],2));

            if i2 == 1
                varargout = {xo1,yo1,xo2,yo2};
            else
                if d1 > d2
                    varargout = {xo1,yo1};
                else
                    varargout = {xo2,yo2};
                end
            end
        end

        function obj = unpad(obj,DEM)

            obj.DEM = DEM;
            obj.pair = [];

            [X,Y] = worldGrid(DEM.georef);
            xmin = min(X(:)); xmax = max(X(:));
            ymin = min(Y(:)); ymax = max(Y(:));

            if isfield(obj.base,'x') && isfield(obj.base,'y')
                x = obj.base.x(:);
                y = obj.base.y(:);
                src = obj.base.src(:);

                ok = x >= xmin & x <= xmax & y >= ymin & y <= ymax;
                x = x(ok);
                y = y(ok);
                src = src(ok);

                ix = coord2ind(DEM,x,y);
                ok = ~isnan(ix);

                x = x(ok);
                y = y(ok);
                ix = ix(ok);
                src = src(ok);

                obj.base.x = x(:);
                obj.base.y = y(:);
                obj.base.ix = ix(:);
                obj.base.dist = [0; cumsum(hypot(diff(x(:)),diff(y(:))))];
                obj.base.src = src(:);
            end

            cs = DEM.cellsize/2;
            Z = DEM.Z;

            for i1 = 1:2
                for i2 = 1:numel(obj.conn{i1})
                    for i3 = 1:numel(obj.conn{i1}(i2).x)

                        x0 = obj.conn{i1}(i2).x{i3};
                        y0 = obj.conn{i1}(i2).y{i3};

                        if numel(x0) < 2
                            continue
                        end

                        d0 = obj.conn{i1}(i2).d{i3};

                        I0 = x0 >= xmin & x0 <= xmax & y0 >= ymin & y0 <= ymax;
                        ci = coord2ind(DEM,x0(I0),y0(I0));
                        ci = unique(ci(~isnan(ci)),'stable');

                        intD = 0:cs:d0(end);
                        xi = interp1(d0,x0,intD,'linear');
                        yi = interp1(d0,y0,intD,'linear');

                        I = xi >= xmin & xi <= xmax & yi >= ymin & yi <= ymax;
                        xi = xi(I);
                        yi = yi(I);

                        ni = coord2ind(DEM,xi,yi);
                        ni = unique(ni(~isnan(ni)),'stable');

                        if isempty(ni)
                            obj.conn{i1}(i2).ix{i3} = [];
                            obj.conn{i1}(i2).x{i3} = [];
                            obj.conn{i1}(i2).y{i3} = [];
                            obj.conn{i1}(i2).z{i3} = [];
                            obj.conn{i1}(i2).d{i3} = [];

                            obj.int{i1}(i2).ix{i3} = [];
                            obj.int{i1}(i2).x{i3} = [];
                            obj.int{i1}(i2).y{i3} = [];
                            obj.int{i1}(i2).z{i3} = [];
                            obj.int{i1}(i2).d{i3} = [];
                            continue
                        end

                        if any(~I) && ~isempty(ci)
                            ci = unique([ci(:); ni(end)],'stable');
                        end

                        [xc,yc] = ind2coord(DEM,ci);
                        [xi,yi] = ind2coord(DEM,ni);

                        obj.conn{i1}(i2).ix{i3} = ci;
                        obj.conn{i1}(i2).x{i3} = xc;
                        obj.conn{i1}(i2).y{i3} = yc;
                        obj.conn{i1}(i2).z{i3} = Z(ci);
                        obj.conn{i1}(i2).d{i3} = [0; cumsum(hypot(diff(xc),diff(yc)))];

                        obj.int{i1}(i2).ix{i3} = ni;
                        obj.int{i1}(i2).x{i3} = xi;
                        obj.int{i1}(i2).y{i3} = yi;
                        obj.int{i1}(i2).z{i3} = Z(ni);
                        obj.int{i1}(i2).d{i3} = [0; cumsum(hypot(diff(xi),diff(yi)))];
                    end
                end
            end
        end
    end
end
