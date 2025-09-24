function T = pairing(T, varargin)
    p = inputParser;
    addRequired(p, 'T');
    isValidFP = @(v) (islogical(v) || (isnumeric(v) && ismatrix(v))) && isequal(size(v), size(T.DEM.Z));
    isValidX = @isnumeric;
    isValidY = @(v) isnumeric(v) && isequal(size(v), size(p.Results.x));

    method_default = 'intersect';
    if any(strcmpi(varargin, 'method'))
        idx = find(strcmpi(varargin, 'method'),1,'first');
        method_default = varargin{idx+1};
    end

    if nargin == 2 || (nargin > 2 && ischar(varargin{2}))
        addRequired(p, 'FP', isValidFP);
        addParameter(p, 'type', 'any', @(v) ismember(v, {'both', 'any'}));
        addParameter(p, 'connectivity', 'D8', @(v) ismember(v, {'D1', 'D8', 'D16'}));
        addParameter(p, 'verbose', false, @islogical);
        addParameter(p, 'parallel', false, @islogical);
        addParameter(p, 'maxlength', false, @islogical);
        addParameter(p, 'gap', 1, @(v) isnumeric(v) && isscalar(v) && v >= 0);
        parse(p, T, varargin{:});
        FP = p.Results.FP;
        met = 'mask';
        parFlag = p.Results.parallel;
        tp = p.Results.type;
        con = p.Results.connectivity;
        vb = p.Results.verbose;
        maxlength = p.Results.maxlength;
        gap = p.Results.gap;
        useFP = true;
        ord = [];
    else
        addRequired(p, 'x', isValidX);
        addRequired(p, 'y', isValidY);
        addParameter(p, 'method', 'intersect', @(v) ismember(v, {'intersect', 'mask'}));
        addParameter(p, 'type', 'any', @(v) ismember(v, {'both', 'any'}));
        addParameter(p, 'connectivity', 'D8', @(v) ismember(v, {'D1', 'D8', 'D16'}));
        addParameter(p, 'verbose', false, @islogical);
        addParameter(p, 'parallel', false, @islogical);
        addParameter(p, 'maxlength', false, @islogical);
        addParameter(p, 'gap', 1, @(v) isnumeric(v) && isscalar(v) && v >= 0);
        if strcmp(method_default, 'intersect')
            addParameter(p, 'order', 'first', @(v) ismember(v, {'first','last'}));
        end
        parse(p, T, varargin{:});
        x = p.Results.x;
        y = p.Results.y;
        if isfield(p.Results, 'method')
            met = p.Results.method;
        else
            met = 'intersect';
        end
        parFlag = p.Results.parallel;
        tp = p.Results.type;
        con = p.Results.connectivity;
        vb = p.Results.verbose;
        maxlength = p.Results.maxlength;
        gap = p.Results.gap;
        if strcmp(met, 'intersect')
            ord = p.Results.order;
        else
            ord = [];
        end
        useFP = false;
        FP = [];
    end

    if parFlag
        T.DEM.Z = gpuArray(T.DEM.Z);
        T.DEM.refmat = gpuArray(T.DEM.refmat);
        if useFP
            FP = gpuArray(FP);
        end
        if exist('x', 'var')
            x = gpuArray(x); y = gpuArray(y);
        end
    end

    nC = numel(T.x);
    if vb, PB = ProgressBar(nC, 'taskname', 'Pairing...', 'ui', 'cli'); end

    if strcmp(met, 'intersect')
        [X,Y] = refmat2XY(T.DEM.refmat, T.DEM.size); X = X(:); Y = Y(:);
        x = x(:); y = y(:);
        dx = X(2)-X(1); dy = Y(2)-Y(1);
        ix1 = round((x-X(1))/dx+1);
        ix2 = round((y-Y(1))/dy+1);
        ob = ix1 > T.DEM.size(2) | ix1 < 1 | ix2 > T.DEM.size(1) | ix2 < 1;
        if any(ob)
            warning('TopoToolbox:outsidegrid','Some pts out-of-bound.');
            x(ob)=[]; y(ob)=[]; 
        end
        if isempty(x)
            warning('No valid pts.');
            return;
        end
        ind = coord2ind(T.DEM, x, y);
    end

    if parFlag
        Ttemp = cell(nC,1);
        parfor i2 = 1:nC
            if useFP
                Ttemp{i2} = pair(T, i2, ind, con, tp, met, FP, 'gap', gap);
            else
                Ttemp{i2} = pair(T, i2, ind, con, tp, met, 'order', ord, 'gap', gap);
            end
            if vb, PB.count(); end
        end
        for i2 = 1:nC
            for i1 = 1:2
                T.int{i1}(i2).x = Ttemp{i2}.int{i1}(i2).x;
                T.int{i1}(i2).y = Ttemp{i2}.int{i1}(i2).y;
                T.int{i1}(i2).ix = Ttemp{i2}.int{i1}(i2).ix;
                T.int{i1}(i2).z = Ttemp{i2}.int{i1}(i2).z;
                T.int{i1}(i2).d = Ttemp{i2}.int{i1}(i2).d;

                T.conn{i1}(i2).x = Ttemp{i2}.conn{i1}(i2).x;
                T.conn{i1}(i2).y = Ttemp{i2}.conn{i1}(i2).y;
                T.conn{i1}(i2).ix = Ttemp{i2}.conn{i1}(i2).ix;
                T.conn{i1}(i2).z = Ttemp{i2}.conn{i1}(i2).z;
                T.conn{i1}(i2).d = Ttemp{i2}.conn{i1}(i2).d;
            end
        end
    else
        for i2 = 1:nC
            if useFP
                ind = [];
                T = pair(T, i2, ind, con, tp, met, 'FP', FP, 'gap', gap);
            else
                T = pair(T, i2, ind, con, tp, met, 'order', ord, 'gap', gap);
            end
            if vb, PB.count(); end
        end
    end

    if strcmp(met, 'mask')
        bw = bwdist(~FP)*T.DEM.cellsize*2;
        T = mask_validation(T, con, bw, 'maxlength', maxlength);
    end

    if parFlag
        T.DEM.Z = gather(T.DEM.Z);
        T.DEM.refmat = gather(T.DEM.refmat);
    end

    %% Compute statistics
            
    xii = cell(2,1);
    yii = cell(2,1);
    zii = cell(2,1);
    dii = cell(2,1);
    
    for i1 = 1:2
        n2 = numel(T.x);
        xii{i1} = cell(n2, 1);
        yii{i1} = cell(n2, 1);
        zii{i1} = cell(n2, 1);
        dii{i1} = cell(n2, 1);
    
        for i2 = 1:n2
            if ~isempty(T.int{i1}(i2).x)
                midx = ceil(numel(T.int{i1}(i2).x)/2);
                xii{i1}{i2} = T.int{i1}(i2).x{midx}(:);
                yii{i1}{i2} = T.int{i1}(i2).y{midx}(:);
                zii{i1}{i2} = T.int{i1}(i2).z{midx}(:);
                dii{i1}{i2} = T.int{i1}(i2).d{midx}(:);
            end
        end
    end

    T.stats.x = cell(numel(T.x),1);
    T.stats.y = cell(numel(T.x),1);
    T.stats.z = cell(numel(T.x),1);
    T.stats.d = cell(numel(T.x),1);

    T.stats.Z1 = nan(numel(T.x),1);
    T.stats.Z2 = nan(numel(T.x),1);
    T.stats.G1 = nan(numel(T.x),1);
    T.stats.G2 = nan(numel(T.x),1);
    T.stats.W = nan(numel(T.x),1);
    T.stats.H = nan(numel(T.x),1);

    
    for i2 = 1:numel(T.x)
    
        try                    
            T.stats.x{i2} = [flip(xii{1}{i2}); xii{2}{i2}(2:end)];
            T.stats.y{i2} = [flip(yii{1}{i2}); yii{2}{i2}(2:end)];
            T.stats.z{i2} = [flip(zii{1}{i2}); zii{2}{i2}(2:end)];
            T.stats.d{i2} = [flip(dii{1}{i2} * (-1)); dii{2}{i2}(2:end)];
        catch
            T.stats.x{i2} = [];
            T.stats.y{i2} = [];
            T.stats.z{i2} = [];
            T.stats.d{i2} = [];
        end
    end
    
    for i1 = 1:numel(T.x)
    
        try
            id1 = find(T.stats.d{i1} < 0);
            id2 = find(T.stats.d{i1} > 0);
            Z1 = T.stats.z{i1}(id1);
            Z2 = T.stats.z{i1}(id2);
            D1 = T.stats.d{i1}(id1);
            D2 = T.stats.d{i1}(id2);
        catch
        end
    
        try
            T.stats.Z1(i1) = max(Z1);
            G1 = abs(atan(diff(Z1) ./ diff(abs(D1))) * (180/pi)); 
            T.stats.G1(i1) = mean(G1);
        catch
        end
             
        try
            T.stats.Z2(i1) = max(Z2);
            G2 = abs(atan(diff(Z2) ./ diff(abs(D2))) * (180/pi)); 
            T.stats.G2(i1) = mean(G2);
        catch
        end
    
        try
            T.stats.W(i1) = T.stats.d{i1}(1) * (-1) + T.stats.d{i1}(end);
            T.stats.H(i1) = mean([T.stats.Z1(i1) T.stats.Z2(i1)]) - min(T.stats.z{i1});
        catch
        end
    end
end

function T = pair(T, i2, ind, con, tp, met, varargin)

    p = inputParser;
    addParameter(p, 'gap', 1, @(v) isnumeric(v) && isscalar(v) && v >= 0);
    addParameter(p, 'order', 'first', @(v) ismember(v, {'first','last'}));   
    addParameter(p, 'FP', [], @(v) islogical(v) || isnumeric(v));
    parse(p, varargin{:});
    gap = p.Results.gap;
    ord = p.Results.order;
    FP = p.Results.FP;

    switch con
        case 'D1'
            dr = 0; dc = 0;
        case 'D8'
            dr = [-1 -1 -1 0 0 0 1 1 1];
            dc = [-1 0 1 -1 0 1 -1 0 1];
        case 'D16'
            dr = [-2 -2 -2 -2 -1 -1 -1 -1 0 0 0 0 0 1 1 1 1 2 2 2 2];
            dc = [-2 -1 0 1 -2 -1 0 1 -2 -1 0 1 2 -2 -1 0 1 -2 -1 0 1];
        otherwise
            error('Unknown connectivity: %s', con);
    end

    [NR, NC] = size(T.DEM.Z);

    for i1 = 1:2
        for i3 = 1:numel(T.int{i1}(i2).x)
            if numel(T.int{i1}(i2).x{i3}) == 1
                continue;
            end

            xi = T.int{i1}(i2).x{i3};
            yi = T.int{i1}(i2).y{i3};
            xi = xi(:); yi = yi(:);
            valid = ~(isnan(xi) | isnan(yi));
            xi = xi(valid); yi = yi(valid);
            if numel(xi) < 2, continue; end
            d = [0; cumsum(sqrt(diff(xi).^2 + diff(yi).^2))];
            di = 0:T.DEM.cellsize/5:d(end);
            xii = interp1(d, xi, di, 'linear');
            yii = interp1(d, yi, di, 'linear');
            ixii = coord2ind(T.DEM, xii, yii);

            if strcmp(met, 'mask')
                msk = FP(ixii);

                n = numel(msk);
                last_end = 0; curr_gap = 0;
                for i = 1:n
                    if msk(i)
                        curr_gap = 0;
                    else
                        if curr_gap < gap
                            curr_gap = curr_gap + 1;
                        else
                            break
                        end
                    end
                    last_end = i;
                end
                if last_end > 0
                    p = last_end;
                else
                    p = [];
                end
                nd = 0; IXn = [];
            else
                [rA, cA] = ind2sub([NR, NC], ixii);
                R = rA(2:end);
                C = cA(2:end);
                nC = numel(R);
                Rn = repmat(R, 1, numel(dr)) + repmat(dr, nC, 1);
                Cn = repmat(C, 1, numel(dc)) + repmat(dc, nC, 1);
                vM = (Rn >= 1) & (Rn <= NR) & (Cn >= 1) & (Cn <= NC);
                nbA = nan(size(Rn));
                nbA(vM) = sub2ind([NR, NC], Rn(vM), Cn(vM));

                dN0 = sqrt(dr.^2 + dc.^2) * T.DEM.cellsize;
                dNA = repmat(dN0, nC, 1);
                mskA = ismember(nbA, ind);
                cf = any(mskA, 2);
                dNA(~mskA) = Inf;
                mDC = min(dNA, [], 2);
                cK = find(cf) + 1;
                cD = mDC(cf);

                if isempty(cK)
                    p = [];
                else
                    if strcmp(ord, 'first')
                        cl = cumsum([1; diff(cK) > gap]);
                        fc = cK(cl == cl(1));
                        cD2 = cD(cl == cl(1));
                    else
                        cl = cumsum([1; diff(cK) > gap]);
                        fc = cK(cl == cl(end));
                        cD2 = cD(cl == cl(end));
                    end

                    [~, im] = min(cD2);
                    p = fc(im);

                    ri = p - 1;
                    dN_row = dNA(ri, :);
                    mR = mskA(ri, :);
                    dN_row(~mR) = Inf;
                    [min_val, ci] = min(dN_row);
                    nd = min_val;
                    IXn = nbA(ri, ci);
                end
            end

            if ~isempty(p)
                xS = xii(1:p)'; yS = yii(1:p)';
                dS = [0; cumsum(sqrt(diff(xS).^2 + diff(yS).^2))];
                if nd == 0 || strcmp(con, 'D1')
                    T.int{i1}(i2).x{i3} = xS(:);
                    T.int{i1}(i2).y{i3} = yS(:);
                    T.int{i1}(i2).ix{i3} = ixii(1:p);
                    T.int{i1}(i2).z{i3} = T.DEM.Z(ixii(1:p));
                    T.int{i1}(i2).d{i3} = dS(:);
                else
                    [xn, yn] = ind2coord(T.DEM, IXn);
                    zn = T.DEM.Z(IXn);
                    dnS = dS(end) + T.DEM.cellsize;
                    T.int{i1}(i2).x{i3} = [xS; xn];
                    T.int{i1}(i2).y{i3} = [yS; yn];
                    T.int{i1}(i2).ix{i3} = [ixii(1:p); IXn];
                    T.int{i1}(i2).z{i3} = [T.DEM.Z(ixii(1:p)); zn];
                    T.int{i1}(i2).d{i3} = [dS; dnS];
                end

                existing = T.conn{i1}(i2).ix{i3};
                comp = ixii(1:p);
                if ~isempty(existing)
                    c0 = existing(ismember(existing, comp));
                else
                    c0 = [];
                end

                if ~isempty(IXn) && ~ismember(IXn, c0)
                    c0 = [c0; IXn];
                end


                [cx, cy] = ind2coord(T.DEM, c0);
                T.conn{i1}(i2).x{i3} = cx;
                T.conn{i1}(i2).y{i3} = cy;
                T.conn{i1}(i2).ix{i3} = c0;
                T.conn{i1}(i2).z{i3} = T.DEM.Z(c0);
                T.conn{i1}(i2).d{i3} = [0; cumsum(sqrt(diff(cx).^2 + diff(cy).^2))];
            else
                T.int{i1}(i2).x{i3} = [];
                T.int{i1}(i2).y{i3} = [];
                T.int{i1}(i2).ix{i3} = [];
                T.int{i1}(i2).z{i3} = [];
                T.int{i1}(i2).d{i3} = [];
                T.conn{i1}(i2).x{i3} = [];
                T.conn{i1}(i2).y{i3} = [];
                T.conn{i1}(i2).ix{i3} = [];
                T.conn{i1}(i2).z{i3} = [];
                T.conn{i1}(i2).d{i3} = [];
            end
        end
    end

    v1 = any(cellfun(@(x) ~isempty(x), T.int{1}(i2).x));
    v2 = any(cellfun(@(x) ~isempty(x), T.int{2}(i2).x));
    if strcmp(tp, 'both')
        vp = v1 && v2;
    else
        vp = v1 || v2;
    end
    if ~vp
        for i1 = 1:2
            for i3 = 1:numel(T.conn{i1}(i2).x)
                T.int{i1}(i2).x{i3} = [];
                T.int{i1}(i2).y{i3} = [];
                T.int{i1}(i2).ix{i3} = [];
                T.int{i1}(i2).z{i3} = [];
                T.int{i1}(i2).d{i3} = [];
                T.conn{i1}(i2).x{i3} = [];
                T.conn{i1}(i2).y{i3} = [];
                T.conn{i1}(i2).ix{i3} = [];
                T.conn{i1}(i2).z{i3} = [];
                T.conn{i1}(i2).d{i3} = [];
            end
        end
    end
end

function T = mask_validation(T, con, bw, varargin)
    p = inputParser;
    addParameter(p, 'maxlength', false, @islogical);
    parse(p, varargin{:});
    mlen = p.Results.maxlength;

    switch con
        case 'D1'
            dr = 0; dc = 0;
        case 'D8'
            dr = [-1 -1 -1 0 0 0 1 1 1];
            dc = [-1 0 1 -1 0 1 -1 0 1];
        case 'D16'
            dr = [-2 -2 -2 -2 -1 -1 -1 -1 0 0 0 0 0 1 1 1 1 2 2 2 2];
            dc = [-2 -1 0 1 -2 -1 0 1 -2 -1 0 1 2 -2 -1 0 1 -2 -1 0 1];
        otherwise
            error('Unknown connectivity: %s', con);
    end

    nB = numel(T.x);
    for i2 = 1:nB
        id = [T.int{1}(i2).ix{1}; T.int{2}(i2).ix{1}];
        if isempty(id), continue; end
        [r,c] = ind2sub(size(bw), id);
        mb = zeros(size(id));
        for k = 1:numel(id)
            r0 = r(k); c0 = c(k);
            nr0 = r0 + dr; nc0 = c0 + dc;
            v = (nr0>=1)&(nr0<=size(bw,1))&(nc0>=1)&(nc0<=size(bw,2));
            nId = sub2ind(size(bw), nr0(v), nc0(v));
            mb(k) = max(bw(nId));
        end
        [mbv, ia] = max(mb);
        mbwI = nan(2,1);
        [~, ib, ~] = intersect(T.int{1}(i2).ix{1}, id(ia));
        if ~isempty(ib)
            mbwI(1) = T.int{1}(i2).d{1}(ib) + mbv;
            mbwI(2) = mbv - T.int{1}(i2).d{1}(ib);
        else
            [~, ib, ~] = intersect(T.int{2}(i2).ix{1}, id(ia));
            mbwI(2) = T.int{2}(i2).d{1}(ib) + mbv;
            mbwI(1) = mbv - T.int{2}(i2).d{1}(ib);
        end
        lens = zeros(2,1);
        for i1 = 1:2
            for i3 = 1:numel(T.int{i1}(i2).x)
                if isempty(T.int{i1}(i2).x{i3}), continue; end
                xii = T.int{i1}(i2).x{i3};
                yii = T.int{i1}(i2).y{i3};
                di = T.int{i1}(i2).d{i3};
                if isempty(di)||isnan(mbwI(i1))
                    T.int{i1}(i2).x{i3} = [];
                    T.int{i1}(i2).y{i3} = [];
                    T.int{i1}(i2).ix{i3} = [];
                    T.int{i1}(i2).z{i3} = [];
                    T.int{i1}(i2).d{i3} = [];
                    T.conn{i1}(i2).x{i3} = [];
                    T.conn{i1}(i2).y{i3} = [];
                    T.conn{i1}(i2).ix{i3} = [];
                    T.conn{i1}(i2).z{i3} = [];
                    T.conn{i1}(i2).d{i3} = [];
                    continue;
                end
                [~, p] = min(abs(mbwI(i1)-di));
                if p < 1 || p > numel(xii)
                    T.int{i1}(i2).x{i3} = [];
                    T.int{i1}(i2).y{i3} = [];
                    T.int{i1}(i2).ix{i3} = [];
                    T.int{i1}(i2).z{i3} = [];
                    T.int{i1}(i2).d{i3} = [];
                    T.conn{i1}(i2).x{i3} = [];
                    T.conn{i1}(i2).y{i3} = [];
                    T.conn{i1}(i2).ix{i3} = [];
                    T.conn{i1}(i2).z{i3} = [];
                    T.conn{i1}(i2).d{i3} = [];
                    continue;
                end
                lens(i1) = p;
                T.int{i1}(i2).x{i3} = xii(1:p);
                T.int{i1}(i2).y{i3} = yii(1:p);
                T.int{i1}(i2).ix{i3} = coord2ind(T.DEM, xii(1:p), yii(1:p));
                T.int{i1}(i2).z{i3} = T.DEM.Z(T.int{i1}(i2).ix{i3});
                T.int{i1}(i2).d{i3} = [0; cumsum(sqrt(diff(xii(1:p)).^2+diff(yii(1:p)).^2))];
            end
        end

        if mlen
            nlim = min(lens(lens > 0));
            for i1 = 1:2
                for i3 = 1:numel(T.int{i1}(i2).x)
                    if isempty(T.int{i1}(i2).x{i3}), continue; end
                    p = min(nlim, numel(T.int{i1}(i2).x{i3}));
                    T.int{i1}(i2).x{i3} = T.int{i1}(i2).x{i3}(1:p);
                    T.int{i1}(i2).y{i3} = T.int{i1}(i2).y{i3}(1:p);
                    T.int{i1}(i2).ix{i3} = T.int{i1}(i2).ix{i3}(1:p);
                    T.int{i1}(i2).z{i3} = T.int{i1}(i2).z{i3}(1:p);
                    T.int{i1}(i2).d{i3} = T.int{i1}(i2).d{i3}(1:p);
                    xi = [T.int{i1}(i2).x{i3}(1); T.int{i1}(i2).x{i3}(end)];
                    yi = [T.int{i1}(i2).y{i3}(1); T.int{i1}(i2).y{i3}(end)];
                    ix = unique(coord2ind(T.DEM, xi, yi), 'stable');
                    [xi, yi] = ind2coord(T.DEM, ix);
                    T.conn{i1}(i2).x{i3} = xi;
                    T.conn{i1}(i2).y{i3} = yi;
                    T.conn{i1}(i2).ix{i3} = ix;
                    T.conn{i1}(i2).z{i3} = T.DEM.Z(ix);
                    T.conn{i1}(i2).d{i3} = [0; cumsum(sqrt(diff(xi).^2+diff(yi).^2))];
                end
            end
        else
            for i1 = 1:2
                for i3 = 1:numel(T.int{i1}(i2).x)
                    if isempty(T.int{i1}(i2).x{i3}), continue; end
                    xi = [T.int{i1}(i2).x{i3}(1); T.int{i1}(i2).x{i3}(end)];
                    yi = [T.int{i1}(i2).y{i3}(1); T.int{i1}(i2).y{i3}(end)];
                    ix = unique(coord2ind(T.DEM, xi, yi), 'stable');
                    [xi, yi] = ind2coord(T.DEM, ix);
                    T.conn{i1}(i2).x{i3} = xi;
                    T.conn{i1}(i2).y{i3} = yi;
                    T.conn{i1}(i2).ix{i3} = ix;
                    T.conn{i1}(i2).z{i3} = T.DEM.Z(ix);
                    T.conn{i1}(i2).d{i3} = [0; cumsum(sqrt(diff(xi).^2+diff(yi).^2))];
                end
            end
        end
    end
end
