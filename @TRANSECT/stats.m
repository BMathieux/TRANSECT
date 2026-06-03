function T = stats(T,varargin)
%STATS Compute transect statistics from a TRANSECT object (single-path only)
%
% Syntax:
%   T = stats(T)
%   T = stats(T,'func',@median)
%
% Options:
%   'func' - function handle used for slope aggregation
%            Default = @mean
%
% Always computes: Z1, Z2, G1, G2, W, H
% Statistics are indexed along the transect nodes stored in T.int/T.conn.

p = inputParser;
p.addParameter('func',@mean,@(f) isa(f,'function_handle'));
p.parse(varargin{:});
aggfun = p.Results.func;

if isempty(T.int)
    n = 0;
else
    n = numel(T.int{1});
end

T.stats.x = cell(n,1);
T.stats.y = cell(n,1);
T.stats.z = cell(n,1);
T.stats.d = cell(n,1);
T.stats.Z1 = nan(n,1);
T.stats.Z2 = nan(n,1);
T.stats.G1 = nan(n,1);
T.stats.G2 = nan(n,1);
T.stats.W  = nan(n,1);
T.stats.H  = nan(n,1);

for i = 1:n
    try
        midL = ceil(numel(T.int{1}(i).x)/2);
        xL = flip(T.int{1}(i).x{midL});
        yL = flip(T.int{1}(i).y{midL});
        zL = flip(T.int{1}(i).z{midL});
        dL = flip(T.int{1}(i).d{midL} * -1);

        midR = ceil(numel(T.int{2}(i).x)/2);
        xR = T.int{2}(i).x{midR};
        yR = T.int{2}(i).y{midR};
        zR = T.int{2}(i).z{midR};
        dR = T.int{2}(i).d{midR};

        xsec = [xL; xR(2:end)];
        ysec = [yL; yR(2:end)];
        zsec = [zL; zR(2:end)];
        dsec = [dL; dR(2:end)];
    catch
        continue
    end

    T.stats.x{i} = xsec;
    T.stats.y{i} = ysec;
    T.stats.z{i} = zsec;
    T.stats.d{i} = dsec;

    id1 = dsec<0; id2 = dsec>0;
    if any(id1), T.stats.Z1(i) = max(zsec(id1)); end
    if any(id2), T.stats.Z2(i) = max(zsec(id2)); end
    if any(id1(:)) && any(id2(:))
        T.stats.W(i) = abs(min(dsec(id1))) + max(dsec(id2));
        T.stats.H(i) = mean([T.stats.Z1(i) T.stats.Z2(i)],'omitnan') - min(zsec);
    end

    if sum(id1)>1
        g1=atand(diff(zsec(id1))./diff(abs(dsec(id1))));
        v=aggfun(abs(g1)); if numel(v)>1, v=nanmean(v(:)); end
        T.stats.G1(i)=v;
    end
    if sum(id2)>1
        g2=atand(diff(zsec(id2))./diff(abs(dsec(id2))));
        v=aggfun(abs(g2)); if numel(v)>1, v=nanmean(v(:)); end
        T.stats.G2(i)=v;
    end
end
end
