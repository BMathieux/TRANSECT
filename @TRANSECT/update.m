function T = update(T,DEM2)
%UPDATE Update transect elevations on a new DEM
%
%   T = UPDATE(T,DEM2) updates the z-values of transect paths (T.int and
%   T.conn) using elevations from another DEM (GRIDobj). If DEM2 is not
%   aligned with T.DEM, UPDATE calls T.resample(DEM2) instead.
%
%   Inputs:
%       T     TRANSECT object
%       DEM2  GRIDobj aligned with T.DEM
%
%   Output:
%       T     Updated TRANSECT object
%
%   See also: TRANSECT.resample, TRANSECT.calcStats, GRIDobj/validatealignment

if ~isa(T,'TRANSECT') || ~isa(DEM2,'GRIDobj')
    error('Inputs must be a TRANSECT object and a GRIDobj.');
end

tf = validatealignment(T.DEM,DEM2);
if ~tf
    warning('DEM grids not aligned. Resampling transect using new DEM...');
    T = T.resample(DEM2);
    return
end

for i1 = 1:2
    for i2 = 1:numel(T.int{i1})
        for i3 = 1:numel(T.int{i1}(i2).ix)
            ix = T.int{i1}(i2).ix{i3};
            if isempty(ix), continue, end
            ix(ix>numel(DEM2.Z)) = [];
            T.int{i1}(i2).z{i3} = DEM2.Z(ix);
        end
    end
end

T.DEM = DEM2;
T = T.calcStats();
end
