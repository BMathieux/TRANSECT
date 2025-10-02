function T = resample(T, DEM)
%RESAMPLE Resample transect connectivity and interpolation data
%
%   T = RESAMPLE(T, DEM) resamples transect connectivity (T.conn) and
%   interpolated paths (T.int) at fixed intervals of DEM.cellsize/2 using
%   linear interpolation. Updates node indices, coordinates, elevations,
%   and cumulative distances based on DEM.
%
%   Inputs:
%       T   - TRANSECT object with fields:
%             * conn: connectivity data (ix, x, y, d)
%             * int : interpolated data (ix, x, y, z, d)
%       DEM - GRIDobj (TopoToolbox) representing the DEM
%
%   Output:
%       T   - Updated TRANSECT object with resampled fields

    %--- input validation
    if ~isa(T,'TRANSECT')
        error('First input must be a TRANSECT object.')
    end
    if ~isa(DEM,'GRIDobj')
        error('DEM must be a GRIDobj.')
    end

    % update DEM reference
    T.DEM = DEM;

    % iterate over both sides of transect
    for i1 = 1:2
        for i2 = 1:numel(T.x)
            for i3 = 1:numel(T.conn{i1}(i2).x)

                % extract path
                xii = T.conn{i1}(i2).x{i3};
                yii = T.conn{i1}(i2).y{i3};
                dii = cumsum(T.conn{i1}(i2).d{i3});

                % skip if empty
                if numel(xii) < 2 || numel(dii) < 2
                    T.conn{i1}(i2).ix{i3} = [];
                    T.conn{i1}(i2).x{i3}  = [];
                    T.conn{i1}(i2).y{i3}  = [];
                    T.int{i1}(i2).ix{i3}  = [];
                    T.int{i1}(i2).x{i3}   = [];
                    T.int{i1}(i2).y{i3}   = [];
                    T.int{i1}(i2).z{i3}   = [];
                    T.int{i1}(i2).d{i3}   = [];
                    continue
                end

                % interpolate at half-cell spacing
                intD = 0:DEM.cellsize/2:dii(end);
                nx = interp1(dii,xii,intD,'linear');
                ny = interp1(dii,yii,intD,'linear');

                % get DEM bounds
                [X,Y] = refmat2XY(DEM.refmat,DEM.size);
                X = X(:); Y = Y(:);

                % remove out-of-bounds
                I = nx > max(X) | nx < min(X) | ny > max(Y) | ny < min(Y);
                nx(I) = []; ny(I) = [];
                I = xii > max(X) | xii < min(X) | yii > max(Y) | yii < min(Y);
                xii(I) = []; yii(I) = [];

                % compute distances + indices
                nd = [0 sqrt(diff(nx).^2+diff(ny).^2)];
                nid = coord2ind(DEM,nx,ny);
                nid(isnan(nid)) = [];
                [nid,idx] = unique(nid,'stable');

                % convert back
                [nx,ny] = ind2coord(DEM,nid);
                nd = nd(idx);
                nz = DEM.Z(nid);

                % update connectivity
                T.conn{i1}(i2).ix{i3} = coord2ind(DEM,xii,yii);
                if any(I)
                    T.conn{i1}(i2).ix{i3} = unique([T.conn{i1}(i2).ix{i3}; nid(end)],'stable');
                end
                [T.conn{i1}(i2).x{i3}, T.conn{i1}(i2).y{i3}] = ...
                    ind2coord(DEM,T.conn{i1}(i2).ix{i3});

                % update interpolated fields
                T.int{i1}(i2).ix{i3} = nid;
                T.int{i1}(i2).x{i3}  = nx;
                T.int{i1}(i2).y{i3}  = ny;
                T.int{i1}(i2).z{i3}  = nz;
                T.int{i1}(i2).d{i3}  = cumsum(nd);
            end
        end
    end

    % compute statistics

    T.calcStats();
end
