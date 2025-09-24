function T = resampleT(T, DEM)
%RESAMPLET Resample transect connectivity and interpolated data
%   T = RESAMPLET(T, DEM) resamples the transect struct T.conn and T.int
%   fields using linear interpolation at intervals of DEM.cellsize/2.
%   Updates T.conn.ix, x, y and T.int.ix, x, y, z, d fields based on DEM.
%
%   Inputs:
%       T     (struct) - Transect struct with fields:
%                        - conn: Connectivity data (ix, x, y, d)
%                        - int: Interpolated data (ix, x, y, z, d)
%       DEM   (struct) - Digital Elevation Model with fields:
%                        - Z: Elevation data
%                        - refmat: Reference matrix
%                        - size: Size of DEM
%                        - cellsize: Cell size for resampling
%
%   Outputs:
%       T     (struct) - Updated transect struct with resampled fields

T.DEM = DEM;
for i1 = 1:2
    for i2 = 1:numel(T.x)
        for i3 = 1:numel(T.conn{i1}(i2).x)
            % Extract path data
            xii = T.conn{i1}(i2).x{i3};
            yii = T.conn{i1}(i2).y{i3};
            dii = cumsum(T.conn{i1}(i2).d{i3});

            % Check if path has at least two points for interpolation
            if numel(xii) < 2 || numel(dii) < 2
                % Skip empty or single-point paths
                T.conn{i1}(i2).ix{i3} = [];
                T.conn{i1}(i2).x{i3} = [];
                T.conn{i1}(i2).y{i3} = [];
                T.int{i1}(i2).ix{i3} = [];
                T.int{i1}(i2).x{i3} = [];
                T.int{i1}(i2).y{i3} = [];
                T.int{i1}(i2).z{i3} = [];
                T.int{i1}(i2).d{i3} = [];
                continue;
            end

            % Interpolate at half cell size intervals
            intD = 0:DEM.cellsize/2:dii(end);
            nx = interp1(dii, xii, intD, 'linear');
            ny = interp1(dii, yii, intD, 'linear');

            % Get coordinate vectors
            [X, Y] = refmat2XY(DEM.refmat, DEM.size);
            X = X(:);
            Y = Y(:);

            % Remove points outside DEM bounds
            I = nx > max(X) | nx < min(X) | ny > max(Y) | ny < min(Y);
            nx(I) = [];
            ny(I) = [];

            I = xii > max(X) | xii < min(X) | yii > max(Y) | yii < min(Y);
            xii(I) = [];
            yii(I) = [];

            % Compute distances and indices
            nd = [0 sqrt(diff(nx).^2 + diff(ny).^2)];
            nid = coord2ind(DEM, nx, ny);
            nid(isnan(nid)) = [];
            [nid, idx] = unique(nid, 'stable');

            % Convert indices back to coordinates
            [nx, ny] = ind2coord(DEM, nid);
            nd = nd(idx);
            nz = DEM.Z(nid);

            % Update T.conn fields
            T.conn{i1}(i2).ix{i3} = coord2ind(DEM, xii, yii);
            if any(I)
                T.conn{i1}(i2).ix{i3} = unique([T.conn{i1}(i2).ix{i3}; nid(end)], 'stable');
            end
            [T.conn{i1}(i2).x{i3}, T.conn{i1}(i2).y{i3}] = ...
                ind2coord(DEM, T.conn{i1}(i2).ix{i3});

            % Update T.int fields
            T.int{i1}(i2).ix{i3} = nid;
            T.int{i1}(i2).x{i3} = nx;
            T.int{i1}(i2).y{i3} = ny;
            T.int{i1}(i2).z{i3} = nz;
            T.int{i1}(i2).d{i3} = cumsum(nd);
        end
    end
end
end