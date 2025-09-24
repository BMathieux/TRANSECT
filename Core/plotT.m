function plotT(T, varargin)
% PLOTT Plots paired connections and end nodes from the structure T.
%
%   PLOTT(T) plots the paired connections stored in structure T,
%   including end nodes and base layer.
%
%   PLOTT(T, 'MarkerSize', VALUE, 'LineWidth', VALUE) allows customization:
%       'MarkerSize' - Size of scatter markers (default: 15)
%       'LineWidth'  - Width of connection lines (default: 1)
%
%   Inputs:
%       T  - Structure containing spatial data with fields:
%            - T.DEM: DEM for background plotting
%            - T.conn: Connection nodes
%            - T.int: Interpolated nodes (used when T.conn has only one node)
%            - T.bl: Base layer coordinates (x, y)
%
%   Example:
%       plotT(T, 'MarkerSize', 20, 'LineWidth', 2);

    % Default parameters
    p = inputParser;
    addRequired(p, 'T');
    addParameter(p, 'MarkerSize', 15, @(v) isnumeric(v) && v > 0);
    addParameter(p, 'LineWidth', 1, @(v) isnumeric(v) && v > 0);
    
    parse(p, T, varargin{:});
    
    markerSize = p.Results.MarkerSize;
    lineWidth = p.Results.LineWidth;

    % Loop through both connection sets
    for i1 = 1:2
        all_x = []; 
        all_y = [];
        end_x = [];
        end_y = [];

        for i2 = 1:numel(T.x)
            for i3 = 1:numel(T.int{i1}(i2).x)

                if ~any(T.conn{i1}(i2).x{i3}), continue, end
                x_data = T.conn{i1}(i2).x{i3};
                y_data = T.conn{i1}(i2).y{i3};

                % If only one node in T.conn, use end members from T.int
                if numel(x_data) == 1
                    x_data = [x_data; T.int{i1}(i2).x{i3}(1); T.int{i1}(i2).x{i3}(end)];
                    y_data = [y_data; T.int{i1}(i2).y{i3}(1); T.int{i1}(i2).y{i3}(end)];
                end

                if numel(x_data) > 1
                    plot(x_data, y_data, 'Color', [0 0 0 .5], 'LineWidth', lineWidth); hold on
                    all_x = [all_x; x_data(:)];
                    all_y = [all_y; y_data(:)];
                    end_x = [end_x; x_data(end)];
                    end_y = [end_y; y_data(end)];
                else
                    scatter(x_data, y_data, 'ko', 'LineWidth', lineWidth); hold on
                    all_x = [all_x; x_data(:)];
                    all_y = [all_y; y_data(:)];
                    end_x = [end_x; x_data(end)];
                    end_y = [end_y; y_data(end)];
                end
            end
        end

        % Plot all scatter points
        if strcmp(T.type,'geometric')
            scatter(all_x, all_y, markerSize/2, [.6 .6 .8], 'filled'); hold on
        end

        % Plot end nodes with distinct colors
        scatter(end_x, end_y, markerSize + 10, [1 .3 .3] * (i1 == 1) + [.3 .3 1] * (i1 == 2), 'filled'); hold on
    end

    % Plot the base layer (T.bl)
    plot(T.x, T.y, '.-', 'color', 'w', 'LineWidth', 4);

end
