function plot(T,varargin)
%PLOT Display a TRANSECT object
%
%   plot(T) displays the transect geometry, paired connections, and end 
%   nodes stored in the TRANSECT object T. The trunk baseline is plotted 
%   together with extracted connection paths.
%
%   plot(T,'Name',Value,...) customizes the appearance:
%       'MarkerSize'     - Scatter marker size (default: 15)
%       'LineWidth'      - Width of connection lines (default: 1)
%       'EndColors'      - 2x3 matrix of RGB colors for left/right ends
%                          Default: [1 .3 .3; .3 .3 1]
%       'LineColor'      - RGB triple (or RGBA) for transect lines
%                          Default: [0 0 0 0.5]
%       'BaselineColor'  - RGB triple for trunk baseline (default: [1 1 1])
%       'BaselineWidth'  - Line width of trunk baseline (default: 4)
%       'BaselineStyle'  - Line style of trunk baseline (default: '-')
%
% Example
%
%       plot(T,'BaselineColor',[0 0 0],'BaselineWidth',2,'BaselineStyle','--')
%
% See also: TRANSECT

    % enforce TRANSECT input
    if ~isa(T,'TRANSECT')
        error('Input must be a TRANSECT object.')
    end

    % parser
    p = inputParser;
    addRequired(p,'T');
    addParameter(p,'MarkerSize',15,@(v) isnumeric(v) && v>0);
    addParameter(p,'LineWidth',1,@(v) isnumeric(v) && v>0);
    addParameter(p,'EndColors',[1 .3 .3; .3 .3 1],@(v) isnumeric(v) && size(v,2)==3);
    addParameter(p,'LineColor',[0 0 0 0.5],@(v) isnumeric(v) && (numel(v)==3 || numel(v)==4));
    addParameter(p,'BaselineColor',[1 1 1],@(v) isnumeric(v) && numel(v)==3);
    addParameter(p,'BaselineWidth',4,@(v) isnumeric(v) && v>0);
    addParameter(p,'BaselineStyle','-',@(v) ischar(v) || isstring(v));
    parse(p,T,varargin{:});

    msz = p.Results.MarkerSize;
    lw  = p.Results.LineWidth;
    ec  = p.Results.EndColors;
    lc  = p.Results.LineColor;
    blc = p.Results.BaselineColor;
    blw = p.Results.BaselineWidth;
    bls = p.Results.BaselineStyle;

    hold on
    for i1 = 1:2
        all_x=[]; all_y=[]; end_x=[]; end_y=[];
        for i2 = 1:numel(T.x)
            for i3 = 1:numel(T.int{i1}(i2).x)
                if ~any(T.conn{i1}(i2).x{i3}), continue, end
                x_data = T.conn{i1}(i2).x{i3};
                y_data = T.conn{i1}(i2).y{i3};

                if numel(x_data)==1
                    x_data=[x_data; T.int{i1}(i2).x{i3}([1 end])];
                    y_data=[y_data; T.int{i1}(i2).y{i3}([1 end])];
                end

                if numel(x_data)>1
                    builtin('plot',x_data,y_data,'Color',lc,'LineWidth',lw);
                    all_x=[all_x; x_data(:)];
                    all_y=[all_y; y_data(:)];
                    end_x=[end_x; x_data(end)];
                    end_y=[end_y; y_data(end)];
                else
                    scatter(x_data,y_data,msz,'ko','LineWidth',lw);
                    all_x=[all_x; x_data(:)];
                    all_y=[all_y; y_data(:)];
                    end_x=[end_x; x_data(end)];
                    end_y=[end_y; y_data(end)];
                end
            end
        end

        if strcmp(T.type,'geometric')
            scatter(all_x,all_y,msz/2,[.6 .6 .8],'filled');
        end

        scatter(end_x,end_y,msz+10,ec(i1,:),'filled');
    end

    % baseline plot
    builtin('plot',T.x,T.y,'.-','Color',blc,'LineWidth',blw,'LineStyle',bls);

end
