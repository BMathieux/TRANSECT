function plot(T,varargin)
%PLOT Display a TRANSECT object
%
%   plot(T) displays the transect geometry, paired connections, and end
%   nodes stored in the TRANSECT object T. The active transect baseline is
%   inferred from the first node of available transect paths.
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
%       'Vectors'        - Transect paths to plot: 'all' (default) or 'middle'
%       'Path'           - Path geometry to plot: 'conn' (default) or 'int'
%
% Example
%
%       plot(T,'BaselineColor',[0 0 0],'BaselineWidth',2,'BaselineStyle','--')
%       plot(T,'Vectors','all')
%       plot(T,'Path','int')
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
    addParameter(p,'Vectors','all',@(v) any(validatestring(v,{'all','middle'})));
    addParameter(p,'Path','conn',@(v) any(validatestring(v,{'conn','int'})));
    parse(p,T,varargin{:});

    msz=p.Results.MarkerSize; lw=p.Results.LineWidth; ec=p.Results.EndColors;
    lc=p.Results.LineColor; blc=p.Results.BaselineColor; blw=p.Results.BaselineWidth; bls=p.Results.BaselineStyle;
    vec=validatestring(p.Results.Vectors,{'all','middle'});
    pth=validatestring(p.Results.Path,{'conn','int'});
    
    hold on
    
    for i1=1:2
        if numel(T.int)<i1 || isempty(T.int{i1})
            continue
        end

        lx={}; ly={}; xc={}; yc={}; end_x=[]; end_y=[];

        for i2=1:numel(T.int{i1})
            id=validpaths(T,i1,i2,pth);

            if isempty(id)
                continue
            end

            if strcmp(vec,'middle')
                id=id(ceil(numel(id)/2));
            end

            for i3=id
                [xd,yd]=getpath(T,i1,i2,i3,pth);

                if isempty(xd) || numel(xd)<2 || all(isnan(xd))
                    continue
                end

                lx{end+1,1}=[xd(:);nan];
                ly{end+1,1}=[yd(:);nan];

                end_x=[end_x;xd(end)];
                end_y=[end_y;yd(end)];

                if strcmp(T.type,'geometric') && ...
                        numel(T.conn)>=i1 && ...
                        i2<=numel(T.conn{i1}) && ...
                        i3<=numel(T.conn{i1}(i2).x) && ...
                        ~isempty(T.conn{i1}(i2).x{i3})

                    xci=T.conn{i1}(i2).x{i3};
                    yci=T.conn{i1}(i2).y{i3};

                    if ~isempty(xci) && ~all(isnan(xci))
                        xc{end+1,1}=xci(:);
                        yc{end+1,1}=yci(:);
                    end
                end
            end
        end

        if ~isempty(lx)
            builtin('plot',vertcat(lx{:}),vertcat(ly{:}),'Color',lc,'LineWidth',lw);
        end

        if strcmp(T.type,'geometric') && ~isempty(xc)
            scatter(vertcat(xc{:}),vertcat(yc{:}),msz/2,[.6 .6 .8],'filled');
        end

        if ~isempty(end_x)
            scatter(end_x,end_y,msz+10,ec(i1,:),'filled');
        end
    end

    % baseline plot
    [xb,yb]=activebase(T);
    if ~isempty(xb)
        builtin('plot',xb,yb,'.-','Color',blc,'LineWidth',blw,'LineStyle',bls);
    end

end

function id = validpaths(T,i1,i2,pth)
%VALIDPATHS Return valid path indices for one baseline node

    id=[];

    if strcmp(pth,'conn') && numel(T.conn)>=i1 && i2<=numel(T.conn{i1})
        S=T.conn{i1}(i2);
    else
        S=T.int{i1}(i2);
    end

    if ~isfield(S,'x') || isempty(S.x)
        return
    end

    for k=1:numel(S.x)
        x=S.x{k};
        if ~isempty(x) && numel(x)>=2 && ~all(isnan(x))
            id=[id k];
        end
    end
end

function [x,y] = getpath(T,i1,i2,i3,pth)
%GETPATH Return connection or interpolated path

    x=[]; y=[];

    if strcmp(pth,'conn') && ...
            numel(T.conn)>=i1 && ...
            i2<=numel(T.conn{i1}) && ...
            i3<=numel(T.conn{i1}(i2).x) && ...
            ~isempty(T.conn{i1}(i2).x{i3})

        x=T.conn{i1}(i2).x{i3};
        y=T.conn{i1}(i2).y{i3};
    elseif i3<=numel(T.int{i1}(i2).x)
        x=T.int{i1}(i2).x{i3};
        y=T.int{i1}(i2).y{i3};
    end
end

function [xb,yb] = activebase(T)
%ACTIVEBASE Reconstruct active baseline from transect path starts

    if isempty(T.int)
        xb=[]; yb=[];
        return
    end

    n=numel(T.int{1});
    xb=nan(n,1); yb=nan(n,1);

    for i=1:n
        for i1=1:2
            if i>numel(T.int{i1})
                continue
            end

            k=find(~cellfun(@isempty,T.int{i1}(i).x),1,'first');

            if ~isempty(k)
                xb(i)=T.int{i1}(i).x{k}(1);
                yb(i)=T.int{i1}(i).y{k}(1);
                break
            end
        end
    end

    I=~isnan(xb);
    xb=xb(I); yb=yb(I);
end
