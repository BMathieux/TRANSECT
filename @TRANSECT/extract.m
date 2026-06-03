function out = extract(T,typ,field,varargin)
%EXTRACT Extract data from a TRANSECT object
%
%   out = extract(T,'int', field, i1,i2,i3) extracts from interpolated paths
%   out = extract(T,'conn',field, i1,i2,i3) extracts from connection paths
%
%       field: 'x','y','z','d','ix','coords'
%       i1: side index (1 or 2)
%       i2: baseline node index
%       i3: path index
%
%   out = extract(T,'stats',field) extracts full stats field
%   out = extract(T,'stats',field,i2) extracts stats at baseline index i2
%
%       field: 'x','y','z','d','Z1','Z2','G1','G2','W','H'
%
%   out = extract(T,'base',field) extracts original baseline information
%       field: 'typ','obj','x','y','ix','coords'
%
%   out = extract(T,'pair',field) extracts paired-object information
%       field: 'typ','obj','bid','ok','opt'
%              if T.pair.obj stores coordinates also: 'x','y','ix'
%
%   Examples:
%       x  = extract(T,'int','x',1,5,2);
%       z  = extract(T,'conn','z',2,3,1);
%       c  = extract(T,'int','coords',1,2,1);
%
%       W  = extract(T,'stats','W');
%       Wi = extract(T,'stats','W',10);
%
%       bt  = extract(T,'base','typ');
%       bxy = extract(T,'base','coords');
%
%       xP  = extract(T,'pair','x');
%       yP  = extract(T,'pair','y');
%       ixP = extract(T,'pair','ix');
%       bid = extract(T,'pair','bid');

    if ~isa(T,'TRANSECT'), error('First input must be a TRANSECT object.'), end
    if ~ischar(typ) && ~isstring(typ), error('typ must be a string.'), end
    if ~ischar(field) && ~isstring(field), error('field must be a string.'), end
    typ=char(typ); field=char(field);

    if ~ismember(typ,{'int','conn','stats','pair','base'})
        error('typ must be ''int'', ''conn'', ''stats'', ''pair'', or ''base''.')
    end

    switch typ

        case {'int','conn'}

            if numel(varargin)~=3
                error('Usage: extract(T,''%s'',field,i1,i2,i3)',typ)
            end
            i1=varargin{1}; i2=varargin{2}; i3=varargin{3};

            if ~ismember(i1,[1 2]), error('i1 must be 1 or 2 (side index).'), end
            if ~isscalar(i2) || i2<1 || i2>numel(T.(typ){i1}), error('i2 out of range.'), end
            if ~isscalar(i3) || i3<1 || i3>numel(T.(typ){i1}(i2).x), error('i3 out of range.'), end
            if ~ismember(field,{'x','y','z','d','ix','coords'})
                error('field must be ''x'',''y'',''z'',''d'',''ix'', or ''coords''.')
            end

            if strcmp(field,'coords')
                x=T.(typ){i1}(i2).x{i3};
                y=T.(typ){i1}(i2).y{i3};
                if isfield(T.(typ){i1}(i2),'z') && numel(T.(typ){i1}(i2).z)>=i3 && ~isempty(T.(typ){i1}(i2).z{i3})
                    z=T.(typ){i1}(i2).z{i3};
                    out=[x(:) y(:) z(:)];
                else
                    out=[x(:) y(:)];
                end
            else
                out=T.(typ){i1}(i2).(field){i3};
            end

        case 'stats'

            if numel(varargin)~=0 && numel(varargin)~=1
                error('Usage: extract(T,''stats'',field) or extract(T,''stats'',field,i2)')
            end
            if ~isfield(T.stats,field), error('Unknown stats field: %s',field), end

            v=T.stats.(field);
            if isempty(varargin)
                out=v;
            else
                i2=varargin{1};
                if iscell(v), n=numel(v); else, n=numel(v); end
                if ~isscalar(i2) || i2<1 || i2>n, error('i2 out of range.'), end
                if iscell(v), out=v{i2}; else, out=v(i2); end
            end

        case 'base'

            if ~isempty(varargin)
                error('Usage: extract(T,''base'',field)')
            end
            if ~isprop(T,'base') || isempty(T.base)
                out=[];
                return
            end

            switch field
                case 'typ'
                    out=T.base.typ;
                case 'obj'
                    if isfield(T.base,'obj'), out=T.base.obj; else, out=[]; end
                case 'x'
                    [out,~]=basecoords(T);
                case 'y'
                    [~,out]=basecoords(T);
                case 'ix'
                    out=baseix(T);
                case 'coords'
                    [x,y]=basecoords(T);
                    out=[x(:) y(:)];
                otherwise
                    error('Unknown base field: %s',field)
            end

        case 'pair'

            if ~isempty(varargin)
                error('Usage: extract(T,''pair'',field)')
            end
            if isempty(T.pair)
                out=[];
                return
            end

            if ismember(field,{'typ','obj','bid','ok','opt'})
                if isfield(T.pair,field)
                    out=T.pair.(field);
                else
                    out=[];
                end
                return
            end

            if isstruct(T.pair.obj) && ismember(field,{'x','y','ix'})
                if isfield(T.pair.obj,field)
                    out=T.pair.obj.(field);
                else
                    out=[];
                end
            else
                error('Unknown pair field: %s',field)
            end
    end
end

function [x,y] = basecoords(T)
%BASECOORDS Return original baseline coordinates from T.base

    switch T.base.typ
        case 'xy'
            x=T.base.x(:); y=T.base.y(:);
        case 'STREAMobj'
            x=T.base.obj.x(:); y=T.base.obj.y(:);
        case 'DIVIDEobj'
            ix=T.base.obj.IX(:);
            [x,y]=ind2coord(T.base.obj,ix);
            x=x(:); y=y(:);
        otherwise
            error('Unknown base type: %s',T.base.typ)
    end
end

function ix = baseix(T)
%BASEIX Return original baseline indices from T.base

    switch T.base.typ
        case 'xy'
            ix=coord2ind(T.DEM,T.base.x(:),T.base.y(:));
        case 'STREAMobj'
            ix=T.base.obj.IXgrid(:);
        case 'DIVIDEobj'
            ix=T.base.obj.IX(:);
        otherwise
            error('Unknown base type: %s',T.base.typ)
    end
end
