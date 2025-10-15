function out = extractT(T,i1,i2,i3,field,varargin)
%extractT Extract data from a TRANSECT object
%
%   out = extractT(T,i1,i2,i3,field) extracts the requested field from the
%   'int' group of the TRANSECT object T at indices (i1,i2,i3).
%
%   out = extractT(T,i1,i2,i3,field,'group','conn') extracts from 'conn'
%   instead of 'int'.
%
%   Inputs:
%       T      - TRANSECT object
%       i1     - side index (1 or 2)
%       i2     - transect index
%       i3     - path index
%       field  - string: 'x','y','z','d','ix','coords'
%
%   Optional Name-Value:
%       'group' - 'int' (default) or 'conn'
%
%   Output:
%       out    - array corresponding to the requested field
%
%   Example:
%       x = extractT(T,1,5,2,'x');                 
%       z = extractT(T,2,3,1,'z','group','conn');  
%       c = extractT(T,1,2,1,'coords');

    % parse optional inputs
    p = inputParser;
    addParameter(p,'group','int',@(v) ismember(v,{'int','conn'}));
    parse(p,varargin{:});
    grp = p.Results.group;

    % input checks
    if ~isa(T,'TRANSECT')
        error('First input must be a TRANSECT object.')
    end
    if ~ismember(i1,[1 2])
        error('i1 must be 1 or 2 (side index).')
    end
    if i2<1 || i2>numel(T.x)
        error('i2 out of range.')
    end
    if i3<1 || i3>numel(T.(grp){i1}(i2).x)
        error('i3 out of range.')
    end
    if ~ismember(field,{'x','y','z','d','ix','coords'})
        error('field must be ''x'',''y'',''z'',''d'',''ix'', or ''coords''.')
    end

    % extract data
    if strcmp(field,'coords')
        x = T.(grp){i1}(i2).x{i3};
        y = T.(grp){i1}(i2).y{i3};
        if isfield(T.(grp){i1}(i2),'z') && numel(T.(grp){i1}(i2).z)>=i3 && ~isempty(T.(grp){i1}(i2).z{i3})
            z = T.(grp){i1}(i2).z{i3};
            out = [x(:) y(:) z(:)];
        else
            out = [x(:) y(:)];
        end
    else
        out = T.(grp){i1}(i2).(field){i3};
    end
end
