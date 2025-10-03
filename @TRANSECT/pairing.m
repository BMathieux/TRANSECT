function T = pairing(T, varargin)
%PAIRING Pair nodes of a TRANSECT object
%
%   T = PAIRING(T,x,y) pairs nodes of a TRANSECT object T with coordinates
%   (x,y). The function searches for intersections or proximity matches
%   between transect cross-sections and the input coordinates.
%
%   T = PAIRING(T,mask) pairs nodes of T using a logical or numeric grid 
%   mask (same size as T.DEM.Z). In this case, intersections are defined by
%   the region provided by mask, which could represent e.g. a floodplain or
%   geomorphic domain.
%
%   This updates T.int and T.conn with interpolated transect nodes and
%   connection nodes, and computes geometric statistics (width, height, 
%   slopes, elevations) stored in T.stats.
%
%   Name-Value options:
%       'type'         - 'both' (require connections on both sides) 
%                        or 'any' (accept one-sided, default).
%       'connectivity' - Neighborhood definition for mask search: 
%                        'D1','D8','D16' (default: 'D8').
%       'verbose'      - true/false (default: false).
%       'parallel'     - true/false (default: false).
%       'maxlength'    - true/false (default: false), truncate to shortest 
%                        side length.
%       'gap'          - scalar (default: 1), tolerance for missing nodes.
%       'order'        - 'first' or 'last' (default: 'first'), used only 
%                        when (x,y) input is given.
%
%   Example
%       % Pair transects with coordinates
%       T = pairing(T,x,y,'order','last');
%
%       % Pair transects with a mask
%       T = pairing(T,mask,'type','both','connectivity','D16');
%
%   See also TRANSECT

    if ~isa(T,'TRANSECT')
        error('Input must be a TRANSECT object.')
    end
    
    % detect input mode
    if numel(varargin)==1 && isequal(size(varargin{1}),size(T.DEM.Z))
        mode='mask';       % mask same size as DEM
    elseif numel(varargin)>=2 && isnumeric(varargin{1}) && isnumeric(varargin{2}) ...
            && isequal(size(varargin{1}),size(varargin{2}))
        mode='intersect';  % coordinate vectors
    else
        error('Second input must be either (x,y) or mask matching DEM size.')
    end
    
    % parser
    p = inputParser;
    addRequired(p,'T');
    switch mode
        case 'mask'
            isMask=@(v)(islogical(v)||isnumeric(v))&&isequal(size(v),size(T.DEM.Z));
            addRequired(p,'mask',isMask);
        case 'intersect'
            addRequired(p,'x',@isnumeric);
            addRequired(p,'y',@isnumeric);
            addParameter(p,'order','first',@(v)ismember(v,{'first','last'}));
    end
    addParameter(p,'type','any',@(v)ismember(v,{'both','any'}));
    addParameter(p,'connectivity','D8',@(v)ismember(v,{'D1','D8','D16'}));
    addParameter(p,'verbose',false,@islogical);
    addParameter(p,'parallel',false,@islogical);
    addParameter(p,'maxlength',false,@islogical);
    addParameter(p,'gap',1,@(v)isnumeric(v)&&isscalar(v)&&v>=0);
    parse(p,T,varargin{:});
    r = p.Results;
    
    % extract parameters
    tp   = r.type;
    con  = r.connectivity;
    vb   = r.verbose;
    parF = r.parallel;
    mlen = r.maxlength;
    gap  = r.gap;
    if strcmp(mode,'mask')
        mask = r.mask; ord = [];
    else
        x = r.x(:); y = r.y(:); ord = r.order;
    end
    
    % send DEM/mask/x,y to GPU if requested
    if parF
        T.DEM.Z      = gpuArray(T.DEM.Z);
        T.DEM.refmat = gpuArray(T.DEM.refmat);
        if strcmp(mode,'mask')
            mask = gpuArray(mask);
        else
            x = gpuArray(x); y = gpuArray(y);
        end
    end
    
    % main loop over transect nodes
    nC = numel(T.x);
    if vb, PB = ProgressBar(nC,'taskname','Pairing...','ui','cli'); end
    
    % precompute intersections if coordinate mode
    if strcmp(mode,'intersect')
        [X,Y]=refmat2XY(T.DEM.refmat,T.DEM.size); X=X(:); Y=Y(:);
        dx=X(2)-X(1); dy=Y(2)-Y(1);
        ix1=round((x-X(1))/dx+1); ix2=round((y-Y(1))/dy+1);
        ob=ix1>T.DEM.size(2)|ix1<1|ix2>T.DEM.size(1)|ix2<1;
        if any(ob)
            warning('TopoToolbox:outsidegrid','Some pts out-of-bound.');
            x(ob)=[]; y(ob)=[];
        end
        if isempty(x), warning('No valid pts.'); return, end
        ind=coord2ind(T.DEM,x,y);
    else
        ind=[];
    end
    
    % pairing (parallel or serial)
    if parF
        Ttmp=cell(nC,1);
        parfor i2=1:nC
            if strcmp(mode,'mask')
                Ttmp{i2}=pair(T,i2,ind,con,tp,'mask','mask',mask,'gap',gap);
            else
                Ttmp{i2}=pair(T,i2,ind,con,tp,'intersect','order',ord,'gap',gap);
            end
            if vb, PB.count(); end
        end
        for i2=1:nC
            for i1=1:2
                T.int{i1}(i2)=Ttmp{i2}.int{i1}(i2);
                T.conn{i1}(i2)=Ttmp{i2}.conn{i1}(i2);
            end
        end
    else
        for i2=1:nC
            if strcmp(mode,'mask')
                T=pair(T,i2,ind,con,tp,'mask','mask',mask,'gap',gap);
            else
                T=pair(T,i2,ind,con,tp,'intersect','order',ord,'gap',gap);
            end
            if vb, PB.count(); end
        end
    end
    
    % post-processing for mask mode
    if strcmp(mode,'mask')
        bw=bwdist(~mask)*T.DEM.cellsize*2;
        T=mask_validation(T,con,bw,'maxlength',mlen);
    end
    
    % gather GPU arrays back to CPU
    if parF
        T.DEM.Z      = gather(T.DEM.Z);
        T.DEM.refmat = gather(T.DEM.refmat);
    end
    
    % compute stats
    T = T.calcStats();


end

function T = pair(T, i2, ind, con, tp, met, varargin)
%PAIR Pair transect paths with intersections or mask constraints
%   T = PAIR(T,i2,ind,con,tp,met,...) updates the TRANSECT object by 
%   identifying valid connection nodes along interpolated transect lines.
%   Supports two methods:
%       - 'intersect': finds nearest matches with external (x,y) points
%       - 'mask':      validates path based on logical/numeric mask
%
%   Inputs:
%       T    - TRANSECT object
%       i2   - index of the base transect node
%       ind  - linear indices of intersect points (if method='intersect')
%       con  - connectivity ('D1','D8','D16')
%       tp   - type ('both' requires both sides, 'any' requires at least one)
%       met  - method string ('mask' or 'intersect')
%
%   Name-Value pairs:
%       'gap'   - allowed gap length in steps (default=1)
%       'order' - intersect order preference ('first'|'last', default='first')
%       'mask'  - logical/numeric mask (same size as DEM)
%
%   Output:
%       T - updated TRANSECT object

    % parse options
    p = inputParser;
    addParameter(p,'gap',1,@(v)isnumeric(v)&&isscalar(v)&&v>=0);
    addParameter(p,'order','first',@(v)ismember(v,{'first','last'}));
    addParameter(p,'mask',[],@(v)islogical(v)||isnumeric(v));
    parse(p,varargin{:});
    gap  = p.Results.gap;
    ord  = p.Results.order;
    mask = p.Results.mask;

    % neighborhood offsets for connectivity
    switch con
        case 'D1'
            dr=0; dc=0;
        case 'D8'
            dr=[-1 -1 -1 0 0 0 1 1 1];
            dc=[-1 0 1 -1 0 1 -1 0 1];
        case 'D16'
            dr=[-2 -2 -2 -2 -1 -1 -1 -1 0 0 0 0 0 1 1 1 1 2 2 2 2];
            dc=[-2 -1 0 1 -2 -1 0 1 -2 -1 0 1 2 -2 -1 0 1 -2 -1 0 1];
        otherwise
            error('Unknown connectivity: %s',con);
    end

    [NR,NC]=size(T.DEM.Z);

    % loop over both transect sides
    for i1=1:2
        for i3=1:numel(T.int{i1}(i2).x)
            if numel(T.int{i1}(i2).x{i3})==1, continue, end

            % interpolate line more densely
            xi=T.int{i1}(i2).x{i3}(:);
            yi=T.int{i1}(i2).y{i3}(:);
            valid=~(isnan(xi)|isnan(yi)); xi=xi(valid); yi=yi(valid);
            if numel(xi)<2, continue, end
            d=[0; cumsum(hypot(diff(xi),diff(yi)))];
            di=0:T.DEM.cellsize/5:d(end);
            xii=interp1(d,xi,di,'linear');
            yii=interp1(d,yi,di,'linear');
            ixii=coord2ind(T.DEM,xii,yii);

            % branch by method
            if strcmp(met,'mask')
                % validate path inside mask allowing short gaps
                msk=mask(ixii);
                n=numel(msk); last_end=0; curr_gap=0;
                for k=1:n
                    if msk(k)
                        curr_gap=0;
                    elseif curr_gap<gap
                        curr_gap=curr_gap+1;
                    else
                        break
                    end
                    last_end=k;
                end
                if last_end>0, p=last_end; else, p=[]; end
                nd=0; IXn=[];
            else
                % neighbor search for intersection indices
                [rA,cA]=ind2sub([NR,NC],ixii);
                R=rA(2:end); C=cA(2:end); nC=numel(R);
                Rn=repmat(R,1,numel(dr))+repmat(dr,nC,1);
                Cn=repmat(C,1,numel(dc))+repmat(dc,nC,1);
                vM=(Rn>=1)&(Rn<=NR)&(Cn>=1)&(Cn<=NC);
                nbA=nan(size(Rn)); nbA(vM)=sub2ind([NR,NC],Rn(vM),Cn(vM));
                dN0=hypot(dr,dc)*T.DEM.cellsize;
                dNA=repmat(dN0,nC,1); mskA=ismember(nbA,ind);
                cf=any(mskA,2); dNA(~mskA)=Inf;
                mDC=min(dNA,[],2); cK=find(cf)+1; cD=mDC(cf);

                if isempty(cK)
                    p=[];
                else
                    if strcmp(ord,'first')
                        cl=cumsum([1; diff(cK)>gap]); 
                        fc=cK(cl==cl(1)); cD2=cD(cl==cl(1));
                    else
                        cl=cumsum([1; diff(cK)>gap]); 
                        fc=cK(cl==cl(end)); cD2=cD(cl==cl(end));
                    end
                    [~,im]=min(cD2); p=fc(im);
                    ri=p-1; dN_row=dNA(ri,:); mR=mskA(ri,:);
                    dN_row(~mR)=Inf; [nd,ci]=min(dN_row);
                    IXn=nbA(ri,ci);
                end
            end

            % update transect and connection info
            if ~isempty(p)
                xS=xii(1:p)'; yS=yii(1:p)';
                dS=[0; cumsum(hypot(diff(xS),diff(yS)))];
                if nd==0||strcmp(con,'D1')
                    T.int{i1}(i2).x{i3}=xS;
                    T.int{i1}(i2).y{i3}=yS;
                    T.int{i1}(i2).ix{i3}=ixii(1:p);
                    T.int{i1}(i2).z{i3}=T.DEM.Z(ixii(1:p));
                    T.int{i1}(i2).d{i3}=dS;
                else
                    [xn,yn]=ind2coord(T.DEM,IXn); zn=T.DEM.Z(IXn);
                    dnS=dS(end)+T.DEM.cellsize;
                    T.int{i1}(i2).x{i3}=[xS;xn];
                    T.int{i1}(i2).y{i3}=[yS;yn];
                    T.int{i1}(i2).ix{i3}=[ixii(1:p);IXn];
                    T.int{i1}(i2).z{i3}=[T.DEM.Z(ixii(1:p));zn];
                    T.int{i1}(i2).d{i3}=[dS;dnS];
                end
                % update connections
                existing=T.conn{i1}(i2).ix{i3};
                comp=ixii(1:p);
                if ~isempty(existing)
                    c0=existing(ismember(existing,comp));
                else
                    c0=[];
                end
                if ~isempty(IXn)&&~ismember(IXn,c0), c0=[c0;IXn]; end
                [cx,cy]=ind2coord(T.DEM,c0);
                T.conn{i1}(i2).x{i3}=cx;
                T.conn{i1}(i2).y{i3}=cy;
                T.conn{i1}(i2).ix{i3}=c0;
                T.conn{i1}(i2).z{i3}=T.DEM.Z(c0);
                T.conn{i1}(i2).d{i3}=[0; cumsum(hypot(diff(cx),diff(cy)))];
            else
                % empty if nothing found
                T.int{i1}(i2).x{i3}=[]; T.int{i1}(i2).y{i3}=[];
                T.int{i1}(i2).ix{i3}=[]; T.int{i1}(i2).z{i3}=[]; T.int{i1}(i2).d{i3}=[];
                T.conn{i1}(i2).x{i3}=[]; T.conn{i1}(i2).y{i3}=[];
                T.conn{i1}(i2).ix{i3}=[]; T.conn{i1}(i2).z{i3}=[]; T.conn{i1}(i2).d{i3}=[];
            end
        end
    end

    % validate that required side(s) exist
    v1=any(cellfun(@(x)~isempty(x),T.int{1}(i2).x));
    v2=any(cellfun(@(x)~isempty(x),T.int{2}(i2).x));
    if strcmp(tp,'both'), vp=v1&&v2; else, vp=v1||v2; end
    if ~vp
        % reset both sides if validation fails
        for i1=1:2
            for i3=1:numel(T.conn{i1}(i2).x)
                T.int{i1}(i2).x{i3}=[]; T.int{i1}(i2).y{i3}=[];
                T.int{i1}(i2).ix{i3}=[]; T.int{i1}(i2).z{i3}=[]; T.int{i1}(i2).d{i3}=[];
                T.conn{i1}(i2).x{i3}=[]; T.conn{i1}(i2).y{i3}=[];
                T.conn{i1}(i2).ix{i3}=[]; T.conn{i1}(i2).z{i3}=[]; T.conn{i1}(i2).d{i3}=[];
            end
        end
    end
end


function T = mask_validation(T, con, bw, varargin)
%MASK_VALIDATION Adjust transect connections using mask distances
%   T = MASK_VALIDATION(T,con,bw,...) refines transect connections by 
%   comparing path lengths to mask-based distances. Ensures consistent 
%   maximum length when traversing across a mask.
%
%   Inputs:
%       T   - TRANSECT object
%       con - connectivity ('D1','D8','D16')
%       bw  - distance transform of mask
%
%   Name-Value pairs:
%       'maxlength' - logical, if true forces both sides to same length

    p=inputParser;
    addParameter(p,'maxlength',false,@islogical);
    parse(p,varargin{:});
    mlen=p.Results.maxlength;

    % connectivity offsets
    switch con
        case 'D1'
            dr=0; dc=0;
        case 'D8'
            dr=[-1 -1 -1 0 0 0 1 1 1];
            dc=[-1 0 1 -1 0 1 -1 0 1];
        case 'D16'
            dr=[-2 -2 -2 -2 -1 -1 -1 -1 0 0 0 0 0 1 1 1 1 2 2 2 2];
            dc=[-2 -1 0 1 -2 -1 0 1 -2 -1 0 1 2 -2 -1 0 1 -2 -1 0 1];
        otherwise
            error('Unknown connectivity: %s',con);
    end

    nB=numel(T.x);
    for i2=1:nB
        id=[T.int{1}(i2).ix{1}; T.int{2}(i2).ix{1}];
        if isempty(id), continue, end

        % compute maximum neighbor distances
        [r,c]=ind2sub(size(bw),id); mb=zeros(size(id));
        for k=1:numel(id)
            r0=r(k); c0=c(k); nr0=r0+dr; nc0=c0+dc;
            v=(nr0>=1)&(nr0<=size(bw,1))&(nc0>=1)&(nc0<=size(bw,2));
            nId=sub2ind(size(bw),nr0(v),nc0(v));
            mb(k)=max(bw(nId));
        end
        [mbv,ia]=max(mb);

        % determine balance index (which side to trim)
        mbwI=nan(2,1);
        [~,ib]=intersect(T.int{1}(i2).ix{1},id(ia));
        if ~isempty(ib)
            mbwI(1)=T.int{1}(i2).d{1}(ib)+mbv;
            mbwI(2)=mbv-T.int{1}(i2).d{1}(ib);
        else
            [~,ib]=intersect(T.int{2}(i2).ix{1},id(ia));
            mbwI(2)=T.int{2}(i2).d{1}(ib)+mbv;
            mbwI(1)=mbv-T.int{2}(i2).d{1}(ib);
        end

        % trim each side according to mask distances
        lens=zeros(2,1);
        for i1=1:2
            for i3=1:numel(T.int{i1}(i2).x)
                if isempty(T.int{i1}(i2).x{i3}), continue, end
                xii=T.int{i1}(i2).x{i3}; yii=T.int{i1}(i2).y{i3};
                di=T.int{i1}(i2).d{i3};
                if isempty(di)||isnan(mbwI(i1))
                    % clear invalid
                    T.int{i1}(i2).x{i3}=[]; T.int{i1}(i2).y{i3}=[];
                    T.int{i1}(i2).ix{i3}=[]; T.int{i1}(i2).z{i3}=[]; T.int{i1}(i2).d{i3}=[];
                    T.conn{i1}(i2).x{i3}=[]; T.conn{i1}(i2).y{i3}=[];
                    T.conn{i1}(i2).ix{i3}=[]; T.conn{i1}(i2).z{i3}=[]; T.conn{i1}(i2).d{i3}=[];
                    continue
                end
                [~,p]=min(abs(mbwI(i1)-di));
                if p<1||p>numel(xii)
                    % clear if out of range
                    T.int{i1}(i2).x{i3}=[]; T.int{i1}(i2).y{i3}=[];
                    T.int{i1}(i2).ix{i3}=[]; T.int{i1}(i2).z{i3}=[]; T.int{i1}(i2).d{i3}=[];
                    T.conn{i1}(i2).x{i3}=[]; T.conn{i1}(i2).y{i3}=[];
                    T.conn{i1}(i2).ix{i3}=[]; T.conn{i1}(i2).z{i3}=[]; T.conn{i1}(i2).d{i3}=[];
                    continue
                end
                lens(i1)=p;
                T.int{i1}(i2).x{i3}=xii(1:p);
                T.int{i1}(i2).y{i3}=yii(1:p);
                T.int{i1}(i2).ix{i3}=coord2ind(T.DEM,xii(1:p),yii(1:p));
                T.int{i1}(i2).z{i3}=T.DEM.Z(T.int{i1}(i2).ix{i3});
                T.int{i1}(i2).d{i3}=[0;cumsum(hypot(diff(xii(1:p)),diff(yii(1:p))))];
            end
        end

        % enforce equal lengths if maxlength is true
        if mlen
            nlim=min(lens(lens>0));
            for i1=1:2
                for i3=1:numel(T.int{i1}(i2).x)
                    if isempty(T.int{i1}(i2).x{i3}), continue, end
                    p=min(nlim,numel(T.int{i1}(i2).x{i3}));
                    T.int{i1}(i2).x{i3}=T.int{i1}(i2).x{i3}(1:p);
                    T.int{i1}(i2).y{i3}=T.int{i1}(i2).y{i3}(1:p);
                    T.int{i1}(i2).ix{i3}=T.int{i1}(i2).ix{i3}(1:p);
                    T.int{i1}(i2).z{i3}=T.int{i1}(i2).z{i3}(1:p);
                    T.int{i1}(i2).d{i3}=T.int{i1}(i2).d{i3}(1:p);
                    xi=[T.int{i1}(i2).x{i3}(1);T.int{i1}(i2).x{i3}(end)];
                    yi=[T.int{i1}(i2).y{i3}(1);T.int{i1}(i2).y{i3}(end)];
                    ix=unique(coord2ind(T.DEM,xi,yi),'stable');
                    [xi,yi]=ind2coord(T.DEM,ix);
                    T.conn{i1}(i2).x{i3}=xi; T.conn{i1}(i2).y{i3}=yi;
                    T.conn{i1}(i2).ix{i3}=ix; T.conn{i1}(i2).z{i3}=T.DEM.Z(ix);
                    T.conn{i1}(i2).d{i3}=[0;cumsum(hypot(diff(xi),diff(yi)))];
                end
            end
        else
            % otherwise just connect endpoints
            for i1=1:2
                for i3=1:numel(T.int{i1}(i2).x)
                    if isempty(T.int{i1}(i2).x{i3}), continue, end
                    % endpoints only
                    xi=[T.int{i1}(i2).x{i3}(1); T.int{i1}(i2).x{i3}(end)];
                    yi=[T.int{i1}(i2).y{i3}(1); T.int{i1}(i2).y{i3}(end)];
                    ix=unique(coord2ind(T.DEM,xi,yi),'stable');
                    [xi,yi]=ind2coord(T.DEM,ix);
                    T.conn{i1}(i2).x{i3}=xi;
                    T.conn{i1}(i2).y{i3}=yi;
                    T.conn{i1}(i2).ix{i3}=ix;
                    T.conn{i1}(i2).z{i3}=T.DEM.Z(ix);
                    T.conn{i1}(i2).d{i3}=[0; cumsum(hypot(diff(xi),diff(yi)))];
                end
            end
        end
    end
end

