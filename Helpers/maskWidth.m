function [w,wi] = maskWidth(mask,ite)
%MASKWIDTH Compute transect width and DEM padding from a binary mask
%
%   [w,wi] = MASKWIDTH(mask,ite) estimates the minimal transect half-width (w)
%   in iteration units and DEM padding (wi) in pixels, such that transects
%   traced across a TRANSECT object will fully span the region defined by mask.
%
%   Inputs:
%     mask - logical/numeric grid for pairing 
%            (e.g. floodplain or basin mask)
%     ite  - iteration step used in TRANSECT construction
%
%   Outputs:
%     w    - half-width in TRANSECT per iteration (value passed to 'w')
%     wi   - half-width in DEM pixels (padding applied to DEM)
%
%   Example
%     % Pair transects with a floodplain mask
%     ite = 5;
%     [w,wi] = maskWidth(FP,ite);
%     T = TRANSECT(DEM,S.x,S.y,'w',w,'ite',ite);
%     T = pairing(T,FP);
%
%     % Pair transects with a basin outline
%     DB = drainagebasins(FD,S.IXgrid(end));
%     mask = logical(DB.Z);
%     ite = 5;
%     [w,wi] = maskWidth(mask,ite);
%     T = TRANSECT(DEM,S.x,S.y,'w',w,'ite',ite);
%     T = pairing(T,mask);

    % outline of mask
    B = bwboundaries(mask,8,'noholes');
    b = B{1}; r = b(:,1); c = b(:,2);
    
    % fast bounds in pixel units
    dmin = max(pdist2([c r],find(mask)','euclidean','Smallest',1)); % lower bound
    dx   = max(c)-min(c); dy = max(r)-min(r);
    dmax = max(dx,dy); 
    lo   = ceil(dmin); hi = ceil(dmax);
    
    % dilation test
    chk = @(w) all(interp2(double(imdilate(mask,strel('disk',w))),c,r)>0.5);
    
    % binary search
    if chk(lo)
        w = lo;
    else
        while lo<hi
            mid = floor((lo+hi)/2);
            if chk(mid), hi=mid; else, lo=mid+1; end
        end
        w = hi;
    end
    
    % outputs
    wi = ceil(w*1.1);
    w  = ceil(wi/ite);
end
