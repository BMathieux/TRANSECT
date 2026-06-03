function [w,wi] = maskWidth(DEM,x,y,mask,ite)
%MASKWIDTH Estimate TRANSECT half-width from all cells inside a mask
%
%   [w,wi] = maskWidth(DEM,x,y,mask,ite) estimates the TRANSECT width
%   parameter w so that the effective search width w*ite includes all
%   nonzero cells of mask.
%
%   w  : width parameter passed to TRANSECT
%   wi : effective width in DEM cells, equal to w*ite

if isa(mask,'GRIDobj')
    mask = mask.Z;
end

mask = mask ~= 0 & ~isnan(mask);

ix = coord2ind(DEM,x,y);
ix = ix(~isnan(ix));

M = false(DEM.size);
M(ix) = true;

D = bwdist(M);
d = D(mask);

if isempty(d)
    w = 1;
    wi = ite;
    return
end

w   = ceil(max(d(:))/ite)*1.02;
wi  = w*ite;

end
