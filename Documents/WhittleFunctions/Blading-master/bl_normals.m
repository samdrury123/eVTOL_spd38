function [n_bl,si,sj,Fn] = bl_normals(xyz,is_blade,smooth_normals)
% BL_NORMALS  Calculate blade normals and parametric coordinates
%
%   [n_bl,si,sj,F] = BL_NORMALS(xyz,is_blade,smooth_normals)
%
%   xyz - 3 dimensional blade cartesian coordinate array
%   is_blade - 0 or 1 to indicate whether to repeat leading edge coordinates
%   smooth_normals - 0 or 1 to smooth normals in streamwise direction
%   n2 - 3 dimensional blade normal array
%   si - 2 dimensional parametric measure of around section distance
%   sj - 2 dimensional parametric measure of spanwise distance
%   F - data structure of coordinates and normals in parametric coords

% Default to calculate normals for a blade
if exist('is_blade','var') == 0
    is_blade = 1;
end

% Default to not smooth
if exist('smooth_normals','var') == 0
    smooth_normals = 0;
end

% Get blade size
nj = size(xyz,2); ni = size(xyz,1);

% Calculate polar coordinates
xrrt = bl_cart2pol(xyz);

% Calculate plaid parametric coordinates on blade surfaces
si = repmat(reshape(linspace(0,1,ni),[ni 1]),[1 nj]);
sj = repmat(reshape(linspace(0,1,nj),[1 nj]),[ni 1]);

% Repeat coordinates around leading edge
if is_blade == 1
    m = 5;
    il = [si(end-m:end-1,:)-1 ; si ; si(2:m,:)+1];
    jl = [sj(end-m:end-1,:) ; sj ; sj(2:m,:)];
    xyzl = [xyz(end-m:end-1,:,:) ; xyz ; xyz(2:m,:,:)];
    xrrtl = [xrrt(end-m:end-1,:,:) ; xrrt ; xrrt(2:m,:,:)];
else
    il = si; jl = sj; xyzl = xyz; xrrtl = xrrt;
end

% Construct cartesian coordinate interpolants
Fn.x = griddedInterpolant(il,jl,xyzl(:,:,1),'cubic');
Fn.y = griddedInterpolant(il,jl,xyzl(:,:,2),'cubic');
Fn.z = griddedInterpolant(il,jl,xyzl(:,:,3),'cubic');
Fn.r = griddedInterpolant(il,jl,xrrtl(:,:,2));

% Maybe coordinates should be repeated before normal computation

% Cross product of grid vectors for blade normals
a2 = xyz(2:end,2:end,:) - xyz(1:end-1,1:end-1,:);
a1 = xyz(1:end-1,2:end,:) - xyz(2:end,1:end-1,:);
n1 = cross(a1,a2); n1 = n1 ./ repmat(sum(n1.^2,3).^0.5,[1 1 3]);

% Fill in normals around the edges
if is_blade == 1
    
    % Fill in the centre of the surface
    n_bl = zeros(size(xyz));
    n_bl(2:end-1,2:end-1,:) = 0.25*(n1(1:end-1,1:end-1,:) + n1(2:end,1:end-1,:) +...
        n1(1:end-1,2:end,:) + n1(2:end,2:end,:));

    % Fill in leading edge
    n_bl(1,2:end-1,:) = 0.25*(n1(1,1:end-1,:) + n1(end,1:end-1,:) +...
        n1(1,2:end,:) + n1(end,2:end,:));
    n_bl(end,2:end-1,:) = 0.25*(n1(1,1:end-1,:) + n1(end,1:end-1,:) +...
        n1(1,2:end,:) + n1(end,2:end,:));

    % Fill in bottom section
    n_bl(2:end-1,end,:) = 0.5 * (n1(1:end-1,end,:) + n1(2:end,end,:));
    n_bl(1,end,:) = 0.5 * (n1(1,end,:) + n1(end,end,:));
    n_bl(end,end,:) = 0.5 * (n1(1,end,:) + n1(end,end,:));

    % Fill in top section
    n_bl(:,1,:) = n_bl(:,2,:);
    
    % Repeat normals about leading edge
    nl = [n_bl(end-m:end-1,:,:) ; n_bl ; n_bl(2:m,:,:)];
    
else
    
    % Construct parametric coordinates for average normals
    si1 = repmat(0.5 * (si(1:end-1,1) + si(2:end,1)),[1 nj-1]);
    sj1 = repmat(0.5 * (sj(1,1:end-1) + sj(1,2:end)),[ni-1 1]);
    
    % Construct interponts for known normals
    nx = griddedInterpolant(si1,sj1,n1(:,:,1),'cubic','cubic');
    ny = griddedInterpolant(si1,sj1,n1(:,:,2),'cubic','cubic');
    nz = griddedInterpolant(si1,sj1,n1(:,:,3),'cubic','cubic');
    
    % Evaluate interpolants including edges
    n_bl = cat(3,nx(si,sj),ny(si,sj),nz(si,sj)); 
    nl = n_bl;
end

% Normalise normals
n_bl = n_bl ./ repmat(sum(n_bl.^2,3).^0.5,[1 1 3]);

% Smooth normals
if smooth_normals == 1
    n_bl = n_bl([4:-1:2 1:end end-1:-1:end-3],:,:);
    n_bl = reshape(smooth_2d(reshape(n_bl,[ni+6 nj*3])),[ni+6 nj 3]);
    n_bl = n_bl(4:end-3,:,:); nl = n_bl;
end

% Construct normal interpolants
Fn.nx = griddedInterpolant(il,jl,nl(:,:,1),'cubic');
Fn.ny = griddedInterpolant(il,jl,nl(:,:,2),'cubic');
Fn.nz = griddedInterpolant(il,jl,nl(:,:,3),'cubic');

% Find whether point is on side 1 or two based upon minimum and maximum x
% at each section
[~,i_min] = min(xyz(:,:,1),[],1); [~,i_max] = max(xyz(:,:,1),[],1);
u = zeros(ni,nj);
for n = 1:length(i_min)
    if i_min(n) > i_max(n)
        u(i_min(n):-1:i_max(n),n) = 2;
        u([i_min(n):ni 1:i_max(n)],n) = 1;
    else
        u([i_min(n):-1:1 ni:-1:i_max(n)],n) = 2;
        u(i_min(n):i_max(n),n) = 1;        
    end
end

Fn.u = griddedInterpolant(si,sj,u);

end