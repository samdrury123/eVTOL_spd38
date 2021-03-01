function b = bl_construct_tip(b,gap_type,tol)
% BL_CONSTRUCT_TIP  Produce blank geometry to machine from

% Specify default to tip
if exist('gap_type','var') == 0
    gap_type = 'cas';
end

% Specify default points resolution
if exist('tol','var') == 0
    tol = 0.1e-3;
end

% Offset meridional line by clearance gap
if strcmp(gap_type,'cas') == 1
    xr_gap = [b.xr_cas(:,1) b.xr_cas(:,2) - b.d_gap];
else
    xr_gap = [b.xr_hub(:,1) b.xr_hub(:,2) + b.d_gap];
end

% Interpolate blade coordinates tip section
sj = linspace(0,1,size(b.xrrt,2));
xrrt_sec = zeros(size(b.xrrt,1),3);
for i = 1:size(b.xrrt,1)
    r_sec = interp1(xr_gap(:,1),xr_gap(:,2),b.xrrt(i,:,1),'pchip');
    sj_gap = interp1(b.xrrt(i,:,2) - r_sec,sj,0,'pchip');
    for n = 1:3
        xrrt_sec(i,n) = interp1(sj,b.xrrt(i,:,n),sj_gap,'pchip');
    end
end

% Record indices of each side of the section
i_te = round(size(b.xrrt,1) / 2); i_1 = 1:i_te; i_2 = size(b.xrrt,1):-1:i_te;

% Get maximum thickness from both sides
thick = sum((xrrt_sec(i_1,:) - xrrt_sec(i_2,:)).^2,2).^0.5;

% Determine numbe of points to meet the required resolution
nk = round(max(thick) / tol);

% Generate linear scaling fraction
sk = linspace(0,1,nk);

% Scale both sides of the section and shape into an array
b.xrrt_tip = permute(repmat((xrrt_sec(i_2,:) - xrrt_sec(i_1,:)),[1 1 nk]) .* ...
    repmat(reshape(sk,[1 1 nk]),[i_te 3 1]) + repmat(xrrt_sec(i_1,:),[1 1 nk]),[1 3 2]);

% Store cartesian coordinates
b.xyz_tip = bl_pol2cart(b.xrrt_tip);


end

