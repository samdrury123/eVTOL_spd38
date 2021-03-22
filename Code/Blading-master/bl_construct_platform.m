function b = bl_construct_platform(b,plat_type,tol)
% BL_CONSTRUCT_PLATFORM  Produce a platform mesh to machine from

% Specify default to hub
if exist('plat_type','var') == 0
    plat_type = 'hub';
end

% Specify default points resolution
if exist('tol','var') == 0
    tol = 0.4e-3;
end

% Default pitch if not specified
if isfield(b,'t_1') == 0
    b.t_1 = - pi / b.N;
    b.t_2 = pi / b.N;
end

% Extract limiting axial coordinates
if strcmp(plat_type,'hub') == 1 && isfield(b,'x_in_hub') == 1
    x_in = b.x_in_hub; x_out = b.x_out_hub; 
elseif strcmp(plat_type,'cas') == 1 && isfield(b,'x_in_cas') == 1
    x_in = b.x_in_cas; x_out = b.x_out_cas; 
else
    x_in = b.x_in; x_out = b.x_out; 
end

% Generate axial coordinates
ni = round((x_out - x_in)/(2*tol))*2+1; 
x = linspace(x_in,x_out,ni)';

% Interpolate radial coordinates
if strcmp(plat_type,'hub') == 1
    r = interp1(b.xr_hub(:,1),b.xr_hub(:,2),x,'pchip');
else
    r = interp1(b.xr_cas(:,1),b.xr_cas(:,2),x,'pchip');
end

% Generate tangential coordinates
nk = round(((b.t_2 - b.t_1) * mean(r))/(2*tol))*2+1; 
t = linspace(b.t_1,b.t_2,nk);

% Assemble coordinate matrices
x = repmat(x,[1 nk]); r = repmat(r,[1 nk]); rt = r .* repmat(t,[ni 1]);
b.xrrt_plat = cat(3,x,r,rt); 

% Store cartesian coordinates
b.xyz_plat = bl_pol2cart(b.xrrt_plat); 


end