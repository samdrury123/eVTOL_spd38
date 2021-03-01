function xr = bl_interp_mer(xr)
% BL_INTERP_MER  Interpolate and increase resolution of meridional curves

% Choose resolution of line
d_tol = 2e-4;

% Calculate distance through curve
d = [0 ; cumsum(sum(diff(xr,1,1).^2,2).^0.5,1)];

% Calculate number of points required
ni = round(d(end) / d_tol);

% New points distribution
d_interp = linspace(0,d(end),ni)';

% Interpolate new meridional curve
xr = [interp1(d,xr(:,1),d_interp,'pchip') interp1(d,xr(:,2),d_interp,'pchip')];


end