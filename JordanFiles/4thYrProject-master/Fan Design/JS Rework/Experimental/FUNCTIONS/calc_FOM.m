function FOM = calc_FOM(T, P, area)
%% FOM = calc_FOM(T, P)
% Given thrust and power calculate figure of merit
% 
%   T     = Thrust double or vector or matrix
%   P     = Power double or vector or matrix with equal dimensions to T
%   area  = double axial flow area
%

rho = 1.225;
FOM = (T ./ P) .* sqrt(T ./ (2 * rho * area));
end