function i_circ = bl_find_circ(xy_1,xy_2,s_cl_1,s_cl_2)
% BL_FIND_CIRC  Find the indices defining the trailing edge circle

% Use surface distances if camberline distances are not specified
if exist('s_cl_1','var') == 0; s_cl_1 = dist_2d(xy_1,1); end;
if exist('s_cl_2','var') == 0; s_cl_2 = dist_2d(xy_2,1); end;

% Calculate curvature of both sides of the blade and limits
dydx_1 = abs(diff(grad_mg(xy_1(:,1),xy_1(:,2))));
lim_1 = interp1(s_cl_1(2:end),dydx_1,0.85);
dydx_2 = abs(diff(grad_mg(xy_2(:,1),xy_2(:,2))));
lim_2 = interp1(s_cl_2(2:end),dydx_2,0.85);

% Threshold curvature limits to find indices of the start of the circle
i_circ(1) = find(dydx_1 > 5*lim_1 & s_cl_1(2:end) > 0.85,1);
i_circ(2) = find(dydx_2 > 5*lim_2 & s_cl_2(2:end) > 0.85,1);


end