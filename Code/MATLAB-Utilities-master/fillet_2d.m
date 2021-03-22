function xy_fil = fillet_2d(xy_cen,r,xy_1,xy_2,th_tol,plot_stuff)
% FILLET_2D Find intersection of a fillet arc with two polylines

% Default angle tolerance
if exist('th_tol','var') == 0
    th_tol = 0.5;
end

% Default no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Open plotting window and plot lines to fillet
if plot_stuff == 1
    figure(); hold on; grid on; box on; axis equal;
    plot(xy_1(:,1),xy_1(:,2),'.-')
    plot(xy_2(:,1),xy_2(:,2),'.-')
    plot(xy_cen(1),xy_cen(2),'o')
end

% Calculate normals to each line
n_1 = norm_2d(xy_1); n_2 = norm_2d(xy_2);

% Roll balls along the first polyline and flip direction if required
xy_off_a = xy_1 + r * n_1; d_a = min(pdist2(xy_off_a,xy_cen)); 
xy_off_b = xy_1 - r * n_1; d_b = min(pdist2(xy_off_b,xy_cen)); 
if d_a < d_b; xy_off_1 = xy_off_a; n_sign_1 = 1; else; xy_off_1 = xy_off_b; n_sign_1 = -1; end;

% Roll balls along the second polyline and flip direction if required
xy_off_a = xy_2 + r * n_2; d_a = min(pdist2(xy_off_a,xy_cen)); 
xy_off_b = xy_2 - r * n_2; d_b = min(pdist2(xy_off_b,xy_cen)); 
if d_a < d_b; xy_off_2 = xy_off_a; n_sign_2 = 1; else; xy_off_2 = xy_off_b; n_sign_2 = -1; end;

% Calculate linear intersections between two offset polylines
[~,~,i,j] = intersections(xy_off_1(:,1),xy_off_1(:,2),xy_off_2(:,1),xy_off_2(:,2));
i = round(i); j = round(j);

% Plot offsets to polylines
if plot_stuff == 1
    plot(xy_off_1(:,1),xy_off_1(:,2),'-')
    plot(xy_off_2(:,1),xy_off_2(:,2),'-')
end

% Get indices in the vicinity of the intersection
i = i-2:i+2; i(i < 1) = []; i(i > size(xy_off_1,1)) = [];
j = j-2:j+2; j(j < 1) = []; j(j > size(xy_off_2,1)) = [];

% Increase resolution in the vicinity of the intersection
nn = 100;
s_1 = dist_2d(xy_1); s_int_1 = linspace(s_1(i(1)),s_1(i(end)),nn)';
xy_1 = [interp1(s_1,xy_1(:,1),s_int_1,'pchip') interp1(s_1,xy_1(:,2),s_int_1,'pchip')];
s_2 = dist_2d(xy_2); s_int_2 = linspace(s_2(j(1)),s_2(j(end)),nn)';
xy_2 = [interp1(s_2,xy_2(:,1),s_int_2,'pchip') interp1(s_2,xy_2(:,2),s_int_2,'pchip')];

% Recaulate normals and offsets
n_1 = n_sign_1 * norm_2d(xy_1); xy_off_1 = xy_1 + r * n_1;
n_2 = n_sign_2 * norm_2d(xy_2); xy_off_2 = xy_2 + r * n_2;

% Recalculate intersection
[x_cen,y_cen,i,j] = intersections(xy_off_1(:,1),xy_off_1(:,2),xy_off_2(:,1),xy_off_2(:,2));

% Plot centre intersection
if plot_stuff == 1
    plot(x_cen,y_cen,'x')
end

% Interpolate of intersecting normals
n_1 = [interp1(1:nn,n_1(:,1),i,'pchip') interp1(1:nn,n_1(:,2),i,'pchip')];
n_2 = [interp1(1:nn,n_2(:,1),j,'pchip') interp1(1:nn,n_2(:,2),j,'pchip')];

% Calculate angles of intersecting normals
th_1 = rad2deg(atan2(n_1(2),n_1(1))); th_2 = rad2deg(atan2(n_2(2),n_2(1))); 

% Choose only small arcs
if th_2 > th_1 && abs(th_1 - th_2) > 180
    th_2 = th_2 - 360;
elseif th_1 > th_2 && abs(th_1 - th_2) > 180
    th_1 = th_1 - 360;
end

% Define arc geometry
na = round(abs((th_1 - th_2) / th_tol));
th = linspace(th_1+90,th_2+90,na)';
xy_fil = [x_cen - r * sind(th) y_cen + r * cosd(th)];

% Plot fillet geometry
if plot_stuff == 1
    plot(xy_fil(:,1),xy_fil(:,2),'-')
end


end

