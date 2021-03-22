function [A,d_th] = bl_calc_throat(b,plot_stuff)
% BL_CALC_THROAT  Compute the areas at inlet to blade row, throat and exit

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Construct blade geometry and add fillet if required
if isfield(b,'xrrt') == 0
    b = bl_construct_blade(b,0);
    b.xrrt = bl_construct_fillet(b.xrrt,b.fillet,b.xr_hub,b.xr_cas,0);
end

% Calculate cartesian coordinates of adjacent blades
xyz_1 = bl_pol2cart(b.xrrt); 
xrrt = cat(3,b.xrrt(:,:,[1 2]),(b.xrrt(:,:,3) ./ b.xrrt(:,:,2) + 2*pi / b.N) .* b.xrrt(:,:,2));
xyz_2 = bl_pol2cart(xrrt);
[ni,nj,~] = size(xyz_1);

% Plot blade meshes
if plot_stuff == 1
    figure(); hold on; grid on; box on; axis equal;
    mesh(xyz_1(:,:,1),xyz_1(:,:,3),xyz_1(:,:,2),zeros(ni,nj,3))
    mesh(xyz_2(:,:,1),xyz_2(:,:,3),xyz_2(:,:,2),zeros(ni,nj,3))
end

% Extract lines either side of the throat
d = sum((repmat(reshape(xyz_1,[ni 1 nj 3]),[1 ni 1 1]) - ...
    repmat(reshape(xyz_2,[1 ni nj 3]),[ni 1 1 1])).^2, 4).^0.5;
[~,i1] = min(min(d,[],2),[],1); [~,i2] = min(min(d,[],1),[],2);
i1 = reshape(i1,[1 nj]); i2 = reshape(i2,[1 nj]); 

% Extract coordinates of lines either side of the throat
xyz_th = zeros(nj,2,3);
for j = 1:nj
    xyz_th(j,:,:) = [xyz_1(i1(j),j,:) ; xyz_2(i2(j),j,:)];
end
xrrt_th = bl_cart2pol(xyz_th);
t = xrrt_th(:,:,3) ./ xrrt_th(:,:,2);

% Plot throat area
if plot_stuff == 1
    plot3(xyz_th(:,:,1)',xyz_th(:,:,3)',xyz_th(:,:,2)','r-')
end

% Calculate throat pitch
d_th = sum(diff(xyz_th,1,2).^2,3).^0.5;

% Calculate area in polar coordinates
dx_1 = xrrt_th(2:end,1:end-1,1) - xrrt_th(1:end-1,2:end,1);
dx_2 = xrrt_th(1:end-1,1:end-1,1) - xrrt_th(2:end,2:end,1);
dr_1 = xrrt_th(2:end,1:end-1,2) - xrrt_th(1:end-1,2:end,2);
dr_2 = xrrt_th(1:end-1,1:end-1,2) - xrrt_th(2:end,2:end,2);

t_av = 0.25 * (t(2:end,1:end-1) + t(1:end-1,1:end-1) + t(1:end-1,2:end) + t(2:end,2:end));

drt_1 = xrrt_th(2:end,1:end-1,2) .* (t(2:end,1:end-1) - t_av) - ...
    xrrt_th(1:end-1,2:end,2) .* (t(1:end-1,2:end) - t_av);
drt_2 = xrrt_th(1:end-1,1:end-1,2) .* (t(1:end-1,1:end-1) - t_av) - ...
    xrrt_th(2:end,2:end,2) .* (t(2:end,2:end) - t_av);

Ax = 0.5 * (dr_1.*drt_2 - dr_2.*drt_1);
Ar = 0.5 * (dx_2.*drt_1 - dx_1.*drt_2);
At = 0.5 * (dx_1.*dr_2 - dx_2.*dr_1);

% dA = (Ax.^2 + Ar.^2 + At.^2).^0.5;
dA = (Ax.^2 + At.^2).^0.5;

% Sum all areas
A = sum(dA);

end