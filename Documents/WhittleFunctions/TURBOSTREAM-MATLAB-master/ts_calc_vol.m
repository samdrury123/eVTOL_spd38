function [V,A_vec,A] = ts_calc_vol(x,r,rt)
% TS_CALC_VOL  Calculate cell volumes of a polar coordinate mesh
%
%   [V,A_vec,A] = TS_CALC_VOL(x,r,rt)
%
%   x - input 3D matrix of axial coordinates
%   r - input 3D matrix of radial coordinates
%   rt - input 3D matrix of radial times tangential coordinates
%   V - output array of cell volume
%   A_vec - output array of area vectors for each cell face
%   A - output array of area magnitudes for each cell face
%
%   Note:
%       Cell faces are in order -i, i, -j, +j, -k, +l
%       Area vectors point into centre of cell

% Assemble coordinate matrix
v = cat(4,x,r,rt); [ni,nj,nk,~] = size(v);

% Extract vertices for each cell
v_cell = cat(5,v(1:end-1,1:end-1,1:end-1,:),v(2:end,1:end-1,1:end-1,:),...
    v(1:end-1,2:end,1:end-1,:), v(1:end-1,1:end-1,2:end,:),...
    v(2:end,2:end,1:end-1,:), v(2:end,1:end-1,2:end,:),...
    v(1:end-1,2:end,2:end,:), v(2:end,2:end,2:end,:));

% Coordinates of cell centre
v_cen = mean(v_cell,5);

% Pitch shift tangential coordinates to zero
t = v_cen(:,:,:,3) ./ v_cen(:,:,:,2);
v_cell(:,:,:,3,:) = v_cell(:,:,:,3,:) - repmat(t,[1 1 1 1 8]) .* v_cell(:,:,:,2,:);

% Shift radial and axial coordinates to zero
v_cell(:,:,:,1:2,:) = v_cell(:,:,:,1:2,:) - repmat(v_cen(:,:,:,1:2),[1 1 1 1 8]);

% Extract vertices for each facet
v_face = cat(6,v_cell(:,:,:,:,[1 4 7 3]),v_cell(:,:,:,:,[2 5 8 6]),v_cell(:,:,:,:,[1 2 6 4]),...
    v_cell(:,:,:,:,[3 7 8 5]),v_cell(:,:,:,:,[1 3 5 2]),v_cell(:,:,:,:,[6 8 7 4]));

% Extract ridges on each facet
r_edge = diff(v_face(:,:,:,:,[1:end 1],:),1,5);
rx = cat(4,-r_edge(:,:,:,3,:,:,:),r_edge(:,:,:,2,:,:,:));
rr = cat(4,r_edge(:,:,:,3,:,:,:),-r_edge(:,:,:,1,:,:,:));
rt = cat(4,-r_edge(:,:,:,2,:,:,:),r_edge(:,:,:,1,:,:,:));

% Midpoints on each ridge 
v_mid = v_face + 0.5 * r_edge;

% Facet area in each direction by divergence theorem
A_vec = zeros(ni-1,nj-1,nk-1,3,6);
A_vec(:,:,:,1,:) = 0.5 * sum(sum(v_mid(:,:,:,2:3,:,:,:) .* rx,5),4);
A_vec(:,:,:,2,:) = 0.5 * sum(sum(v_mid(:,:,:,[1 3],:,:,:) .* rr,5),4);
A_vec(:,:,:,3,:) = 0.5 * sum(sum(v_mid(:,:,:,1:2,:,:,:) .* rt,5),4);

% Calculate total facet area
A = reshape(sum(A_vec.^2,4).^0.5,[ni-1 nj-1 nk-1 6]);

% Coordinates of face centres
v_mid = reshape(mean(v_face,5),[ni-1 nj-1 nk-1 3 6]);

% Calculate volume by divergence theorem
V = (1/3) * sum(sum(v_mid .* -A_vec,4),5);


end

