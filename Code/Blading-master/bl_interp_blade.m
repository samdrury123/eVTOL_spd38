function [xrrt,xrrt_coarse] = bl_interp_blade(xrrt_orig,plot_stuff,s_clust,flip_overide,r_interp)
% BL_INTERP  Interpolate a poorly or coarsely defined blade to a high quality and resolution
%
%   [xrrt,xrrt_coarse] = BL_INTERP_BLADE(xrrt_orig,plot_stuff,s_clust,flip_overide,r_interp)
%
%   xrrt_orig - 3D or cell array of original blade definition
%   plot_stuff - 0 or 1 for showing working
%   s_clust - specification for concentrating point clustering
%   flip_overide - 0 or 1 to ensure turbine blades are interpolated correctly
%   xrrt - interpolated high-res polar blade coordinates
%   xrrt_coarse - interpolated sections at original radii
%
%   s_clust is 3D array
%       1-direction is list of all points
%       2-direction is [s_cen s_wid s_hei]
%       s_cen is non-dimensional coordinate to focus on
%       s_wid is width of clustering
%       s_hei is intensity of clustering

% Default input variables
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end
if exist('s_cluster','var') == 0
    s_clust = [];
end
if exist('flip_overide','var') == 0
    flip_overide = 0;
end

% Specify blade sizes
ni_edge = 100; ni_surf = 200; 
if exist('r_interp','var') == 0
    nj = 43; r_interp = linspace(-0.025,1.025,nj)'; 
else
    nj = length(r_interp);
end

% Record xrrt coordinates for blade depending on input type
if iscell(xrrt_orig) == 1
    xrrt_master = xrrt_orig;
elseif isfloat(xrrt_orig) == 1 && size(xrrt_orig,3) == 3
    xrrt_master = cell(size(xrrt_orig,2),1);
    for n = 1:size(xrrt_orig,2)
        xrrt_master{n} = squeeze(xrrt_orig(:,n,:));
    end
end
n_sec = length(xrrt_master);

% Ensure nose points down by checking rt-coordinates
j_mid = round(n_sec/2);
[~,i_le] = min(xrrt_master{j_mid}(:,1)); [~,i_te] = max(xrrt_master{j_mid}(:,1)); 
if xrrt_master{j_mid}(i_le,3) > xrrt_master{j_mid}(i_te,3)
    for n = 1:n_sec
        xrrt_master{n}(:,3) = - xrrt_master{n}(:,3);
    end
end

% Overide flipping for turbines
if flip_overide == 1
    for n = 1:n_sec
        xrrt_master{n}(:,3) = - xrrt_master{n}(:,3);
    end
end

% Interpolate sections repeatedly
ni = 2*ni_edge+2*ni_surf+1;
xrrt_coarse = zeros(ni,n_sec,3);
for n = 1:n_sec
    xrrt_coarse(:,n,:) = bl_interp_section(xrrt_master{n},0,s_clust,0,0);
end

% Plot original and interpolated blades
if plot_stuff == 1
    figure(); hold on; axis equal;
    for n = 1:n_sec
        plot3(xrrt_master{n}(:,1),xrrt_master{n}(:,3),xrrt_master{n}(:,2),'k-')
        plot3(xrrt_coarse(:,n,1),xrrt_coarse(:,n,3),xrrt_coarse(:,n,2),'r.')
    end
end

% Smooth spacings in the spanwise direction
s = [zeros(1,size(xrrt_coarse,2)) ; cumsum(sum(diff(xrrt_coarse(:,:,[1 3]),1,1).^2,3).^0.5)];
s = s ./ repmat(s(end,:),[size(s,1) 1]);
s_fit = zeros(size(s)); j_smooth = round(size(s,2)/4);
for i = 1:size(s,1)
    s_fit(i,:) = smooth(s(i,:),j_smooth);
end

% Reinterpolate sections
for j = 1:size(s,2)
    for o = 1:3
        xrrt_coarse(:,j,o) = interp1(s(:,j),xrrt_coarse(:,j,o),s_fit(:,j),'pchip');
    end
end

% Interpolate more sections in the r-direction and extrapolate beyond hub and casing lines
xrrt = zeros(ni,nj,3); q = r_interp >= 0 & r_interp <= 1;
for i = 1:ni
    r = (xrrt_coarse(i,:,2) - xrrt_coarse(i,1,2)) / (xrrt_coarse(i,end,2) - xrrt_coarse(i,1,2));
    xrrt(i,q,:) = interp1(reshape(r,[],1),reshape(xrrt_coarse(i,:,:),[],3),r_interp(q),'pchip');
    xrrt(i,q == 0,:) = interp1(reshape(r_interp(q),[],1),reshape(xrrt(i,q,:),[],3),...
        r_interp(q == 0),'pchip','extrap');
end

% Plot blade mesh
if plot_stuff == 1
    mesh(xrrt(:,:,1),xrrt(:,:,3),xrrt(:,:,2),zeros(ni,nj,3));
end

% Flip back if necessary
if flip_overide == 1
    xrrt(:,:,3) = - xrrt(:,:,3);
end

end