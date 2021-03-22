function t = bl_parameterise_blade(xrrt,xr_hub,xr_cas,fast_mode,plot_stuff,flip_overide)
% BL_PARAMETERISE_BLADE  Measure camber surface, thickness and other blade parameters in three dimensions
%
%   t = BL_PARAMETERISE_BLADE(xrrt,xr_hub,xr_cas,fast_mode,plot_stuff)
%
%   xrrt - 3D polar coordinate array specifying blade polar coordinates
%   xr_hub - 2D float list of hub coordinates
%   xr_cas - 2D float list of casing coordinates
%   fast_mode - 0 or 1 for more accurate calculation
%   plot_stuff - 0 or 1 for showing working
%   flip_overide - 0 or 1 for difficult blades
%   t - data structure of output parameters

% Default to not show working
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Default to most accurate mode
if exist('fast_mode','var') == 0
    fast_mode = 0;
end

% Default to not overide flipping
if exist('flip_overide','var') == 0
    flip_overide = 0;
end

% Work in double precision
xrrt = double(xrrt);

% Define resolution of camber surface
ni = 16; nj = 21;

% Open figure window and axes
if plot_stuff == 1
    
    % Full scale window
    h.window = figure('position',[1 26 1280 895]);
    
    % 3D blade plot
    h.xyz = axes('position',[0.0500 0.3667 0.425 0.6]); hold on; grid on; box on; view(3);
    axis equal; xlabel('X'); ylabel('Z'); zlabel('Y');
    
    % Camberline plot
    h.cam = axes('position',[0.05 0.05 0.425 0.2667]); hold on; grid on; box on; 
    xlabel('Chord'); ylabel('Camber'); axis([0 1 0 1]);
    
    % Curvature plot
    h.curv = subplot_tight(3,2,2,0.05); hold on; grid on; box on; 
    xlabel('Chord'); ylabel('Curvature'); axis([0 1 -5 15]);
    
    % Thickness plot
    h.thick = subplot_tight(3,2,4,0.05); hold on; grid on; box on; 
    xlabel('Chord'); ylabel('Thickness'); axis([0 1 0 1]);
    
    % Shape space plot
    h.shape = subplot_tight(3,2,6,0.05); hold on; grid on; box on; 
    xlabel('Chord'); ylabel('Shape'); axis([0 1 1 4]);
    
    % Colours of sections
    j_plot = round(linspace(1,nj,8)); cols = zeros(nj,3);
    cols(j_plot,:) = lines(length(j_plot));
    
else
    
    % Empty plot handle
    h = [];
end


%% Prepare geometry

% Flip blade to point nose down
[~,i_min] = min(xrrt(:,1,1)); [~,i_max] = max(xrrt(:,1,1));
if xrrt(i_min,1,3) > xrrt(i_max,1,3)
    xrrt(:,:,3) = - xrrt(:,:,3);
end

% Ensure blade is specified the right way round
a = squeeze(diff(xrrt(:,1,:),1)); c = cross(a(2:end,:),a(1:end-1,:));
if sum(sign(c(:,2))) < 0
    xrrt = flip(xrrt,1);
end

% Overide the flip
if flip_overide == 1
    xrrt = flip(xrrt,1);
end

% Extract a midspan section
% [~,j_mid] = min(abs(xrrt(1,:,2) - 0.5 * (max(xrrt(1,:,2)) + min(xrrt(1,:,2)))));
% xrrt_mid = squeeze(xrrt(:,j_mid,:));

% Shift to the origin
% t_mean = mean(xrrt_mid(:,3) ./ xrrt_mid(:,2)); xrrt(:,:,3) = xrrt(:,:,3) - t_mean * xrrt(:,:,2);

% Interpolate blade onto an increased number of sections
r_nondim = calc_r(xrrt,xr_hub,xr_cas); r_nondim_new = linspace(-0.05,1.05,89);
xrrt_new = zeros(size(xrrt,1),length(r_nondim_new),3);
for i = 1:size(xrrt,1)
    for n = 1:3
        xrrt_new(i,:,n) = interp1(r_nondim(i,:),xrrt(i,:,n),r_nondim_new,'pchip');
    end
end
xrrt = xrrt_new;

% Calculate cartesian coordinates
xyz = bl_pol2cart(xrrt);

% Mesh the blade surface
if plot_stuff == 1
    axes(h.xyz)
    surf(xyz(:,:,1),xyz(:,:,3),xyz(:,:,2),'facecolor','none','edgecolor',[0.5 0.5 0.5]);
end

% Calculate normals to surface
[nx,ny,nz] = surfnorm(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3)); n_surf = cat(3,nx,ny,nz);


%%%%%% CUT FILLETS



%% Initial camber surface guess

% Calculate radii of sections to split on
r_cas = interp1(xr_cas(:,1),xr_cas(:,2),xrrt(1,1,1),'pchip'); 
r_hub = interp1(xr_hub(:,1),xr_hub(:,2),xrrt(1,end,1),'pchip');
r_guess = linspace(r_hub,r_cas,nj);
[~,js] = min(abs(repmat(xrrt(1,:,2)',[1 nj]) - repmat(r_guess,[size(xrrt,2) 1])),[],1);

% Preallocate cambersurface coordinates
xyz_cam = zeros(ni,nj,3);

% Loop over desired sections
for j = 1:length(js)
    
    % Extract the current section
    xyz_sec = squeeze(xyz(:,js(j),:));
    
    % Calculate the approximate stagger of the current section
    [~,i_max] = max(xyz_sec(:,1)); [~,i_min] = min(xyz_sec(:,1));
    gam = atand(diff(xyz_sec([i_min i_max],3)) / diff(xyz_sec([i_min i_max],1)));
    
    % Rotate the section
    xyz_sec(:,[1 3]) = [xyz_sec(:,1) * cosd(-gam) - xyz_sec(:,3) * sind(-gam) ...
        xyz_sec(:,1) * sind(-gam) + xyz_sec(:,3) * cosd(-gam)];
    
    % Split the section into two lines based upon max and minium axial coordinates
    [x_max,i_max] = max(xyz_sec(:,1)); [x_min,i_min] = min(xyz_sec(:,1));
    [xyz_1,xyz_2] = bl_split(xyz_sec,i_min,i_max);
    
    % Interpolate regularly spaced coordinates
    x = hyperbolic_bunch(ni,0.02,0.02).' * (x_max - x_min) + x_min;
    for n = 2:3
        xyz_cam(:,j,n) = 0.5 * (interp1(xyz_1(:,1),xyz_1(:,n),x,'linear','extrap') +...
            interp1(xyz_2(:,1),xyz_2(:,n),x,'linear','extrap'));
    end
    xyz_cam(:,j,1) = x;
    
    % Rotate camber coordinates back by stagger angle
    xyz_cam(:,j,[1 3]) = [xyz_cam(:,j,1) * cosd(gam) - xyz_cam(:,j,3) * sind(gam) ...
        xyz_cam(:,j,1) * sind(gam) + xyz_cam(:,j,3) * cosd(gam)];
end  


%% Iterate on the camber surface to find its true coordinates

% Calculate camber surface parametric coordinates
srd_cam = cam2param(xyz_cam,xr_hub,xr_cas);

% Correct camber surface to lie in meridional bounds, and have plaid streamwise parametric coords
srd_cam(:,:,2) = repmat(linspace(0,1,nj),[ni 1]);
srd_cam(:,:,1) = repmat(mean(srd_cam(:,:,1),2),[1 nj]); srd_cam(1,:,1) = 0; srd_cam(end,:,1) = 1;
srd_cam([1 end],:,3) = 0;
    
% Decide whether to iterate on the camber surface or not
n_max = 20; n_iter = 1; thick_tol = (x_max - x_min) * 1.25e-4;
if fast_mode == 0
    
    % Extract leading and trailing edge lines
    xrrt_cam = bl_cart2pol(xyz_cam); xrrt_le = xrrt_cam(1,:,:); xrrt_te = xrrt_cam(end,:,:);
    
    % Plot the camber surface guess
    xyz_cam = param2cam(srd_cam,xrrt_le,xrrt_te,xr_hub,xr_cas);
    if plot_stuff == 1
        axes(h.xyz)
        for j = j_plot
            plot3(xyz_cam(:,j,1),xyz_cam(:,j,3),xyz_cam(:,j,2),'--','color',cols(j,:),'linewidth',2)
        end
    end
    
    % Determine direction to iterate in
    dd = sign(mean(mean(srd_cam(:,:,3))));

    % Iterate on the camber surface until thicknesses are equal on both sides
    for a = 1:4

        % Iterate on points between edges
        F = inf;
        while max(abs(F(:))) > thick_tol && n_iter < n_max

            % Calculate cartesian surface from parametric
            xyz_cam = param2cam(srd_cam,xrrt_le,xrrt_te,xr_hub,xr_cas);

            % Extrapolate leading and trailing edge coordinates
            for j = 1:nj
                for n = 1:3
                    xyz_cam([1 end],j,n) = interp1(srd_cam(2:end-1,j,1),xyz_cam(2:end-1,j,n),...
                        srd_cam([1 end],j,1),'spline','extrap');
                end
            end

            % Calculate thickness errors
            F = fun_cs(xyz_cam,xyz,n_surf);

            % Update deviation by residual errors
            srd_cam(:,:,3) = srd_cam(:,:,3) + dd * 0.5 * F;   
            n_iter = n_iter + 1;
            
        end
        
        % Move leading and trailing edge lines
        xyz_cam = refine_edges(xyz_cam,xyz,xr_hub,xr_cas);
        xrrt_le = bl_cart2pol(xyz_cam(1,:,:)); xrrt_te = bl_cart2pol(xyz_cam(end,:,:));
    end
end

% Update size of surface
ni_temp = ni; ni = 200;

% Interpolate deviations at increased resolution
s_interp = interp1(linspace(0,1,ni_temp)',srd_cam(:,1,1),linspace(0,1,ni)','pchip');
d_interp = zeros(ni,nj);
for j = 1:nj
    d_interp(2:end-1,j) = interp1(srd_cam(:,j,1),srd_cam(:,j,3),s_interp(2:end-1),'spline');
end

% Assembe parametric coordinates into arrays
s_interp = repmat(s_interp,[1 nj]);
r_interp = repmat(linspace(0,1,nj),[ni 1]);

% Evaluate camber surface at increased resolution
xyz_cam = param2cam(cat(3,s_interp,r_interp,d_interp),xrrt_le,xrrt_te,xr_hub,xr_cas);
[~,thick_sides] = fun_cs(xyz_cam,xyz,n_surf,0);

% Plot optimised camber surface
if plot_stuff == 1
    axes(h.xyz)
    for j = j_plot
        plot3(xyz_cam(:,j,1),xyz_cam(:,j,3),xyz_cam(:,j,2),'.-','color',cols(j,:))
    end   
end


%% Split the blade into two sides to record geometry

% Calculate camber surface normals
[nx,ny,nz] = surfnorm(xyz_cam(:,:,1),xyz_cam(:,:,2),xyz_cam(:,:,3));
n_cam = cat(3,nx,ny,nz);

% Average the thickness measured on both sides
thick = sum(thick_sides,3);

% Record error in thickness
t.thick_error = max(abs(diff(thick_sides,1,3)),[],1)';

% Construct both sides of blade from thickness and cambersurface
xrrt_1 = bl_cart2pol(xyz_cam + 0.5 * n_cam .* repmat(thick,[1 1 3]));
xrrt_2 = bl_cart2pol(xyz_cam - 0.5 * n_cam .* repmat(thick,[1 1 3]));

% Split the blade into two sides of raw coordinates
[~,~,~,~,i_1,i_2] = split_surfaces(xyz,xyz_cam);
xrrt_1_raw = cell(length(i_1),1); xrrt_2_raw = cell(length(i_1),1);
for j = 1:length(i_1)
    xrrt_1_raw{j} = bl_cart2pol(squeeze(xyz(i_1{j},j,:))); 
    xrrt_2_raw{j} = bl_cart2pol(squeeze(xyz(i_2{j},j,:)));
end

% Record the leading and trailing edge indices
i_le = zeros(length(i_1),1); i_te = zeros(length(i_1),1); 
for j = 1:length(i_1)
    i_le(j) = i_1{j}(1); i_te(j) = i_1{j}(end);
end

% Approximate raw camberline values by surface distances
s_cl_1_raw = cell(nj,1); s_cl_2_raw = cell(nj,1);
for j = 1:length(i_1)
    s_cl_1_raw{j} = [0 ; cumsum(sum(diff(xrrt_1_raw{j},1,1).^2,2).^0.5)]; 
    s_cl_1_raw{j} = s_cl_1_raw{j} / s_cl_1_raw{j}(end);
    s_cl_2_raw{j} = [0 ; cumsum(sum(diff(xrrt_2_raw{j},1,1).^2,2).^0.5)]; 
    s_cl_2_raw{j} = s_cl_2_raw{j} / s_cl_2_raw{j}(end);
end

% Find trailing edge circle indices
i_circ = zeros(length(i_1),2);
for j = 1:length(i_1)
    
    % Calculate curvature of both sides of the blade
    drtdx_1 = abs(diff(grad_mg(xrrt_1_raw{j}(:,1),xrrt_1_raw{j}(:,3))));
    lim_1 = interp1(s_cl_1_raw{j}(2:end),drtdx_1,0.9);
    drtdx_2 = abs(diff(grad_mg(xrrt_2_raw{j}(:,1),xrrt_2_raw{j}(:,3))));
    lim_2 = interp1(s_cl_2_raw{j}(2:end),drtdx_2,0.9);

    % Threshold curvature limits to find indices of the start of the circle
    i_circ(j,1) = i_1{j}(find(drtdx_1 > 5*lim_1 & s_cl_1_raw{j}(2:end) > 0.9,1));
    i_circ(j,2) = i_2{j}(find(drtdx_2 > 5*lim_2 & s_cl_2_raw{j}(2:end) > 0.9,1));
end


%% Characterise the camber surface and thickness distributions

% Calculate polar coordinates of camber surface
xrrt_cam = bl_cart2pol(xyz_cam);

% Pre-allocate camber arrays
t.s_cl = zeros(ni,nj); t.cam = zeros(ni,nj);
t.chi_le = zeros(nj,1); t.chi_te = zeros(nj,1);
t.dcam_le = zeros(nj,1); t.qcam = zeros(nj,1); t.dcam_te = zeros(nj,1);

% Loop over each slice of the camber surface
for j = 1:nj
    
    % Calculate parameters
    t_sec = bl_parameterise_camber(squeeze(xrrt_cam(:,j,[1 3])));
    
    % Record parameters
    t.s_cl(:,j) = t_sec.s_cl; t.cam(:,j) = t_sec.cam;
    t.chi_le(j) = t_sec.chi_le; t.chi_te(j) = t_sec.chi_te;
    t.dcam_le(j) = t_sec.dcam_le; t.dcam_te(j) = t_sec.dcam_te;
    t.qcam(j) = t_sec.qcam;
    
    % Plot camber parameters
    if plot_stuff == 1 && sum(j == j_plot) > 0
    
        % Plot non-dimensional camber line
        axes(h.cam); plot(t_sec.s_cl,t_sec.cam,'-','color',cols(j,:))
        
        % Plot the camberline gradient fits
        s_fit = linspace(0,0.04,100); plot(s_fit,polyval(t_sec.p_le,s_fit),'k-')
        s_fit = linspace(0.96,1,100); plot(s_fit,polyval(t_sec.p_te,s_fit),'k-')
        
        % Plot differential of camber line
        axes(h.curv); [~, d2cds2] = grad_mg(t_sec.s_cl,t_sec.cam); 
        plot(t_sec.s_cl,d2cds2,'-','color',cols(j,:))
    end
end

% Preallocate thickness arrays
t.thick = zeros(ni,nj); t.S = zeros(ni,nj);
t.rad_le = zeros(nj,1);
t.thick_max = zeros(nj,1); t.s_thick_max = zeros(nj,1); t.rad_thick_max = zeros(nj,1);
t.wedge_te = zeros(nj,1); t.thick_te = zeros(nj,1);

% Loop over each slice of the camber surface
for j = 1:nj

    % Calculate parameters
    t_sec = bl_parameterise_thickness(t.s_cl(:,j),thick(:,j));
    
    % Record parameters
    t.S(:,j) = t_sec.S; t.rad_le(j) = t_sec.rad_le; t.thick_max(j) = t_sec.thick_max; 
    t.s_thick_max(j) = t_sec.s_thick_max; t.rad_thick_max(j) = t_sec.rad_thick_max;
    t.wedge_te(j) = t_sec.wedge_te; t.thick_te(j) = t_sec.thick_te; 
    
    % Scale thickness by maximum value
    t.thick(:,j) = thick(:,j) / t.thick_max(j);
    
    % Plot thickness parameters
    if plot_stuff == 1 && sum(j == j_plot) > 0
        
        % Plot the thickness distributions
        axes(h.thick)
        plot(t.s_cl(:,j),t.thick(:,j),'-','color',cols(j,:))
    
        % Plot the trailing edge quadratic fit
        s_fit = linspace(0.96,1,100); plot(s_fit,polyval(t_sec.p_te,s_fit),'k-');
    
        % Plot leading edge circle fit
        th = 170:0.1:180; 
        plot(t_sec.p_le(3)*cosd(th) + t_sec.p_le(1),t_sec.p_le(3)*sind(th) + t_sec.p_le(2),'k-');
    
        % Plot fit at max thickness
        s_fit = linspace(t.s_thick_max(j) - 0.05,t.s_thick_max(j) + 0.05,100);
        plot(s_fit,polyval(t_sec.p_max,s_fit),'k-')
    
        % Plot shape space distribution
        axes(h.shape); plot(t.s_cl(:,j),t.S(:,j),'-','color',cols(j,:))
    end
end
    

%% Characterise stacking parameters

% Radial coordinates
r_nondim = reshape(srd_cam(1,:,2),[nj 1]);

% Chord and true chord measured in polar plane
chord = reshape(sum(diff(xrrt_cam([1 end],:,[1 3]),1,1).^2,3).^0.5,[nj 1]);
tchord = reshape(sum(sum(diff(xrrt_cam(:,:,[1 3]),1,1).^2,3).^0.5,1),[nj 1]);

% Extract leading and trailing edge coordinates
xrrt_le = squeeze(xrrt_cam(1,:,:)); xrrt_te = squeeze(xrrt_cam(end,:,:)); 

% Calculate coordinates of thickness distribution centroids
xrrt_cen = zeros(nj,3);
for j = 1:nj
    s_cen = trapz(t.s_cl(:,j),t.thick(:,j) .* t.s_cl(:,j)) ./ trapz(t.s_cl(:,j),t.thick(:,j));
    for n = 1:3
        xrrt_cen(j,n) = interp1(t.s_cl(:,j),xrrt_cam(:,j,n),s_cen,'pchip');
    end
end

% Reference location at hub
x_ref = xrrt_cen(1,1); t_ref = xrrt_cen(1,3) / xrrt_cen(1,2);

% Straight reference line coordinates
r_cas = interp1(xr_cas(:,1),xr_cas(:,2),x_ref,'pchip'); 
r_hub = interp1(xr_hub(:,1),xr_hub(:,2),x_ref,'pchip');
r = linspace(0,1,nj)' * (r_cas - r_hub) + r_hub; 
xrrt_ref = [ones(nj,1) * x_ref r r*t_ref];

% Chordwise vectors
p = xrrt_te - xrrt_le; p = p ./ repmat(sum(p.^2,2).^0.5,[1 3]);

% Resolve shift of centroids into chordwise and normal directions
dxrrt = (xrrt_cen(:,[1 3]) - xrrt_ref(:,[1 3]));
sweep = sum(dxrrt .* p(:,[1 3]),2);
lean = sum(dxrrt .* [-p(:,3) p(:,1)],2);


%% Output remaining parameters

% Calculated surfaces
t.xrrt_1 = xrrt_1; t.xrrt_2 = xrrt_2; t.xrrt_cam = xrrt_cam;
t.xrrt_le = xrrt_le; t.xrrt_te = xrrt_te;

% Raw surfaces
t.xrrt_1_raw = xrrt_1_raw; t.xrrt_2_raw = xrrt_2_raw; 
t.s_cl_1_raw = s_cl_1_raw; t.s_cl_2_raw = s_cl_2_raw;
t.i_1 = i_1; t.i_2 = i_2;
t.i_le = i_le; t.i_te = i_te;
t.i_circ = i_circ;

% Stacking parameters
t.r_nondim = r_nondim;
t.chord = chord; t.tchord = tchord;
t.sweep = sweep; t.lean = lean; t.x_ref = x_ref; t.t_ref = t_ref;
t.xrrt_cen = xrrt_cen;


end

function [F,thick,xyz_int_1,xyz_int_2] = fun_cs(xyz_cam,xyz,n_surf,fast_mode)
% Calculate errors between thickness on each side

% Default to fast mode
if exist('fast_mode','var') == 0
    fast_mode = 1;
end

% Size of camber surface
[ni,nj,~] = size(xyz_cam);

% Calculate camber surface normals
[nx,ny,nz] = surfnorm(xyz_cam(:,:,1),xyz_cam(:,:,2),xyz_cam(:,:,3));
n_cam = cat(3,nx,ny,nz);

% Split blade surface into two sides
[xyz_1,xyz_2,q_1,q_2] = split_surfaces(xyz,xyz_cam);
n_1 = reshape(n_surf(repmat(q_1,[1 1 3]) == 1),[],3); 
n_2 = reshape(n_surf(repmat(q_2,[1 1 3]) == 1),[],3);

% Calculate radial coordinates
r_cam = (xyz_cam(:,:,3).^2 + xyz_cam(:,:,2).^2).^0.5;
r_1 = (xyz_1(:,3).^2 + xyz_1(:,2).^2).^0.5;
r_2 = (xyz_2(:,3).^2 + xyz_2(:,2).^2).^0.5;

% Define bounding radii and axial coordinates
r_max_cam = [-inf max(r_cam,[],1) inf]; r_min_cam = [-inf min(r_cam,[],1) inf]; 
r_max_cam = ones(size(r_max_cam)) * 100; r_min_cam = zeros(size(r_min_cam));

% Loop over all camber surface sections
xyz_int_1 = zeros(ni,nj,3); xyz_int_2 = zeros(ni,nj,3);
n_int_1 = zeros(ni,nj,3); n_int_2 = zeros(ni,nj,3);
for j = 1:nj

    % Get coordinates of current radial slice
    q_1 = r_1 > r_min_cam(j) & r_1 < r_max_cam(j+2);
    q_2 = r_2 > r_min_cam(j) & r_2 < r_max_cam(j+2);
    xyz_sl_1 = xyz_1(q_1,:); n_sl_1 = n_1(q_1,:);
    xyz_sl_2 = xyz_2(q_2,:); n_sl_2 = n_2(q_2,:);
    
    % Choose intersection point search method
    if fast_mode == 1
        
        % Minimum distance between blade and camber surfaces
        [~,i_1] = min(pdist2(squeeze(xyz_cam(:,j,:)),xyz_sl_1),[],2);
        [~,i_2] = min(pdist2(squeeze(xyz_cam(:,j,:)),xyz_sl_2),[],2);
        
    else
        
        % Most aligned point with camber surface normal
        nsl_1 = size(xyz_sl_1,1); nsl_2 = size(xyz_sl_2,1);
        
        % First side
        p_1 = repmat(reshape(xyz_sl_1,[1 nsl_1 3]),[ni 1 1]) - repmat(xyz_cam(:,j,:),[1 nsl_1 1]);
        p_1 = p_1 ./ repmat(sum(p_1.^2,3).^0.5,[1 1 3]);
        [~,i_1] = max(abs(sum(p_1 .* repmat(n_cam(:,j,:),[1 nsl_1 1]),3)),[],2);
        
        % Second side
        p_2 = repmat(reshape(xyz_sl_2,[1 nsl_2 3]),[ni 1 1]) - repmat(xyz_cam(:,j,:),[1 nsl_2 1]);
        p_2 = p_2 ./ repmat(sum(p_2.^2,3).^0.5,[1 1 3]);
        [~,i_2] = max(abs(sum(p_2 .* repmat(n_cam(:,j,:),[1 nsl_2 1]),3)),[],2);
        
    end
        
    % Record closest points and normals
    xyz_int_1(:,j,:) = reshape(xyz_sl_1(i_1,:),[ni 1 3]); 
    xyz_int_2(:,j,:) = reshape(xyz_sl_2(i_2,:),[ni 1 3]); 
    n_int_1(:,j,:) = reshape(n_sl_1(i_1,:),[ni 1 3]); 
    n_int_2(:,j,:) = reshape(n_sl_2(i_2,:),[ni 1 3]); 
    
end 

% Calculate intersection of camber surface normals with blade surface planes
thick_1 = abs(sum((xyz_int_1 - xyz_cam) .* n_int_1,3) ./ sum(n_cam .* n_int_1,3));
thick_2 = abs(sum((xyz_int_2 - xyz_cam) .* n_int_2,3) ./ sum(n_cam .* n_int_2,3));

% Calculate residual errors
F = thick_2 - thick_1;

% Output thicknesses
thick = cat(3,thick_1,thick_2);

end

function xyz_cam = refine_edges(xyz_cam,xyz,xr_hub,xr_cas)
% Update leading and trailing edge lines to lie on the surface

% Refinement value for maximum alignment
ni_ref = 100;

% Number of sections
[~,nj,~] = size(xyz_cam);

% Extract edge coordinates
xyz_edge = xyz_cam([1 end],:,:);

% Split coordinates into leading and trailing edge fifths
x = xyz(:,1,1); x_max = max(x); x_min = min(x);
i_le = x < 0.8 * x_min + 0.2 * x_max;
i_te = x > 0.2 * x_min + 0.8 * x_max;

% Extract split coordinates
i_mid = round(length(x)/2); xyz_loop = xyz([i_mid:end 2:i_mid],:,:); i_le = i_le([i_mid:end 2:i_mid]);
xyz_surf{1} = xyz_loop(i_le,:,:); xyz_surf{2} = xyz(i_te,:,:);

% Non-dimensional radial coordinates
r_edge = calc_r(bl_cart2pol(xyz_edge),xr_hub,xr_cas);

% Loop over both edges
for m = 1:2
    
    % Interpolate coordinates onto the same non-dimensional radii
    r_surf = calc_r(bl_cart2pol(xyz_surf{m}),xr_hub,xr_cas);
    xyz_sec = zeros(size(xyz_surf{m},1),nj,3);
    for i = 1:size(xyz_surf{m})
        for n = 1:3
            xyz_sec(i,:,n) = interp1(r_surf(i,:),xyz_surf{m}(i,:,n),r_edge(m,:),'pchip');
        end
    end
    
    % Preallocate nearby coordinates
    xyz_near = zeros(ni_ref,3);
    
    % Loop over all sections and refine position of the edge
    for j = 1:nj
    
        % Calculate nearest points in surface to edge coordinates
        xyz_p = reshape(xyz_edge(m,j,:),[1 3]);
        [~,i_near] = min(pdist2(xyz_p,squeeze(xyz_sec(:,j,:))));
    
        % Interpolate more coordinates near the closest point
        i = i_near-2:i_near+2; s = linspace(0,1,length(i)); s_interp = linspace(0,1,ni_ref);
        for n = 1:3
            xyz_near(:,n) = interp1(s,xyz_sec(i,j,n),s_interp,'pchip');
        end
        
        % Find the new closest point
        [~,i_near] = min(pdist2(xyz_p,xyz_near));
        xyz_edge(m,j,:) = xyz_near(i_near,:);
    end
end

% Record updated edges in camber surface
xyz_cam([1 end],:,:) = xyz_edge;

end

function [xyz_1,xyz_2,q_1,q_2,i_1,i_2] = split_surfaces(xyz,xyz_cam)
% Split a blade surface into triangulated lists of both sides at leading and trailing edge

% Record size and reshape blade coordinates into a list
[ni,nj,~] = size(xyz); xyz_list = reshape(xyz,[ni*nj 3]);

% Create lists of indices
I = reshape(repmat(reshape(1:ni,[ni 1]),[1 nj]),[ni*nj 1]); 
J = reshape(repmat(reshape(1:nj,[1 nj]),[ni 1]),[ni*nj 1]); 

% Split list into leading and trailing edge fifths
x = xyz_list(:,1); x_max = max(x); x_min = min(x);
q_le = x < 0.8 * x_min + 0.2 * x_max;
q_te = x > 0.2 * x_min + 0.8 * x_max;
xyz_le = xyz_list(q_le,:); xyz_te = xyz_list(q_te,:);
i_le = I(q_le); j_le = J(q_le); i_te = I(q_te); j_te = J(q_te);

% Calculate nearest points to leading and trailing edges
[~,n_le] = min(pdist2(squeeze(xyz_cam(1,:,:)),xyz_le),[],2);
[~,n_te] = min(pdist2(squeeze(xyz_cam(end,:,:)),xyz_te),[],2);

% Get true indices of points at the leading and trailing edge
i_le = i_le(n_le); j_le = j_le(n_le); i_te = i_te(n_te); j_te = j_te(n_te);

% Interpolate streamwise indices at every location
q = i_le < ni/2; i_le(q) = i_le(q) + ni;
i_le = round(interp1(j_le,i_le,(1:nj)','linear','extrap'));
i_te = round(interp1(j_te,i_te,(1:nj)','linear','extrap'));
q = i_le > ni; i_le(q) = i_le(q) - ni;

% Split blade into two halves
q_1 = zeros(ni,nj); q_2 = zeros(ni,nj); i_1 = cell(nj,1); i_2 = cell(nj,1);
for j = 1:nj
    if i_te(j) > i_le(j)
        i_1{j} = i_le(j):i_te(j); i_2{j} = [i_le(j):-1:1 ni-1:-1:i_te(j)];
    else
        i_1{j} = [i_le(j):ni-1 1:i_te(j)]; i_2{j} = i_le(j):-1:i_te(j);
    end
    q_1(i_1{j},j) = 1; q_2(i_2{j},j) = 1;
end

% Extract lists of both blade surface coordinates
xyz_1 = xyz_list(reshape(q_1 == 1,[ni*nj 1]),:); 
xyz_2 = xyz_list(reshape(q_2 == 1,[ni*nj 1]),:);

end

function xyz_cam = param2cam(srd_cam,xrrt_le,xrrt_te,xr_hub,xr_cas)
% Calculate camber surface cartesian coordinates from parametrics

% Size of camber surface
ni = size(srd_cam,1);

% Calculate non-dimensional radius of leading and trailing edge lines
r_le = calc_r(xrrt_le,xr_hub,xr_cas); r_te = calc_r(xrrt_te,xr_hub,xr_cas); 

% Interpolate xrt coordinates
xrt_le = cat(3,interp1(r_le,xrrt_le(:,:,1),srd_cam(1,:,2),'pchip'),...
    interp1(r_le,xrrt_le(:,:,3),srd_cam(1,:,2),'pchip'));
xrt_te = cat(3,interp1(r_te,xrrt_te(:,:,1),srd_cam(end,:,2),'pchip'),...
    interp1(r_te,xrrt_te(:,:,3),srd_cam(end,:,2),'pchip'));

% Calculate coordinates of points on the chord line
xrt_chord = repmat(xrt_le,[ni 1]) + repmat(xrt_te - xrt_le,[ni 1]) .* repmat(srd_cam(:,:,1),[1 1 2]);

% Calculate deviations from the chord line
n_norm = cat(3,xrt_le(:,:,2)-xrt_te(:,:,2),xrt_te(:,:,1)-xrt_le(:,:,1)); 
n_norm = n_norm ./ repmat(sum(n_norm.^2,3).^0.5,[1 1 2]);
xrt_cam = xrt_chord + repmat(n_norm,[ni 1 1]) .* repmat(srd_cam(:,:,3),[1 1 2]);

% Calculate radius from meridional lines
r_hub = interp1(xr_hub(:,1),xr_hub(:,2),xrt_cam(:,:,1),'pchip');
r_cas = interp1(xr_cas(:,1),xr_cas(:,2),xrt_cam(:,:,1),'pchip');
r_cam = (r_cas - r_hub) .* srd_cam(:,:,2) + r_hub;

% Calculate cartesian coordinates
xyz_cam = bl_pol2cart(cat(3,xrt_cam(:,:,1),r_cam,xrt_cam(:,:,2)));

end

function srd_cam = cam2param(xyz_cam,xr_hub,xr_cas)
% Calculate camber surface parametric coordinates from cartesians

% Polar coordinates of camber surface
xrrt_cam = bl_cart2pol(xyz_cam);

% Size of camber surface
ni = size(xrrt_cam,1);

% Extract leading and trailing edge coordinates
xrrt_le = xrrt_cam(1,:,:); xrrt_te = xrrt_cam(end,:,:); 

% Calculate chord and stagger angles
chord = sum((xrrt_te(:,:,[1 3]) - xrrt_le(:,:,[1 3])).^2,3).^0.5; 
gam = atand((xrrt_te(:,:,3) - xrrt_le(:,:,3)) ./ (xrrt_te(:,:,1) - xrrt_le(:,:,1)));
n_chord = cat(3,cosd(gam),sind(gam));

% Calculate camber surface radial coordinates
r = calc_r(xrrt_cam,xr_hub,xr_cas);

% Calculate camber surface streamwise coordinates
s = dot((xrrt_cam(:,:,[1 3]) - repmat(xrrt_le(:,:,[1 3]),[ni 1 1])),repmat(n_chord,[ni 1 1]),3) ...
    ./ repmat(chord,[ni 1]);

% Calclate camber surface deviation from parametric coordinates
xrt_chord = repmat(xrrt_le(:,:,[1 3]),[ni 1]) + repmat(chord,[ni 1 2]) .* repmat(s,[1 1 2]) .*...
    repmat(n_chord,[ni 1 1]);
d = dot(xrrt_cam(:,:,[1 3]) - xrt_chord,repmat(cat(3,-n_chord(:,:,2),n_chord(:,:,1)),[ni 1 1]),3);

% Assemble parametric coordinates
srd_cam = cat(3,s,r,d);

end

function r = calc_r(xrrt,xr_hub,xr_cas)

% Calculate radius of hub and casing lines at desired axial coordinates
r_hub = interp1(xr_hub(:,1),xr_hub(:,2),xrrt(:,:,1),'pchip');
r_cas = interp1(xr_cas(:,1),xr_cas(:,2),xrrt(:,:,1),'pchip');

% Non-dimensional radial coordinate
r = (xrrt(:,:,2) - r_hub) ./ (r_cas - r_hub);

end