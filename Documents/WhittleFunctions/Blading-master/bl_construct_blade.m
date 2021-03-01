function b = bl_construct_blade(b,plot_stuff,nj,dev,gen_sections,smooth_te,tang_lean)
% BL_CONSTRUCT_BLADE  Create coordinates of a blade from a spline definition
%
%   b = BL_CONSTRUCT_BLADE(b,plot_stuff,nj,gen_sections,smooth_te,tang_lean)
%
%   b - data structure containing spline definition
%   plot_stuff - 0 or 1 for showing working
%   nj - optional argument to specify number of sections
%   dev - optional argument for thickness deviations
%   gen_sections - 0 or 1 to build on a camber surface or camber lines
%   smooth_te - 0 or 1 to smooth discontinuity at trailing edge circle
%   tang_lean - 0 or 1 to construct lean and sweep in blade or absolute frame

% Set default to not plot
if exist('plot_stuff','var') == 0 || isempty(plot_stuff) == 1
    plot_stuff = 0;
end

% Set default to 43 sections
if exist('nj','var') == 0 || isempty(nj) == 1
    nj = 43;
end

% Set default to no measured deviations
if exist('dev','var') == 0 || isempty(dev) == 1
    dev = [];
end

% Set default to build on a camber surface
if exist('gen_sections','var') == 0 || isempty(gen_sections) == 1
    gen_sections = 0;
end

% Set default to not smooth the trailing edge
if exist('smooth_te','var') == 0 || isempty(smooth_te) == 1
    smooth_te = 0;
end

% Set default to have lean and sweep in blade chord reference frame 
if exist('tang_lean','var') == 0 || isempty(tang_lean) == 1
    tang_lean = 0;
end

% Calculate non-dimensional radius to evaluate sections at
if nj == 21
    r_nondim = linspace(0,1,nj);
elseif nj == 47
    r_nondim = linspace(-0.075,1.075,nj);
else
    r_nondim = linspace(-0.025,1.025,nj);
end

% Extend deviations beyond hub and casing
if isempty(dev) == 0
    dev.r_nondim = [-0.025 dev.r_nondim 1.025];
    dev.thick = [dev.thick(:,1) dev.thick dev.thick(:,end)];
end

% Set number of streamwise coordinates
if nj == 21
    ni = 141;
else
    ni = 301;
end

% Set default to shape space leading edge
if isfield(b,'le') == 0
    b.le = 0;    
end


%% Evaluate camber and thickness distributions on all sections

% Preallocate thickness and camber arrays
xrrt_cam = zeros(ni,nj,3); thick = zeros(ni,nj); chi = zeros(ni,nj);

% Loop over all sections
for j = 1:nj
    
    % Evaluate splines at current section
    c = bl_spline_eval(b,r_nondim(j));
    
    % Calculate thickness distribution in shape space with two part cubic spline
    thick_param = bl_construct_thickness(c,ni,-smooth_te,0);
    s_cl = thick_param.s_cl; thick(:,j) = thick_param.thick * c.thick_max;

    % Add on deviations to thickness distribution
    if isempty(dev) == 0
        thick(:,j) = thick(:,j) + interp2(dev.r_nondim,dev.s,dev.thick,r_nondim(j),s_cl,...
            'spline') / thick_max;
    end
    
    % Evaluate camberline at same chordwise points as thickness distribution
    cam_param = bl_construct_camber(c,s_cl);
    xrt_cam = cam_param.xy_cam; xrt_chord = cam_param.xy_chord;
    chi(:,j) = cam_param.chi;
    
    % Calculate position of thickness distribution centroid
    s_cen = trapz(s_cl,thick(:,j) .* s_cl) ./ trapz(s_cl,thick(:,j));
    xrt_cen = [interp1(s_cl,xrt_cam(:,1),s_cen,'pchip') interp1(s_cl,xrt_cam(:,2),s_cen,'pchip')];
    
    % Calculate normalised chordwise vector and sweep and lean displacements
    if tang_lean == 0
        p = diff(xrt_chord([1 end],:),1,1); p = p / sum(p.^2).^0.5;
    else
        p = [1,0];
    end
    p = c.sweep * [p(1) p(2)] + c.lean * [-p(2) p(1)];
    
    % Shift camber line to desired location
    xrrt_cam(:,j,1) = xrt_cam(:,1) + p(1) - xrt_cen(1) + b.x_ref; 
    xrrt_cam(:,j,3) = xrt_cam(:,2) + p(2) - xrt_cen(2);
    
    % Set r-coordinate to non-dimensional radius until meridional lines are created
    xrrt_cam(:,j,2) = r_nondim(j);
end


%% Construct annulus lines from parameters and calculate camber surface radii

% Use camber surface for leading and trailing edge axial coordinates
b.xrrt = zeros(ni*2-1,nj,3); b.xrrt(1,:,:) = xrrt_cam(1,:,:); b.xrrt(ni,:,:) = xrrt_cam(end,:,:);

% Construct new annulus lines
b = bl_construct_mer(b);
    
% Hub and casing radii at camber surface axial coordinates
r_hub = interp1(b.xr_hub(:,1),b.xr_hub(:,2),xrrt_cam(:,:,1),'pchip');
r_cas = interp1(b.xr_cas(:,1),b.xr_cas(:,2),xrrt_cam(:,:,1),'pchip');

% Calculate radii of blade camber surface
xrrt_cam(:,:,2) = repmat(r_nondim,[ni 1]) .* (r_cas - r_hub) + r_hub;


%% Add thickness distribution onto camber to generate blade

% Check whether to build the blade from the camber surface or on sections
if gen_sections == 0

    % Calculate camber surface normals
    xyz_cam = bl_pol2cart(xrrt_cam);
    [nx,ny,nz] = surfnorm(xyz_cam(:,:,1),xyz_cam(:,:,2),xyz_cam(:,:,3));
    n_cam = cat(3,nx,ny,nz);

    % Construct both sides of blade from thickness and camber surface
    xrrt_1 = bl_cart2pol(xyz_cam + 0.5 * n_cam .* repmat(thick,[1 1 3]));
    xrrt_2 = bl_cart2pol(xyz_cam - 0.5 * n_cam .* repmat(thick,[1 1 3]));
    
else
    
    % Add on the thickness in polar coordinates on sections
    n_cam = cat(3,sind(chi),-cosd(chi));
    xrrt_1 = cat(3,xrrt_cam(:,:,1) + 0.5 * n_cam(:,:,1) .* thick,xrrt_cam(:,:,2),...
        xrrt_cam(:,:,3) + 0.5 * n_cam(:,:,2) .* thick);
    xrrt_2 = cat(3,xrrt_cam(:,:,1) - 0.5 * n_cam(:,:,1) .* thick,xrrt_cam(:,:,2),...
        xrrt_cam(:,:,3) - 0.5 * n_cam(:,:,2) .* thick);    
    
end

% Assemble both sides into one array
b.xrrt = [xrrt_1 ; flip(xrrt_2(1:end-1,:,:),1)];

% Record camber surface and non-dimensional radii
b.xrrt_cam = xrrt_cam; b.r_nondim = r_nondim;


%% Plot selected blade sections

% Check if plotting is turned on
if plot_stuff == 1
    
    % Open figure window
    figure(); hold on; grid on; box on; axis equal;
    
    % Choose section radii to plot
    r = [0 0.25 0.5 0.75 1];
    for j = 1:length(r)

        % Find nearest section to desired radius and plot in 3D
        [~,n] = min(abs(r_nondim - r(j)));
        plot3(b.xrrt(:,n,1),b.xrrt(:,n,3),b.xrrt(:,n,2),'k-')
    end
    
    % Plot hub and casing lines
    plot3(b.xr_hub(:,1),zeros(size(b.xr_hub,1),1),b.xr_hub(:,2),'k-');
    plot3(b.xr_cas(:,1),zeros(size(b.xr_cas,1),1),b.xr_cas(:,2),'k-');
    
end


end
