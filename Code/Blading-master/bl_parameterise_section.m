function [t,h] = bl_parameterise_section(xy,fast_mode,plot_stuff,find_circ)
% BL_PARAMETERISE_SECTION  Determine thickness and camberline from a blade section
%
%   t = BL_PARAMETERISE_SECTION(xy,plot_stuff)
%
%   xy - 2D array of section coordinates [x1 y1 ; x2 y2 ; ... ; xn yn]
%   fast_mode - 0 or 1 for more accurate calculation
%   plot_stuff - 0 or 1 for showing working
%   t - data structure of output parameters

% Default to not show working
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Default to not use fast mode
if exist('fast_mode','var') == 0
    fast_mode = 0;
end

% Default to find trailing edge circle
if exist('find_circ','var') == 0
    find_circ = 1;
end

% Work in double precision
xy = double(xy);

% Define resolution of camber line
ni = 25;

% Number of iterations to use
n_outer = 10; n_inner = 5;

% Thickness step size
sf_thick = 0.01;

% Open figure window and axes
if plot_stuff == 1
    
    % Full scale window
    h.window = figure('position',[1 26 1280 895]);
    
    % Whole section plot
    h.xy = axes('position',[0.0500 0.3667 0.425 0.6]); hold on; grid on; box on;
    axis equal; xlabel('X'); ylabel('Y');
    
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
    
    % Colour to use for derived values
    col = lines(1);
    
else
    
    % Empty plot handle
    h = [];
end


%% Prepare geometry

% Specify the section from the approximate leading edge
[~,i_off] = min(xy(:,1)); xy = xy([i_off:end 2:i_off],:);
i_orig = 1:size(xy,1); i_orig = i_orig([i_off:end 2:i_off]);

% Flip section to point nose down
[~,i_min] = min(xy(:,1)); [~,i_max] = max(xy(:,1));
if xy(i_min,2) > xy(i_max,2)
    xy(:,2) = - xy(:,2); flip_sec = 1;
else
    flip_sec = 0;
end

% Ensure blade is specified the right way round
a = squeeze(diff(xy,1,1)); a = a ./ repmat(sum(a.^2,2).^0.5,[1 2]);
a = [a zeros(size(a,1),1)]; c = cross(a(2:end,:),a(1:end-1,:));
if sum(sign(c(:,3))) < 0
    xy = flip(xy,1); rev_sec = 1;
else
    rev_sec = 0;
end

% Plot the section line
if plot_stuff == 1
    axes(h.xy)
    plot(xy(:,1),xy(:,2),'color',col,'linewidth',2);
end

% Calculate normals to section
xyz = repmat(cat(3,xy(:,1),zeros(size(xy,1),1),xy(:,2)),[1 3 1]); xyz(:,2,2) = 1; xyz(:,3,2) = 2;
[nx,~,nz] = surfnorm(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3));
n_surf = [nx(:,1,1),nz(:,1,1)];


%% Initial camber line guess

% Calculate the approximate stagger
[~,i_max] = max(xy(:,1)); [~,i_min] = min(xy(:,1));
gam = atand(diff(xy([i_min i_max],2)) / diff(xy([i_min i_max],1)));

% Rotate the section
xy_rot = [xy(:,1) * cosd(-gam) - xy(:,2) * sind(-gam) xy(:,1) * sind(-gam) + xy(:,2) * cosd(-gam)];

% Split the section into two lines based upon max and minium axial coordinates
[x_max,i_max] = max(xy_rot(:,1)); [x_min,i_min] = min(xy_rot(:,1));
[xy_1,xy_2] = bl_split(xy_rot,i_min,i_max);

% Measure trailing edge radius
s_1 = (xy_1(:,1) - x_min) / (x_max - x_min);
s_2 = (xy_2(:,1) - x_min) / (x_max - x_min);
i_circ = bl_find_circ(xy_1,xy_2,s_1,s_2);
% i_circ = [length(xy_1) length(xy_2)];
rad_te = 0.5 * sum((xy_1(i_circ(1),:) - xy_2(i_circ(2),:)).^2)^0.5;
x_te = 0.5 * (xy_1(i_circ(1),1) + xy_2(i_circ(2),1));

% Interpolate hyperbolic spaced coordinates
xy_cam(:,1) = poly_bunch(x_min,x_max,2*rad_te,2*(x_max-x_te),ni);
xy_cam(:,2) = 0.5 * (interp1(xy_1(:,1),xy_1(:,2),xy_cam(:,1),'spline','extrap') +...
    interp1(xy_2(:,1),xy_2(:,2),xy_cam(:,1),'spline','extrap'));

% Rotate camber coordinates back by stagger angle
xy_cam = [xy_cam(:,1) * cosd(gam) - xy_cam(:,2) * sind(gam) ...
    xy_cam(:,1) * sind(gam) + xy_cam(:,2) * cosd(gam)];


%% Iterate on the camber line to find its true coordinates

% Calculate camber line parametric coordinates
sd_cam = cam2param(xy_cam);

% Decide whether to iterate on the camber line or not
if fast_mode == 0
    
    % Extract leading and trailing edge points
    xy_le = xy_cam(1,:); xy_te = xy_cam(end,:);
    
    % Plot the camber line guess
    xy_cam = param2cam(sd_cam,xy_le,xy_te);
    if plot_stuff == 1
        axes(h.xy); plot(xy_cam(:,1),xy_cam(:,2),'o--','color',col);
    end

    % Determine direction to iterate in based upon whether it is a turbine or compressor
    dd = sign(mean(sd_cam(:,2)));
    chi_le = atand(diff(xy_cam(1:2,2),1,1) ./ diff(xy_cam(1:2,1),1,1));
    chi_te = atand(diff(xy_cam(end-1:end,2),1,1) ./ diff(xy_cam(end-1:end,1),1,1));
    if abs(chi_te) > abs(chi_le)
        dd = - dd;
    end

    % Iterate on the camber line until thicknesses are equal on both sides
    for a = 1:n_outer
        
        % Move edges to be perpendicular to the section surface
        [xy_le,xy_te] = move_edges(xy_cam,xy);

        % Iterate on points between edges
        F = inf; n = 1;
        while max(abs(F)) > 1e-6 && n < n_inner
            
            % Calculate cartesian line from parametric
            xy_cam = param2cam(sd_cam,xy_le,xy_te);

            % Calculate thickness errors
            [F,thick] = fun_cl(xy_cam,xy,n_surf,0);

            % Update deviation by residual errors
            thick_max = max(min(thick,[],2)); F_max = max(abs(F));
            sd_cam(:,2) = sd_cam(:,2) + dd * F * min(0.5,sf_thick * thick_max / F_max);   
            
            % Fix leading and trailing edge locations
            sd_cam([1 end],2) = 0; F = F(2:end-1);
            
            % Plot convergence of camber line
            if plot_stuff == 1
                axes(h.xy); plot(xy_cam(:,1),xy_cam(:,2),':','color',col);
            end

            % Update inner iteration number
            n = n + 1;
        end

%         % Find new leading and trailing edge lines
%         xy_cam = refine_edges(xy_cam,xy);
%         
%         % Update normal parametric coordinate
%         d_le = sum(n_norm .* (xy_cam(1,:) - xy_le)); d_te = sum(n_norm .* (xy_cam(end,:) - xy_te));
%         sd_cam(:,2) = sd_cam(:,2) - interp1([0 1],[d_le d_te],sd_cam(:,1),'spline','extrap');
%         xy_le = xy_cam(1,:); xy_te = xy_cam(end,:);

        % Redistribute parametric coordinates towards uniform turning
        if a > 1
            
            % Calculate turning distribution along camberline
            n_cam = [diff(xy_cam(:,2)) -diff(xy_cam(:,1))];
            n_cam = n_cam ./ repmat(sum(n_cam.^2,2).^0.5,[1 2]);
            dchi = acosd(sum(n_cam(1:end-1,:) .* n_cam(2:end,:),2));
            chi = [0 ; cumsum(dchi)]; chi = [chi(1) ; 0.5 * (chi(1:end-1) + chi(2:end)) ; chi(end)];

            % Interpolate unform turning and uniform distance distributions
            ni = min(2*round(chi(end)),25);
            s_turn = interp1(chi,sd_cam(:,1),linspace(chi(1),chi(end),ni)');
            s_dist = dist_2d(xy_cam); tchord = s_dist(end);
            s_dist = poly_bunch(0,1,2*rad_te/tchord,2*rad_te/tchord,ni)';
            s_interp = smooth(0.5 * (s_turn + s_dist));
            
            % Reinterpolate parametric coordinates
            s_interp = smooth([sd_cam(1:2,1) ; s_interp(3:end-2) ; sd_cam(end-1:end,1)]);
            sd_cam = [s_interp interp1(sd_cam(:,1),sd_cam(:,2),s_interp,'spline')];
            xy_cam = param2cam(sd_cam,xy_le,xy_te);

        end
    end
else
    
    % Record leading and trailing edge points
    xy_le = xy_cam(1,:); xy_te = xy_cam(end,:);
    
end

% Extrapolate edges for a smoother result
% sd_cam([1 end],2) = fnval(fnxtr(spapi(optknt(sd_cam(2:end-1,1)',5),sd_cam(2:end-1,1)',...
%     sd_cam(2:end-1,2)'),5),[0 1]);
% % sd_cam([1 end],2) = interp1(sd_cam(2:end-1,1),sd_cam(2:end-1,2),[0 1],'spline','extrap');
% xy_cam = param2cam(sd_cam,xy_le,xy_te);
% xy_cam = refine_edges(xy_cam,xy);

% Plot final camberline before interpolation
if plot_stuff == 1
    axes(h.xy); plot(xy_cam(:,1),xy_cam(:,2),'x-','color',col,'markersize',15);
end

% Update size of line
ni = 200;

% Interpolate deviations at increased resolution
% s_interp = interp1(linspace(0,1,ni_temp)',sd_cam(:,1),linspace(0,1,ni)','pchip');
s_interp = hyperbolic_bunch(ni,0.001,0.001)';
d_interp = interp1(sd_cam(:,1),sd_cam(:,2),s_interp,'spline');
% d_interp = fnval(spapi(optknt(sd_cam(:,1)',5),sd_cam(:,1)',sd_cam(:,2)'),s_interp);

% Evaluate camber line at increased resolution
% for n = 1:3
xy_cam = param2cam([s_interp d_interp],xy_le,xy_te);
[~,thick_sides] = fun_cl(xy_cam,xy,n_surf,0);
% d_interp = d_interp + dd * 0.5 * F;
% end

% Plot optimised camber line
if plot_stuff == 1
    axes(h.xy)
    plot(xy_cam(:,1),xy_cam(:,2),'.-','color',col)  
end


%% Characterise the camber line and thickness distributions
    
% Calculate camber parameters
t_sec = bl_parameterise_camber(xy_cam);

% Record parameters
t.s_cl = t_sec.s_cl; t.cam = t_sec.cam;
t.chi_le = t_sec.chi_le; t.chi_te = t_sec.chi_te;
t.dcam_le = t_sec.dcam_le; t.dcam_te = t_sec.dcam_te;
t.qcam = t_sec.qcam;

% Plot camber parameters
if plot_stuff == 1

    % Plot non-dimensional camber line
    axes(h.cam); plot(t_sec.s_cl,t_sec.cam,'-','color',col)

    % Plot the camberline gradient fits
    s_fit = linspace(0,0.04,100); plot(s_fit,polyval(t_sec.p_le,s_fit),'k-')
    s_fit = linspace(0.96,1,100); plot(s_fit,polyval(t_sec.p_te,s_fit),'k-')

    % Plot differential of camber line
    axes(h.curv); [~, d2cds2] = grad_mg(t_sec.s_cl,t_sec.cam); 
    plot(t_sec.s_cl,d2cds2,'-','color',col)
end

% Average the thickness measured on both sides and enforce zero at leading and trailing edge
thick = sum(thick_sides,2);
thick([1 end]) = 0;

% Record error in thickness
t.thick_error = max(abs(diff(thick_sides,1,2)));

% Calculate thickness parameters
t_sec = bl_parameterise_thickness(t.s_cl,thick);

% Record parameters
t.S = t_sec.S; t.rad_le = t_sec.rad_le; t.thick_max = t_sec.thick_max; 
t.s_thick_max = t_sec.s_thick_max; t.rad_thick_max = t_sec.rad_thick_max;
t.wedge_te = t_sec.wedge_te; t.thick_te = t_sec.thick_te; 

% Scale thickness by maximum value
t.thick = thick / t.thick_max;

% Plot thickness parameters
if plot_stuff == 1

    % Plot the thickness distributions
    axes(h.thick)
    plot(t.s_cl,t.thick,'-','color',col)

    % Plot the trailing edge quadratic fit
    s_fit = linspace(0.96,1,100); plot(s_fit,polyval(t_sec.p_te,s_fit),'k-');

    % Plot leading edge circle fit
    th = 170:0.1:180; 
    plot(t_sec.p_le(3)*cosd(th) + t_sec.p_le(1),t_sec.p_le(3)*sind(th) + t_sec.p_le(2),'k-');

    % Plot fit at max thickness
    s_fit = linspace(t.s_thick_max - 0.05,t.s_thick_max + 0.05,100);
    plot(s_fit,polyval(t_sec.p_max,s_fit),'k-')

    % Plot shape space distribution
    axes(h.shape); plot(t.s_cl,t.S,'-','color',col)
    
    % Plot thickness on both sides
    axes(h.thick); plot(s_interp,2 * thick_sides / t.thick_max,'--');
    
end


%% Split the blade into two sides to record geometry

% Calculate camber line normals
xyz_cam = repmat(cat(3,xy_cam(:,1),zeros(ni,1),xy_cam(:,2)),[1 3 1]); 
xyz_cam(:,2,2) = 1; xyz_cam(:,3,2) = 2;
[nx,~,nz] = surfnorm(xyz_cam(:,:,1),xyz_cam(:,:,2),xyz_cam(:,:,3));
n_cam = [nx(:,1,1),nz(:,1,1)];

% Construct both sides of blade from thickness and camberline
xy_1 = xy_cam + 0.5 * n_cam .* repmat(thick,[1 2]);
xy_2 = xy_cam - 0.5 * n_cam .* repmat(thick,[1 2]);

% Split the blade into two sides of raw coordinates
[xy_1_raw,xy_2_raw,i_1,i_2] = split_lines(xy,xy_cam);

% Record the leading and trailing edge indices
i_le = i_1(1); i_te = i_1(end);

% Distance through surfaces
s_1 = [0 ; cumsum(sum(diff(xy_1,1,1).^2,2).^0.5)]; s_1 = s_1 / max(s_1);
s_2 = [0 ; cumsum(sum(diff(xy_2,1,1).^2,2).^0.5)]; s_2 = s_2 / max(s_2);

% Camberline distances on raw surfaces
s_1_raw = [0 ; cumsum(sum(diff(xy_1_raw,1,1).^2,2).^0.5)]; s_1_raw = s_1_raw / s_1_raw(end);
s_cl_1_raw = interp1(s_1,t.s_cl,s_1_raw,'pchip');
s_2_raw = [0 ; cumsum(sum(diff(xy_2_raw,1,1).^2,2).^0.5)]; s_2_raw = s_2_raw / s_2_raw(end);
s_cl_2_raw = interp1(s_2,t.s_cl,s_2_raw,'pchip');

% Find trailing edge circle indices
if find_circ == 1
    i_circ = bl_find_circ(xy_1_raw,xy_2_raw,s_cl_1_raw,s_cl_2_raw);
    i_circ(1) = i_1(i_circ(1)); i_circ(2) = i_2(i_circ(2)); 
else
    i_circ = [length(xy_1_raw) length(xy_2_raw)];
end


%% Characterise chord and stacking parameters

% Chord and true chord measured in polar plane
chord = sum(diff(xy_cam([1 end],:),1,1).^2,2).^0.5;
tchord = sum(sum(diff(xy_cam,1,1).^2,2).^0.5,1);

% Extract leading and trailing edge coordinates
xy_le = xy_cam(1,:); xy_te = xy_cam(end,:); 

% Calculate coordinates of thickness distribution centroid
xy_cen = zeros(1,2);
s_cen = trapz(t.s_cl,t.thick .* t.s_cl) ./ trapz(t.s_cl,t.thick);
for n = 1:2
    xy_cen(1,n) = interp1(t.s_cl,xy_cam(:,n),s_cen,'pchip');
end

% Axial reference coordinate
x_ref = xy_cen(1,1);


%% Output remaining parameters

% Calculated lines
t.xy_1 = xy_1; t.xy_2 = xy_2; t.xy_cam = xy_cam;
t.xy_le = xy_le; t.xy_te = xy_te;

% Raw lines
t.xy_1_raw = xy_1_raw; t.xy_2_raw = xy_2_raw; 
t.s_cl_1_raw = s_cl_1_raw; t.s_cl_2_raw = s_cl_2_raw;
t.i_1 = i_1; t.i_2 = i_2;
t.i_le = i_le; t.i_te = i_te;
t.i_circ = i_circ;

% Chord and stacking parameters
t.chord = chord; t.tchord = tchord;
t.x_ref = x_ref; t.xy_cen = xy_cen;

% Flip the tangential coordinates back
varnames = {'xy_1','xy_2','xy_cam','xy_le','xy_te','xy_1_raw','xy_2_raw','xy_cen'};
if flip_sec == 1
    for v = 1:length(varnames)
        t.(varnames{v})(:,2) = - t.(varnames{v})(:,2);
    end
end

% Flip the metal angles
if flip_sec == 1
    t.chi_le = - t.chi_le; t.chi_te = - t.chi_te;
end

% Reverse indices
varnames = {'i_1','i_2','i_le','i_te','i_circ'};
if rev_sec == 1
    for v = 1:length(varnames)
        t.(varnames{v}) = size(xy,1) - t.(varnames{v}) + 1;
    end
end

% Offset indices
for v = 1:length(varnames)
    t.(varnames{v}) = i_orig(t.(varnames{v}));
end


end

function [F,thick,xy_int_1,xy_int_2] = fun_cl(xy_cam,xy,n_surf,fast_mode)
% Calculate errors between thickness on each side

% Default to fast mode
if exist('fast_mode','var') == 0
    fast_mode = 1;
end

% Length of camber line
[ni,~] = size(xy_cam);

% Calculate camber line normals
xyz_cam = repmat(cat(3,xy_cam(:,1),zeros(ni,1),xy_cam(:,2)),[1 3 1]); 
xyz_cam(:,2,2) = 1; xyz_cam(:,3,2) = 2;
[nx,~,nz] = surfnorm(xyz_cam(:,:,1),xyz_cam(:,:,2),xyz_cam(:,:,3));
n_cam = [nx(:,1,1),nz(:,1,1)];
% n_cam = norm_2d(xy_cam);

% Split blade section into two sides
[xy_1,xy_2,i_1,i_2] = split_lines(xy,xy_cam);
n_1 = n_surf(i_1,:); n_2 = n_surf(i_2,:);

% Choose intersection point search method
if fast_mode == 1

    % Minimum distance between blade and camber lines
    [~,i_1] = min(pdist2(xy_cam,xy_1),[],2);
    [~,i_2] = min(pdist2(xy_cam,xy_2),[],2);

else

    % Most aligned point with camber line normal
    nsl_1 = size(xy_1,1); nsl_2 = size(xy_2,1);

    % First side
    p_1 = repmat(reshape(xy_1,[1 nsl_1 2]),[ni 1 1]) - repmat(reshape(xy_cam,[ni 1 2]),[1 nsl_1 1]);
    p_1 = p_1 ./ repmat(sum(p_1.^2,3).^0.5,[1 1 2]);
    [~,i_1] = max(abs(sum(p_1 .* repmat(reshape(n_cam,[ni 1 2]),[1 nsl_1 1]),3)),[],2);

    % Second side
    p_2 = repmat(reshape(xy_2,[1 nsl_2 2]),[ni 1 1]) - repmat(reshape(xy_cam,[ni 1 2]),[1 nsl_2 1]);
    p_2 = p_2 ./ repmat(sum(p_2.^2,3).^0.5,[1 1 2]);
    [~,i_2] = max(abs(sum(p_2 .* repmat(reshape(n_cam,[ni 1 2]),[1 nsl_2 1]),3)),[],2);

end

% Record closest points and normals
xy_int_1 = xy_1(i_1,:); xy_int_2 = xy_2(i_2,:); 
n_int_1 = n_1(i_1,:); n_int_2 = n_2(i_2,:); 

% Calculate intersection of camber line normals with blade line planes
thick_1 = abs(sum((xy_int_1 - xy_cam) .* n_int_1,2) ./ sum(n_cam .* n_int_1,2));
thick_2 = abs(sum((xy_int_2 - xy_cam) .* n_int_2,2) ./ sum(n_cam .* n_int_2,2));

% Calculate residual errors
F = thick_2 - thick_1;

% Output thicknesses
thick = [thick_1 thick_2];

end

function [xy_le,xy_te] = move_edges(xy_cam,xy)
% Update leading and trailing edge points so camber line intersects at right angles

% Extract edge coordinates
xy_edge = xy_cam([1 end],:);

% Extrapolate normals at leading edges
s_cam = dist_2d(xy_cam,1); s_interp = linspace(0,1,1000)';
xy_interp = [interp1(s_cam(2:end-1),xy_cam(2:end-1,1),s_interp,'spline','extrap') ...
    interp1(s_cam(2:end-1),xy_cam(2:end-1,2),s_interp,'spline','extrap')];
n_interp = norm_2d(xy_interp); n_edge = n_interp([1 end],:);

% Split coordinates into leading and trailing edge fifths
x = xy(:,1); x_max = max(x); x_min = min(x);
i_le = x < 0.8 * x_min + 0.2 * x_max;
i_te = x > 0.2 * x_min + 0.8 * x_max;

% Rearrange coordinates from trailing edge first
i_mid = round(length(x)/2); xy_loop = xy([i_mid:end 2:i_mid],:,:); 

% Extract split coordinates
if xy(1,1) < xy(i_mid,1)
    i_le = i_le([i_mid:end 2:i_mid]);
    xy_surf{1} = xy_loop(i_le,:); xy_surf{2} = xy(i_te,:);
else
    i_te = i_te([i_mid:end 2:i_mid]);
    xy_surf{1} = xy(i_le,:); xy_surf{2} = xy_loop(i_te,:);
end

% Loop over both edges
for m = 1:2
    
    % Calculate normals to section surface
    n = norm_2d(xy_surf{m});
    
    % Angle between edge normal and surface
    th = acosd(sum(n .* repmat(n_edge(m,:),[size(n,1) 1]),2));
    
    % Find the closest point to right angles
    [~,i_near] = min(abs(th - 90));
    
    % Interpolate coordinates near the closest point
    i = i_near-2:i_near+2;
    for n = 1:2
        xy_edge(m,n) = interp1(th(i),xy_surf{m}(i,n),90,'spline');
    end
end

% Return edge coordinates
xy_le = xy_edge(1,:); xy_te = xy_edge(2,:);
    
end


function xy_cam = refine_edges(xy_cam,xy)
% Update leading and trailing edge lines to lie on the section

% Refinement value for maximum alignment
ni_ref = 1000;

% Extract edge coordinates
xy_edge = xy_cam([1 end],:);

% Split coordinates into leading and trailing edge fifths
x = xy(:,1); x_max = max(x); x_min = min(x);
i_le = x < 0.8 * x_min + 0.2 * x_max;
i_te = x > 0.2 * x_min + 0.8 * x_max;

% Rearrange coordinates from trailing edge first
i_mid = round(length(x)/2); xy_loop = xy([i_mid:end 2:i_mid],:,:); 

% Extract split coordinates
if xy(1,1) < xy(i_mid,1)
    i_le = i_le([i_mid:end 2:i_mid]);
    xy_surf{1} = xy_loop(i_le,:); xy_surf{2} = xy(i_te,:);
else
    i_te = i_te([i_mid:end 2:i_mid]);
    xy_surf{1} = xy(i_le,:); xy_surf{2} = xy_loop(i_te,:);
end

% Loop over both edges
xy_near = zeros(ni_ref,2);
for m = 1:2

    % Calculate nearest points in line to edge coordinates
    [~,i_near] = min(pdist2(xy_edge(m,:),xy_surf{m}));

%     plot(xy_surf{m}(i_near,1),xy_surf{m}(i_near,2),'ko');
    
    % Interpolate more coordinates near the closest point
    i = i_near-2:i_near+2; s = linspace(0,1,length(i)); s_interp = linspace(0,1,ni_ref);
    for n = 1:2
        xy_near(:,n) = interp1(s,xy_surf{m}(i,n),s_interp,'pchip');
    end

    % Find the new closest point
    [~,i_near] = min(pdist2(xy_edge(m,:),xy_near));
    xy_edge(m,:) = xy_near(i_near,:);
end

% Record updated edges in camber line
xy_cam([1 end],:,:) = xy_edge;

end

function [xy_1,xy_2,i_1,i_2] = split_lines(xy,xy_cam)
% Split a blade line into triangulated lists of both sides at leading and trailing edge

% Length of surface array
ni = size(xy,1);

% Record list of indices
I = 1:ni;

% Split coordinates into leading and trailing edge fifths
x = xy(:,1); x_max = max(x); x_min = min(x);
q_le = x < 0.8 * x_min + 0.2 * x_max;
q_te = x > 0.2 * x_min + 0.8 * x_max;
xy_le = xy(q_le,:); xy_te = xy(q_te,:);
i_le = I(q_le); i_te = I(q_te); 

% Calculate nearest points to leading and trailing edges
[~,n_le] = min(pdist2(xy_cam(1,:),xy_le),[],2);
[~,n_te] = min(pdist2(xy_cam(end,:),xy_te),[],2);

% Get true indices of points at the leading and trailing edge
i_le = i_le(n_le); i_te = i_te(n_te); 

% Extract two halves of the section
[xy_1,xy_2,i_1,i_2] = bl_split(xy,i_le,i_te);

end

function xy_cam = param2cam(sd_cam,xy_le,xy_te)
% Calculate camber line cartesian coordinates from parametrics

% Size of camber line
ni = size(sd_cam,1);

% Calculate coordinates of points on the chord line
xy_chord = repmat(xy_le,[ni 1]) + repmat(xy_te - xy_le,[ni 1]) .* repmat(sd_cam(:,1),[1 2]);

% Calculate camber line from deviations from the chord line
n_norm = [xy_le(:,2)-xy_te(:,2) xy_te(:,1)-xy_le(:,1)]; n_norm = n_norm ./ sum(n_norm.^2).^0.5;
xy_cam = xy_chord + repmat(n_norm,[ni 1]) .* repmat(sd_cam(:,2),[1 2]);

end

function sd_cam = cam2param(xy_cam)
% Calculate camber line parametric coordinates from cartesians

% Size of camber line
ni = size(xy_cam,1);

% Extract leading and trailing edge coordinates
xy_le = xy_cam(1,:); xy_te = xy_cam(end,:); 

% Calculate chord and stagger angles
chord = sum((xy_te - xy_le).^2).^0.5; 
gam = atand((xy_te(:,2) - xy_le(:,2)) ./ (xy_te(:,1) - xy_le(:,1)));
n_chord = [cosd(gam) sind(gam)];

% Calculate camber line streamwise coordinates
s = dot((xy_cam - repmat(xy_le,[ni 1])),repmat(n_chord,[ni 1]),2) ./ repmat(chord,[ni 1]);

% Calclate camber line deviation from parametric coordinates
xy_chord = repmat(xy_le,[ni 1]) + chord * repmat(s,[1 2]) .* repmat(n_chord,[ni 1]);
d = dot(xy_cam - xy_chord,repmat([-n_chord(:,2) n_chord(:,1)],[ni 1]),2);

% Assemble parametric coordinates
sd_cam = [s d];

end

