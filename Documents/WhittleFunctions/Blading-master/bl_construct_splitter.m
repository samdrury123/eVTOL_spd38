function b = bl_construct_splitter(b,plot_stuff,nj)
% BL_CONSTRUCT_SPLITTER  Create coordinates of splitter blades from spline definitions
%
%   b = BL_CONSTRUCT_SPLITTER(b,plot_stuff,nj)
%
%   b - data structure containing spline definition in cell array called split
%   plot_stuff - 0 or 1 for showing working
%   nj - optional argument to specify number of sections

% Set default to not plot
if exist('plot_stuff','var') == 0 || isempty(plot_stuff) == 1
    plot_stuff = 0;
end

% Set default to 45 sections
if exist('nj','var') == 0 || isempty(nj) == 1
    nj = 45;
end

% Calculate non-dimensional pitch to evaluate sections at
if nj == 21
    t_nondim = linspace(0,1,nj);
else
    t_nondim = linspace(-0.05,1.05,nj);
end

% Set number of streamwise coordinates
if nj == 21
    ni = 141;
else
    ni = 301;
end
% ni = 201;


%% Evaluate camber surface and thickness distributions, then add together for splitter geometry

% Loop over all splitter definitions
for m = 1:length(b.split)
    
    % Section parameters at the join with the blade
    b.split{m}.tchord = 1; c = bl_spline_eval(b.split{m},0); 
    
    % Set display options and tight tolerancing
    if plot_stuff == 1
        options = optimoptions('fsolve','Display','iter','TolFun',1e-10); 
    else
        options = optimoptions('fsolve','Display','off','TolFun',1e-10);
    end
    
    % Set initial guess for true chord
    q = fnval(b.tchord,b.split{m}.r_nondim);
    
    % Blade leading edge and trailing edge coordinates from camber surface
    xr_bl_le = squeeze(b.xrrt_cam(1,:,:)); xr_bl_te = squeeze(b.xrrt_cam(end,:,:)); 
    
    % Splitter leading edge coordinates
    x_sec = interp1(b.r_nondim',b.xrrt_cam(:,:,1)',b.split{m}.r_nondim,'pchip');
    r_sec = interp1(b.r_nondim',b.xrrt_cam(:,:,2)',b.split{m}.r_nondim,'pchip');
    sx = (x_sec - x_sec(1)) / (x_sec(end) - x_sec(1));
    xr_le = interp1(sx',[x_sec' r_sec'],c.x_nondim(1),'pchip');
    
    % Optimise to find the correct value of true chord
    f = @(q) fun_chord(q,c,xr_le,b.xr_hub,b.xr_cas,xr_bl_le,xr_bl_te); 
    q = fsolve(f,q,options);
    
    % Preallocate thickness and camber arrays
    xrrt_cam = zeros(ni,nj,3); thick = zeros(ni,nj);   
   
    % Loop over all sections
    for j = 1:nj
        
        % Evaluate splines at current section
        c = bl_spline_eval(b.split{m},t_nondim(j)); c.tchord = c.fchord * q;
        
        % Calculate thickness distribution in shape space with two part cubic spline
        thick_param = bl_construct_thickness(c,ni,0,0);
        s_cl = thick_param.s_cl; thick(:,j) = thick_param.thick * c.thick_max;
        
        % Evaluate camberline at same chordwise points as thickness distribution
        t_cam = bl_construct_camber(c,s_cl); xr_cam = t_cam.xy_cam;
        
        % Calculate position of thickness distribution centroid
        s_cen = trapz(s_cl,thick(:,j) .* s_cl) ./ trapz(s_cl,thick(:,j));
        xr_cen = [interp1(s_cl,xr_cam(:,1),s_cen,'pchip') interp1(s_cl,xr_cam(:,2),s_cen,'pchip')];
        
        % Shift camberline with lean and sweep
        xrrt_cam(:,j,1) = xr_cam(:,1) + c.sweep - xr_cen(1);
        xrrt_cam(:,j,2) = xr_cam(:,2) + c.lean - xr_cen(2);
        
    end
    
    % Align the camber surface at blade join to specified non-dimensional radial and axial coordianates
    j = find(t_nondim == 0,1);
    xrrt_cam(:,:,1:2) = xrrt_cam(:,:,1:2) + ...
        repmat(reshape(xr_le,[1 1 2]) - xrrt_cam(1,j,1:2),[ni nj 1]);

    % Hub and casing radii at joining axial coordinates
    r_hub = interp1(b.xr_hub(:,1),b.xr_hub(:,2),xrrt_cam(:,j,1),'pchip');
    r_cas = interp1(b.xr_cas(:,1),b.xr_cas(:,2),xrrt_cam(:,j,1),'pchip');
    
    % Non-dimensional radius at join
    r_nondim = (xrrt_cam(:,j,2) - r_hub) ./ (r_cas - r_hub);
    
    % Blade camber surface coordinates at join
    x_bl = interp1(b.r_nondim',b.xrrt_cam(:,:,1)',r_nondim,'pchip');
    rt_bl = interp1(b.r_nondim',b.xrrt_cam(:,:,3)',r_nondim,'pchip');
    for i = 1:ni; xrrt_cam(i,j,3) = interp1(x_bl(i,:),rt_bl(i,:),xrrt_cam(i,j,1),'pchip'); end;
    
    % Interpolate tangential offset for all axial coordinates
    dt = interp1(xrrt_cam(:,j,1),xrrt_cam(:,j,3) ./ xrrt_cam(:,j,2),xrrt_cam(:,:,1),'linear','extrap');
    
    % Tangential coordinates from offset values
    xrrt_cam(:,:,3) = (dt + repmat(t_nondim,[ni 1]) * 2 * pi / b.N) .* xrrt_cam(:,:,2);
    
    % Calculate camber surface normals
    xyz_cam = bl_pol2cart(xrrt_cam);
    [nx,ny,nz] = surfnorm(xyz_cam(:,:,1),xyz_cam(:,:,2),xyz_cam(:,:,3));
    n_cam = cat(3,nx,ny,nz);

    % Construct both sides of blade by adding thickness to camber surface
    xrrt_1 = bl_cart2pol(xyz_cam + 0.5 * n_cam .* repmat(thick,[1 1 3]));
    xrrt_2 = bl_cart2pol(xyz_cam - 0.5 * n_cam .* repmat(thick,[1 1 3]));

    % Assemble both sides into one array
    b.split{m}.xrrt = [xrrt_1 ; flip(xrrt_2(1:end-1,:,:),1)];

    % Record camber surface and non-dimensional pitch
    b.split{m}.xrrt_cam = xrrt_cam; b.split{m}.t_nondim = t_nondim;
    
end


%% Plot selected splitter sections

% Check if plotting is turned on
if plot_stuff == 1
    
    % Open figure window
    figure(); hold on; grid on; box on; axis equal;
    
    % Plot the blades
    xrrt = b.xrrt_cam; xyz = bl_pol2cart(xrrt); C = 0.5 * ones(size(xrrt));
    mesh(xyz(:,:,1),xyz(:,:,3),xyz(:,:,2),C);
    xrrt(:,:,3) = xrrt(:,:,3) + xrrt(:,:,2) * 2 * pi / b.N; xyz = bl_pol2cart(xrrt);
    mesh(xyz(:,:,1),xyz(:,:,3),xyz(:,:,2),C);
    
    % Choose sections to plot
    t_cam = [0 0.25 0.5 0.75 1]; cols = lines(length(b.split));
    for j = 1:length(t_cam)

        % Find nearest section to desired radius and plot in 3D
        [~,n] = min(abs(t_nondim - t_cam(j)));
        for m = 1:length(b.split)
            xyz = bl_pol2cart(b.split{m}.xrrt);
            plot3(xyz(:,n,1),xyz(:,n,3),xyz(:,n,2),'-','color',cols(m,:))
            xyz = bl_pol2cart(b.split{m}.xrrt_cam);
            mesh(xyz(:,:,1),xyz(:,:,3),xyz(:,:,2),repmat(reshape(cols(m,:),[1 1 3]),[ni nj 1]));
        end
    end
    
end


end


function F = fun_chord(q,c,xr_le,xr_hub,xr_cas,xr_bl_le,xr_bl_te);
% Calculate error in trailing edge coordinates for given splitter chord

% Set true chord
c.tchord = q;

% Evaluate camberline
t_cam = bl_construct_camber(c);

% Calculate coordinates centred on the leading edge
xr_cam = t_cam.xy_cam + repmat(xr_le,[length(t_cam.s) 1]);

% Calcualate non-dimensional coordinates of trailing edge
[~,sx] = bl_passage_coords([xr_cam zeros(length(t_cam.s),1)],xr_hub,xr_cas,xr_bl_le,xr_bl_te);

% Error in axial coordinates at trailing edge
F = sx(t_cam.i_te) - c.x_nondim(2);


end


