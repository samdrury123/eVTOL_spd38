function [b,h] = bl_fit_blade(xrrt,xr_hub,xr_cas,plot_stuff,param_sections)
% BL_FIT_BLADE  Fit a blade and reproduce with a spline definition
%
%   [b,h] = BL_FIT_BLADE(xrrt,xr_hub,xr_cas,plot_stuff,param_sections)
%   
%   xrrt - 3D float array of original blade polar coordinates
%   xr_hub - 2D float list of hub coordinates
%   xr_cas - 2D float list of casing coordinates
%   plot_stuff - 0 or 1 for showing working
%   param_sections - 0 or 1 for parameterising on sections
%   b - spline definition output
%   h - figure handle structure

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Default to parameterise the blade in one go
if exist('param_sections','var') == 0
    param_sections = 0;
end

% Define knots and spline order for each parameter
knots.sweep = [0.2 0.8]; order.sweep = 5;
knots.lean = [0.2 0.8]; order.lean = 5;
knots.tchord = [0.2 0.8]; order.tchord = 5;
knots.chi_le = [0.3 0.7]; order.chi_le = 4;
knots.chi_te = 0.5; order.chi_te = 4;
knots.dcam_le = 0.5; order.dcam_le = 4;
knots.dcam_te = 0.5; order.dcam_te = 4;
knots.qcam = 0.5; order.qcam = 4;
knots.thick = 0.5; order.thick = 4;
knots.le_te = 0.5; order.le_te = 4;

% Add on starting and ending knots
r_start = -0.025; r_end = 1.025;
varnames = fieldnames(knots);
for v = 1:length(varnames)
    varname = varnames{v};
    knots.(varname) = [r_start * ones(1,order.(varname))  knots.(varname) ...
        r_end * ones(1,order.(varname))];
end

% Number of sections
nj = 43; r_nondim_fit = linspace(-0.025,1.025,nj)'; r_nondim_fit(2) = 0; r_nondim_fit(end-1) = 1;

% Choose sections to fit
q = 2:nj-1;


%% Interpolate all sections onto lines of constant non-dimensional radius and extract parameters

% Extract meridional parameters
b = bl_parameterise_mer(xrrt,xr_hub,xr_cas);

% Calculate current non-dimensional radius
r_nondim = (xrrt(:,:,2) - interp1(xr_hub(:,1),xr_hub(:,2),xrrt(:,:,1),'pchip')) ./ ...
    (interp1(xr_cas(:,1),xr_cas(:,2),xrrt(:,:,1),'pchip') -...
    interp1(xr_hub(:,1),xr_hub(:,2),xrrt(:,:,1),'pchip'));

% Preallocate new coordinate array
xrrt_new = zeros(size(xrrt,1),nj,3);

% Loop around section and interpolate coordinates
for i = 1:size(xrrt,1)
    xrrt_new(i,:,1) = interp1(r_nondim(i,:),xrrt(i,:,1),r_nondim_fit,'pchip','extrap');
    xrrt_new(i,:,2) = interp1(r_nondim(i,:),xrrt(i,:,2),r_nondim_fit,'pchip','extrap');
    xrrt_new(i,:,3) = interp1(r_nondim(i,:),xrrt(i,:,3),r_nondim_fit,'pchip','extrap');
end

% Assign new coordinates
xrrt = xrrt_new(:,q,:);

% Parameterise blade geometry in one go or in sections
if param_sections == 0
    t = bl_parameterise_blade(xrrt,xr_hub,xr_cas,0,0);
    r_nondim_par = t.r_nondim;
else
    r_nondim_par = r_nondim_fit(q(1:2:end));
    t = parameterise_sections(xrrt(:,1:2:end,:),xr_hub,xr_cas,r_nondim_par);
end

% Re-use hub reference coordinate
b.x_ref = t.x_ref;


%% Basic radial fits to lean, sweep and chord

% Fit sweep and ensure zero value at hub
b.sweep = spap2(knots.sweep,order.sweep,t.r_nondim,t.sweep);
b.sweep.coefs = b.sweep.coefs - fnval(b.sweep,0);

% Fit lean and ensure zero value at hub
b.lean = spap2(knots.lean,order.lean,t.r_nondim,t.lean);
b.lean.coefs = b.lean.coefs - fnval(b.lean,0);

% Fit true chord
b.tchord = spap2(knots.tchord,order.tchord,t.r_nondim,t.tchord);


%% Fit camber and thickness distributions and smooth radially

% Size of parameters
[ni_par,nj_par] = size(t.cam);

% Preallocate variables
cam_fit = zeros(ni_par,nj_par); s_fit = zeros(ni_par,nj_par);
chi_le_temp = zeros(nj_par,1); chi_te_temp = zeros(nj_par,1);
dcam_le_temp = zeros(nj_par,1); dcam_te_temp = zeros(nj_par,1);
qcam_temp = zeros(nj_par,1);

% Fit camber distribution section by section
for j = 1:nj_par
    
    % Fit dimensional camber distribution with specified function
    cam_param = bl_fit_camber(t.s_cl(:,j),t.cam(:,j),t.chi_le(j),t.chi_te(j));
    
    % Store camber parameters
    chi_le_temp(j) = cam_param.chi_le; 
    chi_te_temp(j) = cam_param.chi_te; 
    dcam_le_temp(j) = cam_param.dcam_le; 
    dcam_te_temp(j) = cam_param.dcam_te; 
    qcam_temp(j) = cam_param.qcam; 
    cam_fit(:,j) = cam_param.cam;
    s_fit(:,j) = cam_param.s_cl;
end

% Fit camber parameters radially
b.chi_le = spap2(knots.chi_le,order.chi_le,t.r_nondim,chi_le_temp);
b.chi_te = spap2(knots.chi_te,order.chi_te,t.r_nondim,chi_te_temp);
b.dcam_le = spap2(knots.dcam_le,order.dcam_le,t.r_nondim,dcam_le_temp);
b.dcam_te = spap2(knots.dcam_te,order.dcam_te,t.r_nondim,dcam_te_temp);
b.qcam = spap2(knots.qcam,order.qcam,t.r_nondim,qcam_temp);

% Fit trailing edge thickness
b.thick_te = spap2(knots.le_te,order.le_te,t.r_nondim,t.thick_te);
spline_temp = fnxtr(b.thick_te,order.le_te);
thick_te_fit = fnval(spline_temp,r_nondim_par);

% Preallocate variables
thick_fit = zeros(ni_par,nj_par); s_thick_fit = zeros(ni_par,nj_par);
rad_le_temp = zeros(nj_par,1); s_thick_max_temp = zeros(nj_par,1);
rad_thick_max_temp = zeros(nj_par,1); wedge_te_temp = zeros(nj_par,1);

% Loop over all sections
for j = 1:nj_par
    
    % Fit thickness distribution with specified function
    thick_param = bl_fit_thickness(t.s_cl(:,j),t.S(:,j),thick_te_fit(j),0);

    % Record thickness parameters
    rad_le_temp(j) = thick_param.rad_le; 
    s_thick_max_temp(j) = thick_param.s_thick_max; 
    rad_thick_max_temp(j) = thick_param.rad_thick_max; 
    wedge_te_temp(j) = thick_param.wedge_te; 

    thick_fit(:,j) = thick_param.thick;
    s_thick_fit(:,j) = thick_param.s_cl;
end

% Fit thickness parameters radially
b.s_thick_max = spap2(knots.thick,order.thick,t.r_nondim,s_thick_max_temp);
b.rad_thick_max = spap2(knots.thick,order.thick,t.r_nondim,rad_thick_max_temp);
b.rad_le = spap2(knots.le_te,order.le_te,t.r_nondim,rad_le_temp);
b.wedge_te = spap2(knots.le_te,order.le_te,t.r_nondim,t.wedge_te);

% Fit max thickness
b.thick_max = spap2(knots.thick,order.thick,t.r_nondim,t.thick_max);


%% Assemble blade and back out parameters for comparison between intent and actual

% Construct new blade
b = bl_construct_blade(b,0,nj,[],param_sections);

% Reparameterise
if param_sections == 0
    t_rep = bl_parameterise_blade(b.xrrt(:,q,:),b.xr_hub,b.xr_cas,0,0);
else
%     t_rep = parameterise_sections(b.xrrt(:,q(1:2:end),:),b.xr_hub,b.xr_cas,r_nondim_par);
    t_rep = t;
end
% t_rep = bl_parameterise_blade(b.xrrt(:,q,:),b.xr_hub,b.xr_cas,0,0,flip_overide);

% Record reference geometry and raw parameters
b.ref = t; b.ref.xrrt = xrrt; b.ref.xr_hub = xr_hub; b.ref.xr_cas = xr_cas;
    

%% Plot comparison of blade fits

if plot_stuff == 1
    
    % Legend names
    leg = {'Old','Check','New'}; C = [1 0 0 ; 0 0 1 ; 0 0 0]; l = zeros(2,1);
    
    % Plot radial profiles of original blade, fitted number and reparameterised blade
    h.spline = figure();
    bl_spline_plot(t,h.spline,[1 0 0]);
    bl_spline_plot(t_rep,h.spline,[0 0 1]);
    bl_spline_plot(b,h.spline,[0 0 0]);
    legend(leg,'location','best')
    
    % Collect data on old and new geometries
    r_target = [0 0.25 0.5 0.75 1.0];
    xr = {b.ref.xr_hub b.ref.xr_cas ; b.xr_hub b.xr_cas};
    T{1} = t; T{2} = t_rep;
        
    % Open figure windows for distributions and sections    
    h.dist = figure('Position',[1 82 1280 895]); hold on;
    if param_sections == 0
        h.blade = figure('Position',[1 82 1280 895]); hold on; grid on; axis equal; view(3);
    end
    
    % Initialise axes
    figure(h.dist)
    for j = 1:length(r_target)
        h.ax(j,1) = subplot_tight(length(r_target),2,(j*2-1),0.05); hold on; grid on; box on;
        ylabel('Camber Line Angle');  axis([0 1 0 1]);
        h.ax(j,2) = subplot_tight(length(r_target),2,(j*2),0.05); hold on; grid on; box on;
        ylabel('Thickness Distribution'); axis([0 1 0 1]);
    end
    
    % Plot data on select sections for comparison
    for a = 1:length(T)
        for j = 1:length(r_target)

            % Get index of closest section
            [~,n] = min(abs(T{a}.r_nondim - r_target(j))); 

            % Plot camber distributions
            figure(h.dist)
            axes(h.ax(j,1)); plot(T{a}.s_cl(:,n),T{a}.cam(:,n),'color',C(a,:));
            
            % Plot thickness distributions
            axes(h.ax(j,2)); plot(T{a}.s_cl(:,n),T{a}.thick(:,n),'color',C(a,:));
            
            % Plot sections
            if param_sections == 0
                figure(h.blade);
                rt_1 = T{a}.xrrt_1(:,n,3) - T{a}.xrrt_1(:,n,2) * T{a}.t_ref;
                rt_2 = T{a}.xrrt_2(:,n,3) - T{a}.xrrt_2(:,n,2) * T{a}.t_ref;
                l(a) = plot3(T{a}.xrrt_1(:,n,1),rt_1,T{a}.xrrt_1(:,n,2),'-','color',C(a,:));
                plot3(T{a}.xrrt_2(:,n,1),rt_2,T{a}.xrrt_2(:,n,2),'-','color',C(a,:));

                % Plot camber lines
                rt = T{a}.xrrt_cam(:,n,3) - T{a}.xrrt_cam(:,n,2) * T{a}.t_ref;
                plot3(T{a}.xrrt_cam(:,n,1),rt,T{a}.xrrt_cam(:,n,2),'-','color',C(a,:));
            end
        end
        
        % Meridional lines
        if param_sections == 0
            plot3(xr{a,1}(:,1),zeros(size(xr{a,1},1),1),xr{a,1}(:,2),'-','color',C(a,:))
            plot3(xr{a,2}(:,1),zeros(size(xr{a,2},1),1),xr{a,2}(:,2),'-','color',C(a,:))
        end
        
        % Leading, trailing and centre lines
        if param_sections == 0
            rt_le = T{a}.xrrt_le(:,3) - T{a}.xrrt_le(:,2) * T{a}.t_ref;
            rt_te = T{a}.xrrt_te(:,3) - T{a}.xrrt_te(:,2) * T{a}.t_ref;
            rt_cen = T{a}.xrrt_cen(:,3) - T{a}.xrrt_cen(:,2) * T{a}.t_ref;
            plot3(T{a}.xrrt_le(:,1),rt_le,T{a}.xrrt_le(:,2),'-','color',C(a,:));
            plot3(T{a}.xrrt_te(:,1),rt_te,T{a}.xrrt_te(:,2),'-','color',C(a,:));
            plot3(T{a}.xrrt_cen(:,1),rt_cen,T{a}.xrrt_cen(:,2),'-','color',C(a,:));
        end
    end
    
    % Plot desired fitted geometry
    for j = 1:length(r_target)
        [~,n] = min(abs(r_nondim_par - r_target(j))); 
        axes(h.ax(j,1)); plot(s_fit(:,n),cam_fit(:,n),'k-');
        axes(h.ax(j,2)); plot(s_thick_fit(:,n),thick_fit(:,n),'k-');
    end
    
    % Add legends
    if param_sections == 0
        figure(h.blade); legend(l,leg(1:2),'location','best');
    end
    figure(h.dist); legend(leg,'location','best');
end

end

function t = parameterise_sections(xrrt,xr_hub,xr_cas,r_nondim)
% Call bl_parameterise_sections to build up the blade parameterisation

% Size of array to parameterise
nj = size(xrrt,2); t.r_nondim = r_nondim;

% Preallocate arrays
t.cam = zeros(200,nj); t.s_cl = zeros(200,nj);
t.chi_le = zeros(nj,1); t.chi_te = zeros(nj,1); t.dcam_le = zeros(nj,1); 
t.qcam = zeros(nj,1); t.dcam_te = zeros(nj,1); 
t.thick = zeros(200,nj); t.S = zeros(200,nj);
t.rad_le = zeros(nj,1); t.thick_max = zeros(nj,1); t.s_thick_max = zeros(nj,1); 
t.rad_thick_max = zeros(nj,1); t.wedge_te = zeros(nj,1); 
t.thick_te = zeros(nj,1); t.tchord = zeros(nj,1);
t.sweep = zeros(nj,1); t.lean = zeros(nj,1); 

% Loop over sections in the blade
for j = 1:nj
    
    % Parameterise current section
    t_sec = bl_parameterise_section(xrrt(:,j,[1 3]));
    
    % Assemble parameters into the array
    t.cam(:,j) = t_sec.cam; t.s_cl(:,j) = t_sec.s_cl;
    t.chi_le(j) = t_sec.chi_le; t.chi_te(j) = t_sec.chi_te; t.dcam_le(j) = t_sec.dcam_le;
    t.qcam(j) = t_sec.qcam; t.dcam_te(j) = t_sec.dcam_te;
    t.thick(:,j) = t_sec.thick; t.S(:,j) = t_sec.S;
    t.rad_le(j) = t_sec.rad_le; t.thick_max(j) = t_sec.thick_max; t.s_thick_max(j) = t_sec.s_thick_max; 
    t.rad_thick_max(j) = t_sec.rad_thick_max; t.wedge_te(j) = t_sec.wedge_te; 
    t.thick_te(j) = t_sec.thick_te; t.tchord(j) = t_sec.tchord;
    
    % Take reference coordinates from the hub
    if j == 1
        t.x_ref = t_sec.x_ref; r_ref = interp1(xr_hub(:,1),xr_hub(:,2),t.x_ref);
        t.t_ref = t_sec.xy_cen(2) / r_ref;
    end
    
    % Chordwise vector
    p = t_sec.xy_te - t_sec.xy_le; p = p / sum(p.^2)^0.5;
    
    % Reference line coordinates
    r_hub = interp1(xr_hub(:,1),xr_hub(:,2),t_sec.x_ref,'pchip','extrap');
    r_cas = interp1(xr_cas(:,1),xr_cas(:,2),t_sec.x_ref,'pchip','extrap');
    r = (r_cas - r_hub) * r_nondim(j) + r_hub;
    xrt_ref = [t.x_ref t.t_ref * r];
    
    % Resolve shift of centroid into chordiwse and normal direction
    dxrt = t_sec.xy_cen - xrt_ref;
    t.sweep(j) = sum(dxrt .* p);
    t.lean(j) = sum(dxrt .* [-p(2) p(1)]);
end


end

