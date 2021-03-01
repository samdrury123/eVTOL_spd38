% A script to demonstrate the use of blading tools

% Clear the workspace
clc
clear
close all


%% Define an example blade

% Blade parameters
r = [0 0.5 1]; 
b.sweep = [0 4e-4 4e-4 0e-4]; b.lean = [0 -5e-4 -5e-4 0]; b.tchord = 0.01 * ones(1,4); 
b.chi_le = [55 60 65 70]; b.chi_te = 20 * ones(1,4); b.dcam_le = -1.1 * ones(1,4); 
b.dcam_te = -0.4 * ones(1,4); b.qcam = 2 * ones(1,4); b.thick_te = 0.25 * ones(1,4); 
b.s_thick_max = 0.3 * ones(1,4); b.rad_thick_max = 0.2 * ones(1,4); b.rad_le = 7 * ones(1,4); 
b.wedge_te = 60 * ones(1,4); b.thick_max = 7e-4 * ones(1,4);

% Fit splines to numeric blade parameters
b = bl_spline_fit(b);

% Meridional parameters
b.A_le = 0.01; b.A_te = 0.0098; b.r_le = 0.13; b.r_te = 0.13;
b.drdx_hub_le = 0.03; b.drdx_hub_te = 0.018; b.drdx_cas_le = 0; b.drdx_cas_te = 0;
b.x_ref = 0.1;

% Blade count
b.N = 101;


%% Construct a blade from a spline definition

% Plot the spline definition
h = figure(); bl_spline_plot(b,h,[0 0 0]);

% Construct surface coordinates
b = bl_construct_blade(b,1);
gtg
% Construct longer meridional lines
b.x_in = 0.096; b.x_out = 0.1065; 
x = linspace(b.x_in,b.x_out,100)';
b.xr_hub = [x polyval(polyfit(b.xr_hub(:,1),b.xr_hub(:,2),3),x)];
b.xr_cas = [x polyval(polyfit(b.xr_cas(:,1),b.xr_cas(:,2),3),x)];

% Add a hub fillet
rad_fillet = 0.8e-3;
b.fillet.hub.s1 = [0 1]; b.fillet.hub.rad1 = [1 1] * rad_fillet;
b.fillet.hub.s2 = [0 1]; b.fillet.hub.rad2 = [1 1] * rad_fillet;
b.fillet.flag = [1 0]; b.fillet.tol = 1e-6;
xrrt = bl_construct_fillet(b.xrrt,b.fillet,b.xr_hub,b.xr_cas,1);

% Construct coordinates for a single section with closed trailing edge
c = bl_spline_eval(b,0.5); ote = 0;
c = bl_construct_section(c,1,ote);

% Fit an existing blade definition with splines
b = bl_fit_blade(b.xrrt,b.xr_hub,b.xr_cas,1);


%% Fit a RR blade definition file
%{
% Read RR file
b_rr = rr_read_jh05('example_directory/Blade_Definition');

% Get meridional lines
xr_hub = b_rr.xr{1}(b_rr.xr{1} > x_in & b_rr.xr{1} < x_out,:);
xr_cas = b_rr.xr{end}(b_rr.xr{end} > x_in & b_rr.xr{end} < x_out,:);

% Interpolate and increase resolutions
xrrt = bl_interp_blade(b_rr.xrrt);
xr_hub = bl_interp_mer(xr_hub); 
xr_cas = bl_interp_mer(xr_cas);

% Fit blade geometry
b = bl_fit_blade(xrrt,xr_hub,xr_cas,1);

% Blade counts
b.N = b_rr.N;

% Save spline definition
save('example_directory/spline.mat','b')
%}

