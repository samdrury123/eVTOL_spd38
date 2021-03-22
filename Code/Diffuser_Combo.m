%   Diffuser_Combo
%                               
%   Generate an example fan stator and diffuser combo
%
%   Written by James Taylor              
%   May 2020

% Clear workspace and load eVTOL environment
clear; close all; clc;
eVTOL

% Load stator example
load([dr.geom 'NRF_3b_splines.mat']); b = B{2};

% Modify inlet angle
b.chi_le.coefs = b.chi_le.coefs + 8;
b.sweep.coefs = b.sweep.coefs * 1.5;
b.lean.coefs = b.lean.coefs * 2;

% Modify hub and casing lines
b.x_in = 0.1; b.x_out = 0.35;
b.A_te = 1.5;
b.drdx_hub_le = -0.6;
b.drdx_hub_te = -0.7;
b.drdx_cas_le = 0.6;
b.drdx_cas_te = 0.6;

% Plot the blade
b = bl_construct_blade(b,1);
m = bl_construct_mer({b}); b.xr_hub = m.xr_hub; b.xr_cas = m.xr_cas;
plot3(m.xr_hub(:,1),zeros(size(m.xr_hub,1),1),m.xr_hub(:,2),'k-')
plot3(m.xr_cas(:,1),zeros(size(m.xr_cas,1),1),m.xr_cas(:,2),'k-')

% Generate splitter blades
c = bl_spline_eval(b,0.5);
t = bl_construct_thickness(c);
js = [18 30];
for j = 1:length(js)
    xr_cam = squeeze(b.xrrt_cam(:,js(j),[1 2]));
    s_cam = dist_2d(xr_cam,1); n_cam = norm_2d(xr_cam,0);
    thick = repmat(interp1(s_cam,t.s_cl,t.thick,'pchip') * c.thick_max * 0.5,[1 2]);
    xr_splitter = [xr_cam + n_cam .* thick ; ...
        flip(xr_cam(1:end-1,:) - n_cam(1:end-1,:) .* thick(1:end-1,:),1)];
    plot3(xr_splitter(:,1),zeros(size(xr_splitter,1),1),xr_splitter(:,2),'k-');
    b.xr_splitter{j} = xr_splitter;
end

% Export the blades for CAD
bl_write_ibl([dr.cad 'diffuser_combo.ibl'],b);

