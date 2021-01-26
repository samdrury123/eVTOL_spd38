%%% Processing Fan Characteristic Data %%%
% 01/12/2020 - Sam Drury %
% This takes the RPM, Phi, Psi_ts, and Psi_tt data from the EDF, and
% outputs - Phi is varied through changing sigma, exit area ratio

clear; clc; close;
%% INPUTS

% Geometry
rhub = 20/1000; rcas = 60/1000; c = 20/1000;
rm = rms([rhub, rcas]);

% Vortex Exponents
A = -1;
B = -1.5;

%% DATA PREPROCESSING

% Adjusting values to mean line instead of casing values, according to
% vortex exponent
k = rcas/rm;
kphi = k;
kpsi = k^2;

%% DATA ENTRY

phi_design = 0.8; sigma = 1.1;
psi_design = 0.5*(phi_design/sigma)^2;

% Casing values - note psi is total-total
%%% Separate values according to stall, varied sigma, and design %%%
% Stall
% phi_st = [0.07 0.08 0.31 0.31 0.21 0.22 0.31 0.31]*kphi;
% psi_st = [0.037 0.042 0.073 0.074 0.058 0.059 0.081 0.085]*kpsi;
phi_st = [.335 .075 .217]*kphi;
psi_st = [.09 .044 .062]*kpsi;

% Varied sigma
% phi_var = [0.48 0.48 0.52 0.51 0.52 0.54 0.55 0.55 0.43 0.43 0.45 0.45]*kphi;
% psi_var = [0.118 0.118 0.119 0.112 0.106 0.109 0.115 0.115 0.136 0.137 ...
%     0.145 0.143]*kpsi;
phi_var = [.606 .56 .52 .445 .47 .65]*kphi;
psi_var = [.139 .135 .135 .147 .159 .128]*kpsi;

% Design
% phi = [0.52 0.54]*kphi;
% psi = [0.106 0.109]*kpsi;
phi = [.598]*kphi;
psi = [.137]*kpsi;

%% DATA PROCESSING

% Getting characteristic curve!
x = polyfit(phi_var, psi_var, 1);
x1 = linspace(min(phi_var), max(phi_var), length(phi_var));
y = polyval(x, x1);

% % Generating ideal curve
% phi_ideal = linspace(0.1,1,20);
% % 0.7 based on psi = 1 - A*phi, A = tan(beta2)+tan(alpha1)
% psi_ideal = 1 - 0.7.*phi_ideal;%0.5.*(phi_ideal./sigma).^2;

%% PLOTTING

figure(1); hold on; box on; grid on;
plot(phi_design, psi_design, 'mx', 'LineWidth', 3, 'MarkerSize', 10);
plot(phi, psi, 'go', 'LineWidth', 4, 'MarkerSize', 3);
plot(phi_var, psi_var, 'bo', 'LineWidth', 4, 'MarkerSize', 3);
plot(phi_st, psi_st, 'ro', 'LineWidth', 4, 'MarkerSize', 3);

legend('Theoretical Design Point', 'Experimental Design Point', 'Varying \sigma', 'Stall Region', 'Location', 'northwest');
%plot(x1,y, '-k', 'LineWidth', 1);
% plot(phi_ideal, psi_ideal, '-m', 'LineWidth', 3);
xlabel('\phi_m'); ylabel('\psi_m'); xlim([0,.9]); ylim([0,.3]);
set(gca,'FontSize',12)
% x0, y0, width, height
set(gcf, 'units','normalized','outerposition',[0.2 0.2 0.5 0.7]);

exportgraphics(gcf, 'Fan_chic.png', 'Resolution', 600);
