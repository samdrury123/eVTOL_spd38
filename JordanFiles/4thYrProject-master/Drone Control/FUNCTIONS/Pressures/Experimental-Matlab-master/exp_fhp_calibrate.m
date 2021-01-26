function [c,h] = exp_fhp_calibrate(Iota,Tau,P_raw,probe,T,coeffs,plot_stuff)
% Calculate probe coefficients from measured pressures to get calibration
% map

% Default to Hodson & Dominy style coefficients
if exist('coeffs','var') == 0
    coeffs = 'hodson';
end

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Record probe angles
c.Iota = Iota; c.Tau = Tau;

% Extract free stream and centre hole data
%   Po = P_raw(:,:,:,1); P = P_raw(:,:,:,2); P_cn = P_raw(:,:,:,probe.N.P_cn+2); 

% Extract free stream and centre hole data without +2
  Po = P_raw(:,:,:,1); P = P_raw(:,:,:,2); P_cn = P_raw(:,:,:,probe.N.P_cn); 

% Calculate viscosity and density for probe Reynolds number
exp_air
mu = air.mu_ref*((air.T_ref+air.C)./(T+air.C)).*(T./air.T_ref).^(3/2);
ro = P ./ (air.R * T); 

% Get pressures from each probe hole
if strcmp(probe.arrange,'plus') == 1
    
      % "Plus" probe configuration
%       P_up = P_raw(:,:,:,probe.N.P_up+2); P_dn = P_raw(:,:,:,probe.N.P_dn+2); 
%       P_rt = P_raw(:,:,:,probe.N.P_rt+2); P_lf = P_raw(:,:,:,probe.N.P_lf+2);
%       P_av = 0.25 * (P_up + P_dn + P_rt + P_lf);
%     
      % "Plus" probe configuration without the +2?
        P_up = P_raw(:,:,:,probe.N.P_up); P_dn = P_raw(:,:,:,probe.N.P_dn); 
        P_rt = P_raw(:,:,:,probe.N.P_rt); P_lf = P_raw(:,:,:,probe.N.P_lf);
        P_av = 0.25 * (P_up + P_dn + P_rt + P_lf);
%     


elseif strcmp(probe.arrange,'cross') == 1
    
    % "Cross" probe configuration
    P_ur = P_raw(:,:,:,probe.N.P_ur+2); P_ul = P_raw(:,:,:,probe.N.P_ul+2); 
    P_dl = P_raw(:,:,:,probe.N.P_dl+2); P_dr = P_raw(:,:,:,probe.N.P_dr+2);
    P_av = 0.25 * (P_ur + P_ul + P_dl + P_dr);

end

% Calculate probe coefficients
if strcmp(probe.arrange,'plus') == 1 && strcmp(coeffs,'hodson') == 1
    
    % Coefficients - Hodson & Dominy - plus probe
    c.C_Po = (Po - P_cn) ./ (P_cn - P_av);
    c.C_P = (Po - P) ./ (P_cn - P_av);
    c.C_Alpha = (P_rt - P_lf) ./ (P_cn - P_av);
    c.C_Beta = (P_up - P_dn) ./ (P_cn - P_av);
    
elseif strcmp(probe.arrange,'cross') == 1 && strcmp(coeffs,'hodson') == 1
    
    % Coefficients for tilted probe - Hodson & Dominy adapted
    c.C_Po = (Po - P_cn) ./ (P_cn - P_av);
    c.C_P = (Po - P) ./ (P_cn - P_av);
    c.C_Beta = 0.5 * (P_ur + P_ul - P_dr - P_dl) ./ (P_cn - P_av);
    c.C_Alpha = 0.5 * (P_ur + P_dr - P_ul - P_dl) ./ (P_cn - P_av);
%     c.C_Beta = 0.5 * (P_ur + P_ul - 2*P_dl) ./ (P_cn - P_av);
%     c.C_Alpha = 0.5 * (2*P_ur - P_ul - P_dl) ./ (P_cn - P_av);
    
elseif strcmp(probe.arrange,'plus') == 1 && strcmp(coeffs,'curtis') == 1
    
    % Coefficients - Curtis & Pullan - plus probe
    P_sort = sort(cat(4,P_up,P_dn,P_rt,P_lf),4);
    P_ref = (2/3) * P_cn + (1/3) * mean(P_sort(:,:,:,[3 4]),4) - mean(P_sort(:,:,:,[1 2]),4);
    
    c.C_Po = (Po - P_cn) ./ P_ref;
    c.C_P = (Po - P) ./ P_ref;
    c.C_Alpha = (P_up - P_dn) ./ P_ref;
    c.C_Beta = (P_rt - P_lf) ./ (P_cn - P_av);
    
elseif strcmp(probe.arrange,'cross') == 1 && strcmp(coeffs,'curtis') == 1
    
    % Coefficients for tilted probe - Curtis & Pullan adapted
    P_sort = sort(cat(4,P_ur,P_ul,P_dl,P_dr),4);
    P_ref = (2/3) * P_cn + (1/3) * mean(P_sort(:,:,:,[3 4]),4) - mean(P_sort(:,:,:,[1 2]),4);
    
    c.C_Po = (Po - P_cn) ./ P_ref;
    c.C_P = (Po - P) ./ P_ref;
    c.C_Beta = 0.5 * (P_ur + P_ul - P_dr - P_dl) ./ P_ref;
    c.C_Alpha = 0.5 * (P_ur + P_dr - P_ul - P_dl) ./ P_ref;
    
end

% Compressible velocity calculation
M = ((2/(air.ga-1)) * ((P ./ Po) .^ ((1-air.ga)/air.ga) - 1)).^0.5;
To = (1 + (M.^2)*(air.ga-1)/2) .* T;
V = ((air.cp * To).^0.5) .* ((air.ga-1).^0.5) .* M .* ...
    ((1 + (M.^2)*(air.ga-1)/2).^(-0.5));

% Calculate probe reynolds number
c.Re = ro .* V * probe.D / mu;

% Calculate "mach number" coefficient
% c.C_M = P_av ./ P_cen;

% Record raw pressures
c.P_raw = P_raw;

% Plot calibration map
if plot_stuff == 1
    exp_fhp_plot(c)
end

end

