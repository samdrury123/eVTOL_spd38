function [c,h,e] = exp_chic_incidence(directory,chic_names,calib_name,N,r_mid,h,C,plot_stuff)
% EXP_CHIC_INCIDENCE  Calculate and plot characteristics against inlet incidence

% Default to plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end

% Load air properties
exp_air

% Check operation
if iscell(chic_names) == 1
    % Two separate characteristics, one casing, one fhp
    c_cas = exp_chic_statics(directory,chic_names{1},N,r_mid,[],[],0);
    load([directory chic_names{2}]);
else
    load([directory chic_names])
end

% Load in five hole probe incidence chic and calibration map
if isempty(calib_name) == 0
    load([directory calib_name]); c_fhp = c;
else
    c_fhp = e.c;
end

% Interpolate all data in time to desired rate and calculate rig operating
% point
e = exp_chic_interp(e);
e.N.P = N.P_in; e.N.Po = N.Po_in; e.N.Pa = N.Pa; e.N.T = N.To_in;
e.P = e.P_interp; e.T = e.T_interp;
e = exp_rig_calculate(e,r_mid);
c_inc.nrt = e.nrt; c_inc.mrtp = e.mrtp; c_inc.phi = e.phi;

% Calculate flow angles from five hole probe 
P_abs = e.P_interp - repmat(e.P_interp(:,N.Pa),[1 size(e.P_interp,2)]) + e.Pa;
P_fhp = P_abs(:,N.FHP); T_fhp = e.T(:,N.T_fhp);
d = exp_fhp_squares(c_fhp,P_fhp,e.probe,e.Iota(1),0,T_fhp);
c_inc.Alpha = d.Alpha; c_inc.t = e.t; c_inc.A = d.A;

% Cut incidence chic into different locations
ro = P_abs(:,N.P_in) ./ (air.R * e.T(:,N.To_in));
c_inc.psi = (P_abs(:,N.P_out) - P_abs(:,N.Po_in)) ./ (ro * e.U_mid^2);
c_inc.psi_fhp = (P_abs(:,N.P_out_fhp) - P_abs(:,N.P_in_fhp)) ./ (ro * e.U_mid^2);
c_inc = exp_chic_cut(c_inc);

% Determine operation type
if iscell(chic_names) == 1
    % Interpolate incidence onto pressure rise rates
    c_cas.Alpha = zeros(length(c_cas.phi),length(e.Iota));
    i = find(isnan(c_inc.Alpha)); i = [i ; length(c_inc.Alpha)+1];
    for n = 1:length(i)-1
        c_cas.Alpha(:,n) = interp1(c_inc.phi(i(n)+1:i(n+1)-1),c_inc.Alpha(i(n)+1:i(n+1)-1),c_cas.phi);
    end

    % Area average gas angles
    c_cas.Alpha = mean(c_cas.Alpha,2);
    c = c_cas;
else
    c = c_inc;
    
    % Interpolate incidence onto one line
    i = find(isnan(c.Alpha)); i = [i ; length(c.Alpha)+1];
    phi_av = c.phi(i(1)+1:i(2)-1);
    c.Alpha_interp = zeros(length(phi_av),length(e.Iota));
    c.psi_interp = zeros(length(phi_av),length(e.Iota));
    c.A_interp = zeros(length(phi_av),length(e.Iota));
    for n = 1:length(i)-1
        c.Alpha_interp(:,n) = interp1(c.phi(i(n)+1:i(n+1)-1),c.Alpha(i(n)+1:i(n+1)-1),phi_av);
        c.psi_interp(:,n) = interp1(c.phi(i(n)+1:i(n+1)-1),c.psi_fhp(i(n)+1:i(n+1)-1),phi_av);
        c.A_interp(:,n) = interp1(c.phi(i(n)+1:i(n+1)-1),c.A(i(n)+1:i(n+1)-1),phi_av);
    end
    c.Alpha_av = mean(c.Alpha_interp,2);
    c.psi_av = mean(c.psi_interp,2);
    c.phi_av = phi_av;
    c.A_av = mean(c.A_interp,2);
end
    
    
% Determine whether to plot or not
if plot_stuff == 1
    
    % Open a new figure window if handle is not given
    if exist('h','var') == 0 || isempty(h) == 1
        h = figure(); hold on; grid on; box on;
        xlabel('Incidence'); ylabel('Pressure Rise Coefficient');
    end

    % Plot results
    figure(h)
    plot(c.Alpha,c.psi,'-','Color',C,'LineWidth',2);
%     plot(c.Alpha,c.psi_rotor,'-','Color',C,'LineWidth',2);
%     plot(c.Alpha,c.psi_stator,'-','Color',C,'LineWidth',2);
end

end