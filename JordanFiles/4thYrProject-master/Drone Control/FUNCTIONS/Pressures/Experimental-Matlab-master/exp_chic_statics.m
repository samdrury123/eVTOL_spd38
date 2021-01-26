function [c,h,e,l] = exp_chic_statics(directory,chic_name,N,r_mid,h,col,plot_stuff,cut_chic)
% EXP_CHIC_STATICS  Calculate and plot casing static characteristics
%
%   [c,h,e] = EXP_CHIC_STATICS(chic_name,N,r_mid,h,C,plot_stuff,cut_chic)
%
%   directory - string of chic directory location
%   chic_name - string of chic filename
%   N - channel structure pointing to useful channels
%   r_mid - midspan rig radius in metres
%   h - figure handle for multiple plots
%   col - 3 value RGB colour array
%   plot_stuff - 0 or 1 for showing working
%   cut_chic - 0 or 1 for cutting chic to useful area

% Default to show plots
if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end

% Default to open new window
if plot_stuff == 1 && (exist('h','var') == 0 || isempty(h) == 1)
    h = figure('position',[21 85 806 883]); hold on; grid on; box on;
    xlabel('Flow Coefficient'); ylabel('Pressure Rise Coefficient');
end

% Default colour to black
if exist('col','var') == 0 || isempty(col) == 1
    col = [0 0 0];
end

% Default to cutting
if exist('cut_chic','var') == 0
    cut_chic = 1;
end

% Determine type of characteristic to plot and load data
load([directory chic_name]);

% Load in air properties
exp_air

% Interpolate all data in time to desired rate
e = exp_chic_interp(e);

% Join together temperature measurements from RTD and themocouples
e.T_interp = [e.T_interp e.Q];

% Calculate rig operating point
e.N = N.rig; e.P = e.P_interp; e.T = e.T_interp;
e = exp_rig_calculate(e,r_mid);

% Record time vector
c.time = e.time_interp;

% Get operating points
c.nrt = e.nrt; c.mrtp = e.mrtp; c.phi = e.phi;
ro = e.ro; U_mid = e.U_mid; 

% Calcualate pressure rise coefficients
P_abs = e.P_interp - repmat(e.P_interp(:,N.Pa),[1 size(e.P_interp,2)]) + e.Pa;
c.psi = (mean(P_abs(:,N.P_out),2) - mean(P_abs(:,N.Po_in),2)) ./ (0.5 * ro .* U_mid.^2);
if isfield(N,'P_Sin') == 1
    c.psi_rotor= (P_abs(:,N.P_Sin) - P_abs(:,N.P_Rin)) ./ (0.5 * ro .* U_mid.^2);
    c.psi_stator = (P_abs(:,N.P_out) - P_abs(:,N.P_Sin)) ./ (0.5 * ro .* U_mid.^2);
end
c.Rp = mean(P_abs(:,N.P_out),2) ./ mean(P_abs(:,N.Po_in),2);

% Arrange chics into individual stable lines padded by NaNs
if cut_chic == 1
    c = exp_chic_cut(c);
end

% Save raw data
c.P = P_abs; c.ro = ro; c.U_mid = U_mid; c.T = e.T_interp; c.ro = ro; 
if isfield(e,'V_interp') == 1; c.V = e.V_interp; end;

% Determine whether to plot or not
if plot_stuff == 1

    % Plot results
    figure(h)
    l = plot(c.phi,c.psi,'-','Color',col,'LineWidth',2);
    if isfield(N,'P_Sin') == 1
        plot(c.phi,c.psi_rotor,'-','Color',col,'LineWidth',2);
        plot(c.phi,c.psi_stator,'-','Color',col,'LineWidth',2);
    end
end

end

