function [V,ph] = exp_imbalance_log(s,N,rpm_target,edge_type,plot_stuff,fourier_mode)
% Log accelerometer and shaft signal for a given duration and calculate imblance force and phase

% Default to log falling edges
if exist('edge_type','var') == 0
    edge_type = 'falling';
end

% Default to show plots
if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end

% Default to use fourier coefficients
if exist('fourier_mode','var') == 0
    fourier_mode = 1;
end

% Set the shaft speed
rpm = 0;
while abs(rpm - rpm_target) / rpm_target > 0.05
    rpm = exp_shaft_freq(s,1,0) * 30 / pi;
    disp(['RPM = ' num2str(rpm)]);
end

% Log pxie for desired duration
[q,t] = s.startForeground();

% Record shaft and accelerometer voltages
q_shaft = q(:,N.V_shaft);
q_accel = q(:,N.V_accel);

% Normalise once per rev voltage
q_shaft = ((q_shaft - min(q_shaft)) / (max(q_shaft) - min(q_shaft))) > 0.5;

% Calculate shaft speed from rising or falling edges
if strcmp(edge_type,'falling') == 1
    i = find(q_shaft(1:end-1) .* (q_shaft(2:end)-1) == -1); t_edge = t(i);
elseif strcmp(edge_type,'rising') == 1
    i = find((q_shaft(1:end-1)-1) .* q_shaft(2:end) == -1); t_edge = t(i);
end

% Find average shaft frequency
f = 1 / mean(diff(t_edge));

% Ensemble average acceleration data
ni = 3600; q_interp = zeros(ni,length(i)-1);
for n = 1:length(i)-1
    q_interp(:,n) = interp1(linspace(0,1,i(n+1)-i(n)+1),q_accel(i(n):i(n+1)),linspace(0,1,ni));
end
q_av = mean(q_interp,2);

% Find phase and amplitudes
if fourier_mode == 0
    
    % Find values from turning points of ensemble average
    [q_max,i_max] = max(q_av); q_min = min(q_av);
    V = q_max - q_min;
    ph = 360 * (ni - i_max) / ni;
    
else
    
    % Find values from fourier coefficients
    x = linspace(0,2*pi,ni)';
    a = trapz(x,q_av .* cos(x) / pi); b = trapz(x,q_av .* sin(x) / pi);
    V = 2 * hypot(a,b); 
    ph = rad2deg(atan2(b,a)); q = ph < 0; ph(q) = ph(q) + 360;
    ph = 360 - ph;
    
end

% Plot time series data
if plot_stuff == 1
    figure(); hold on; grid on; box on;
    l(1) = plot(t,q_accel - mean(q_accel)); l(2) = plot(t,q_shaft);
end

% Generate fitted sine wave from data
t_new = t(i(1):i(end));
q_fit = 0.5 * V * cos(2 * pi * (t_new - t_new(1)) * f + deg2rad(ph)); 
[~,i_max] = max(q_fit);

% Plot the fitted sine wave and the peak 
if plot_stuff == 1
    v = axis; 
    l(3) = plot(t_new,q_fit);
    legend(l,'Raw Accel','Shaft Signal',' Fitted Accel');
    plot([t_edge t_edge],v(3:4),'k-')
    plot([t_new(i_max) t_new(i_max)],v(3:4),'k-')
    axis(v)
end

% Print output values
if plot_stuff == 1
    disp(['P-P imbalance = ' num2str(V)])
    disp(['Phase imbalance = ' num2str(ph)])
end

end