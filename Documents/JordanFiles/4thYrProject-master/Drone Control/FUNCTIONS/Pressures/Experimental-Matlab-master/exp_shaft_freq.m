function [w,V,dtime] = exp_shaft_freq(s,n,plot_stuff,V,dtime)
% EXP_SHAFT_FREQ  Calculate shaft frequency in radians per second
%
%   s - connection data structure
%   plot_stuff - 0 or 1 for showing working
%   w - shaft frequency in radians per second
%   V - raw readings of voltage
%   dtime - raw readings in time

% Default no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Open figure window
if plot_stuff == 1
    figure(); hold on; grid on; box on;
    xlabel('Time / s'); ylabel('Voltage / V');
end

% Log the pxie if the voltage isn't already supplied
if exist('V','var') == 0 && exist('dtime','var') == 0
    [V,dtime] = exp_pxie_read(s);
end

% Isolate once per rev signal
V_shaft = V(:,n);

% Calculate guess shaft frequency from fft
% f = shaft_freq(s,t,V,plot_stuff);

% Check for maximum difference in voltage is greater than 2 volts
if max(V_shaft) - min(V_shaft) < 2
%     warning('Rig may be turned off or sample period is too small')
    w = 0;
    return
end

% Find shaft edges
% V_shaft = V_shaft - 0.5*max(V_shaft);
V_shaft = V_shaft - 0.5 * (max(V_shaft) + min(V_shaft));
i = find(V_shaft(2:end) .* V_shaft(1:end-1) < 0 & V_shaft(2:end) < 0);

% Plot time domain signal
if plot_stuff == 1
    plot(dtime,V_shaft,'b-');
end

% Calculate guess shaft frequency from maximum period between falling edges
if numel(i) > 1
    f = 1 / max(diff(dtime(i)));

    % Remove falling edges whose period is too short
    dt = diff(dtime(i)); 
    i(dt < 0.5 * 1/f) = [];

    % Calculate mean falling edge period and shaft speed
    f = 1 / trimmean(diff(dtime(i)),12);
    w = f * 2 * pi;

    % Plot falling edge points
    if plot_stuff == 1
    %     subplot(2,1,1);
        plot(dtime(i),V_shaft(i),'ko');
    end
   
%     figure(); hold on; grid on; box on;
%     plot(dtime(i),w);
   
    
    
    
else
    w = 0;
%     warning('Rig may be turned off or sample period is too small')


end

end
