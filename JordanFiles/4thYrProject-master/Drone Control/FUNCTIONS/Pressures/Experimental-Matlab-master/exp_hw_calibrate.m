function c = exp_hw_calibrate(c,N,plot_stuff,AR,n)
% Calculate King's law coefficients from hotwire calibration data

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Default to no contraction ratio
if exist('AR','var') == 0
    AR = 1;
end

% Default to velocity exponent of 0.4
if exist('n','var') == 0
    n = 0.4;
end

% Prepare calibration run for velocity calculation
c.N = N; r_mid = 1; c.w = 0;

% Calculate velocity in rig
[e,Vx] = exp_rig_calculate(c,r_mid); Vx = Vx * AR; c.Vx = Vx;

% Mean voltages
V = mean(c.V(:,:,2),2);

% Correct for atmospheric temperature
To = mean(c.T(:,N.To),2);
Rt_fac = abs(c.T_over ./ (c.T_over + c.T_std - To)).^0.5;
V = V .* Rt_fac; c.V = V;

% Calculate powers of velocity and voltage
% x = V.^2; y = (e.ro .* Vx).^n;
x = V.^2; y = Vx.^n;

% Fix line through data
c.p_king = polyfit(x,y,3); c.n = n;
x_fit = linspace(0,max(interp1(y,x,130^n,'linear','extrap'),max(x)),100);

% Calculate log-log data (King's Law)
% x = log(V(2:end)-V(1));
% y = log(Vx(2:end)); 

% Calculate polynomial fit of log-log data
% c.p_king = polyfit(x,y,3);
% x_fit = linspace(min(x)-0.2,max(x) + 0.2,100);

% Plot calibration points and fit
if plot_stuff == 1
    figure(); hold on; grid on; box on; xlabel('Voltage'); ylabel('Velocity');
    plot(x,y,'.-');
    plot(x_fit,polyval(c.p_king,x_fit),'-');
end


end

