function e = exp_hw_calculate(e,c,N,x,plot_stuff)
% EXP_HW_CALCULATE  Apply calibration curve to hot-wire voltages

% Default to not show calibration curve
if exist('plot_stuff','var') == 0; plot_stuff = 0; end;

% Interpolate rig temperature at traverse grid
time_rig = datenum(e.rig.time); time_trav = datenum(reshape(e.time,[],6));
T = interp1(time_rig,mean(e.rig.T(:,N.T),2),time_trav,'linear','extrap');
T = reshape(T,size(e.r));

% Correct amplified AC signal back to standard level
e.V_hw = e.V_hw / e.probe.gain + repmat(e.V_av,[1 1 1 size(e.V_hw,4)]);

% Correct hotwire voltages for temperature drift
Rt_fac = abs(c.T_over ./ (c.T_over + c.T_std - T)).^0.5;
e.V_hw = e.V_hw .* repmat(Rt_fac,[1 1 1 size(e.V_hw,4)]);

% Calculate velocity from voltage reading and calibration curve
e.V = exp(polyval(c.p_king,log(e.V_hw - mean(e.V0_start))));
e.V = polyval(c.p_king,e.V_hw.^2).^(1 / c.n);

% Record axial coordinates
e.x = ones(size(e.r)) * x;

% Plot velocities on calibration curve
if plot_stuff == 1
    figure(); hold on; grid on; box on; xlabel('Voltage / V'); ylabel('Velocity / ms^{-1}');
    plot(mean(c.V(:,:,2),2),c.Vx,'.-');
    i = unique(round(linspace(1,numel(e.V_hw),1000)));
    plot(e.V_hw(i),e.V(i),'.');
end


end

