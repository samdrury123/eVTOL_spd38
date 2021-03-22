function t = bl_parameterise_thickness(s_cl,thick)
% BL_PARAMETERISE_THICKNESS  Decompose thickness distribution into governing parameters

% Get initial maximum value for parameters at maximum thickness
[~, i_max] = max(thick); s_fit = s_cl(i_max);

% Fit quartic through points near maximum
q = s_cl > s_fit - 0.05 & s_cl < s_fit + 0.05;
s_fit = linspace(s_fit - 0.05,s_fit + 0.05,1000).';
thick_fit = polyval(polyfit(s_cl(q),thick(q),4),s_fit);

% Recalculate maximum
[t.thick_max, i_max] = max(thick_fit); t.s_thick_max = s_fit(i_max);
    
% Scale thickness by maximum value
thick = thick / t.thick_max;

% Find where the trailing edge circle begins
dtds = abs(diff(grad_mg(s_cl,thick)));
lim = interp1(s_cl(2:end),dtds,0.85);
s_circ = s_cl(find(dtds > 5*lim & s_cl(2:end) > 0.85,1));
if isempty(s_circ) == 1; s_circ = 1; end;

% Get points in the vicinity of the trailing edge and fit quadratic
q = s_cl > s_circ - 0.025 & s_cl < s_circ - 0.005;
t.p_te = polyfit(s_cl(q),thick(q),2);

% Evaluate quadratic at the trailing edge for thickness and wedge angle
t.thick_te = polyval(t.p_te,1);
t.wedge_te = atand(-polyval(polyder(t.p_te),1));

% Interpolate chordwise spacing to be uniform around leading edge
ni_le = 10; i_le = 5;
d = [0 ; cumsum((diff(thick).^2 + diff(s_cl).^2).^0.5)];
s_le = interp1(d(1:i_le),s_cl(1:i_le),linspace(0,d(i_le),ni_le)','pchip');
thick_le = interp1(s_cl(1:i_le),thick(1:i_le),s_le,'pchip');

% Circle fit leading edge for radius of curvature
t.p_le = CircleFitByTaubin([flipud([s_le -thick_le]) ; s_le(2:end) thick_le(2:end)]);
t.rad_le = t.p_le(3);
    
% Fit quartic through points near maximum and evaulate curvature at max thickness
q = s_cl > t.s_thick_max - 0.05 & s_cl < t.s_thick_max + 0.05;
t.p_max = polyfit(s_cl(q),thick(q),4);
t.rad_thick_max = -1/polyval(polyder(polyder(t.p_max)),t.s_thick_max);

% Calculate shape space function
t.S = (thick - s_cl * t.thick_te) ./ (s_cl.^0.5 .* (1 - s_cl));

end