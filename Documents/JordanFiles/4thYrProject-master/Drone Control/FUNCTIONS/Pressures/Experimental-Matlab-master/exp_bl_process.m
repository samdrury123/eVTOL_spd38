function g = exp_bl_process(filename,r_hub,r_cas,T_name,N_ref,gas,plot_stuff)
% Calculate boundary layer parameters from a pitot probe measurement

if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end

% Load in file
load(filename);

% Calculate pressure coefficient
T_ref = e.T_raw(:,N_ref.T); P_ref = e.P_raw(:,N_ref.P) + e.Pa;
ro_ref = P_ref ./ (gas.R * T_ref);
U_mid = 0.5 * (r_hub + r_cas) * 2 * pi * e.rpm / 60;
Cp = (e.P_raw(:,e.N.PoP)+e.Pa - P_ref) ./ (0.5 * ro_ref * U_mid^2);

% Ask for point closest to the wall if not already saved
if isfield(e,'l_wall') == 0
    figure(); hold on; grid on; box on;
    plot(Cp,-e.Iota,'b.-'); v = axis; axis([v(1) v(2) 0 2]);
    l_wall = str2double(inputdlg('Input Wall Index'));
    e.l_wall = l_wall;
    save(filename,'e')
end

% Calculate wall distance and extrapolate pressure near wall
Cp = Cp(e.l_wall:end);
d = -[0 ; cumsum(2*pi*diff(e.Iota(e.l_wall:end)) * e.l_probe / 360)];
d = d + e.probe.D/2;
p = polyfit(d(1:3),Cp(1:3),1);
d = [0 ; d];
Cp = [polyval(p,0) ; Cp];

% Plot pressure distribution
% if plot_stuff == 1
%     figure(); hold on; grid on; box on;
%     plot(Cp,d,'b.-');
% end

% Calculate velocity profile
Po = Cp * (0.5 * mean(ro_ref) * U_mid^2) + mean(P_ref);
P = Po(1);
ro = P / (gas.R * mean(e.T_raw(e.l_wall:end,e.N.(T_name))));
V = ((Po - P) / (0.5 * ro)).^0.5;

% Plot boundary layer velocity profile
if plot_stuff == 1
    figure(); hold on; grid on; box on;
    plot(V,d,'b.-')
end

% Calculate boundary layer parameters
V_max = mean(V(end-2:end));
dst = sum( (1 - 0.5 * (V(1:end-1) + V(2:end)) / V_max) .* diff(d));
th = sum( (0.5 * (V(1:end-1) + V(2:end)) / V_max) .* (1 - 0.5 * (V(1:end-1) + V(2:end)) / V_max)...
    .* diff(d));
de = sum( (0.5 * (V(1:end-1) + V(2:end)) / V_max) .* (1 - (0.5 * (V(1:end-1) + V(2:end)) / V_max).^2)...
    .* diff(d));
H = dst / th;

% Record results
g.Cp = Cp;
g.d = d;
g.V = V;
g.dst = dst;
g.th = th;
g.H = H;
g.de = de;
g.H32 = de /th;
g.r_nondim = (e.r(1) - r_hub) / (r_cas - r_hub);

end