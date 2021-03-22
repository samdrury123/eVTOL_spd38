function g = ts_calc_guess(g,rpm_max)
% TS_CALC_GUESS  Generate a 1D initial guess based on metal angles and flow properties

% Get machine rpm from maximum of all rows
if exist('rpm_max','var') == 0
    rpm_max = 0; rpm_names = {'rpm' 'rpmi1' 'rpmi2' 'rpmj1' 'rpmj2' 'rpmk1' 'rpmk2'};
    for bid = ts_get_bids(g)
        for v = 1:length(rpm_names)
            rpm_max = max(rpm_max,g{bid+1}.bv.(rpm_names{v}));
        end
    end
end

% Initialise compressor parameters
phi = 0.5; 

% Calculate gas properties
cp = g{1}.av.cp; ga = g{1}.av.ga;
cv = cp / ga; R = cp - cv;

% Get inlet boundary condition and axial coordinates of important planes
x_mp = []; A_mp = [];
for bid = ts_get_bids(g)
    for pid = ts_get_pids(g,bid)
        if g{bid+1}.patch{pid+1}.attribute.kind == 0
            nj = g{bid+1}.attribute.nj;
            r = g{bid+1}.r(1,:,1); r_mid = 0.5 * (r(1) + r(end));
            Po_in = interp1(r,g{bid+1}.patch{pid+1}.pstag(1:nj),r_mid);
            To_in = interp1(r,g{bid+1}.patch{pid+1}.tstag(1:nj),r_mid);
            x_in = min(g{bid+1}.x(:)); A_in = pi * (r(end)^2 - r(1)^2);
        elseif g{bid+1}.patch{pid+1}.attribute.kind == 1
            x_out = max(g{bid+1}.x(:)); 
            A_out = pi * (g{bid+1}.r(end,end,1)^2 - g{bid+1}.r(end,1,1)^2);
        elseif g{bid+1}.patch{pid+1}.attribute.kind == 2 && g{bid+1}.patch{pid+1}.attribute.ist == 0
            x_mp = [x_mp min(g{bid+1}.x(:))];
            A_mp = [A_mp pi * (g{bid+1}.r(1,end,1)^2 - g{bid+1}.r(1,1,1)^2)];
        end
    end
end

% Assemble axial coordinates and areas
x = [x_in x_mp x_out]; A = [A_in A_mp A_out];

% Get blade metal angles
[~,bid_os,chi] = ts_find_chord(g,r([1 end])); nn = size(chi,2);

% Get rpms of spinning rows
rpm = zeros(size(bid_os));
for n = 1:length(bid_os)
    rpm(n) = g{bid_os(n)+1}.bv.rpm;
end

% Calculate relative angular velocities
Vx = phi * r_mid * 2 * pi * rpm_max / 60;
Vt_rel(1) = tand(chi(1,1)) * Vx;
Vt_rel(2:nn+1) = tand(chi(2,:)) * Vx;

% Calculate absolute angular velocities and whirl angles
U_mid = [rpm(1) rpm] * r_mid * 2 * pi / 60;
Vt = Vt_rel + U_mid; V = (Vx^2 + Vt.^2).^0.5;
Alpha = atand(Vt./Vx);

% Stagnation enthalpy using Euler's work equation
ho = zeros(1,nn); ho(1) = cp * To_in;
for n = 2:nn+1
    ho(n) = ho(n-1) + U_mid(n) * (Vt(n) - Vt(n-1));
end

% Mach number from velocity compressible relation
To = ho ./ cp; Mi = 0:0.01:2; 
Fi = (ga-1)^0.5 * Mi .* ((1 + 0.5*(ga-1) * Mi.^2) .^ (-0.5));
F = V ./ ((cp * To).^0.5);
M = interp1(Fi,Mi,F);

% Total pressure from mass flow compressible relation
Fi = (ga * (ga-1)^(-0.5)) * Mi .* ((1 + 0.5 * (ga -1) * Mi.^2).^(-0.5 * (ga+1) / (ga-1)));
m = interp1(Mi,Fi,M(1)) * A(1) * Po_in / ((cp * To_in)^0.5);
Po = m * (cp * To).^0.5 ./ (A .* interp1(Mi,Fi,M .* cosd(Alpha)));

% Static pressure and temperature from compressible relations
Fi = (1 + 0.5 * (ga - 1) * Mi.^2).^(-ga/(ga-1));
P = Po .* interp1(Mi,Fi,M);
Fi = (1 + 0.5 * (ga - 1) * Mi.^2).^(-1);
T = To .* interp1(Mi,Fi,M);

% Calculate density from ideal gas law
ro = P ./ (T * R);

% Primary flow variables
rovx = ro * Vx;
rovr = ones(size(ro));
rorvt = ro .* r_mid .* Vt;
roe = ro .* (cv * T + 0.5 * V.^2);

% Apply guess to every block
for bid = ts_get_bids(g)
    g{bid+1}.ro = interp1(x,ro,g{bid+1}.x,'linear','extrap');
    g{bid+1}.rovx = interp1(x,rovx,g{bid+1}.x,'linear','extrap');
    g{bid+1}.rovr = interp1(x,rovr,g{bid+1}.x,'linear','extrap');
    g{bid+1}.rorvt = interp1(x,rorvt,g{bid+1}.x,'linear','extrap');
    g{bid+1}.roe = interp1(x,roe,g{bid+1}.x,'linear','extrap');
    g{bid+1}.trans_dyn_vis = ones(size(g{bid+1}.ro)) * 0.0002;
end

% Set backpressure
% for bid = ts_get_bids(g)
%     for pid = ts_get_pids(g,bid)
%         if g{bid+1}.patch{pid+1}.attribute.kind == 1
%             g{bid+1}.patch{pid+1}.pout = P(end);
%         end
%     end
% end


end