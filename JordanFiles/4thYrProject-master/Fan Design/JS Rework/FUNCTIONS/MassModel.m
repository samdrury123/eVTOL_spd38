function [mass, diff] = MassModel(phi, psi, rc, rh, chassis)
%% Mass Model v2.0
% - Inlet
% - Blade passage
% - Blades
% - Diffuser hub
% - Diffuser casing
%
% ESDU DATA
% L/h1 = [1.21 1.28 1.37 1.46 1.56 1.70 1.90];
% Ae/A = [2.604 3.432 4.383 5.412 6.583 8.11 10];
% C_pr = [0.3 0.35 0.4 0.45 0.5 0.55 0.6];

rhoPLA = 1240; % kg/m^3

sigma = sqrt(phi.^2 ./ (2*psi));

rce = ((rc+rh)/2) + (sigma/2) .* (rc-rh);
rhe = ((rc+rh)/2) - (sigma/2) .* (rc-rh);

if sigma > 1
%     a = 0.3313;
%     b = 0.5275;
%     c = 1;
%     L = (rc-rh) .* (11.23*sigma -10.99);
    L = (rc-rh) .* (-5.56*sigma.^3 + 24.16*sigma.^2 - 23.41*sigma + 5.413);
else
    L = 0;
end

t = 2e-3; % ie 2mm

%% Inlet
mass.inlet = 3.90 * (rc - rh) .* (3*rc - rh + 0.001);

%% Blade Passage
mass.annul = 2.73*rh + 0.818*rc - 0.0062;

%% Blades
mass.blade = 17.5*(rc.^2 - rh.^2);

%% Diffuser Hub
R = rh;
r = rhe;
S = rh - t;
s = rhe - t;
h = L;

mass.dhub = rhoPLA * (h*pi/3) .* (R.^2 + R.*r + r.^2 - S.^2 - S.*s - s.^2);

%% Diffuser Casing
R = rce + t;
r = rc + t;
S = rce;
s = rc;
h = L;

mass.case = rhoPLA * (h*pi/3) .* (R.^2 + R.*r + r.^2 - S.^2 - S.*s - s.^2);

%% Motor
mass.motor = 0.087;

%% Total
if chassis == 0
    mass.chassis = 0.529 + 0.466; % Dry weight + battery
else 
    mass.chassis = chassis;
end

mass.total = mass.chassis + 4*(mass.inlet + mass.annul + mass.blade + mass.dhub + mass.case + mass.motor);

mass.prop = mass.chassis + 4*(mass.motor + mass.blade);

diff.L = L;
diff.rce = rce;
diff.rhe = rhe;

mass.total(mass.total < 0) = NaN;
mass.dhub(mass.total < 0) = NaN;
mass.case(mass.total < 0) = NaN;

diff.L(mass.total < 0) = NaN;
diff.rce(mass.total < 0) = NaN;
diff.rhe(mass.total < 0) = NaN;