function [mass, diff] = MassModel_SIG(sig, rc, rh, chassis)
%% Mass Model v2.0
% [mass, diff] = MassModel_SIG(sig, rc, rh, chassis)
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

rhoPLA = 1600; % 1240 kg/m^3
r = rh./rc;
siglim = (1+r) ./ (1-r);

rce = ((rc+rh)/2) + (sig/2) .* (rc-rh);
rhe = ((rc+rh)/2) - (sig/2) .* (rc-rh);

[a,b] = size(rc);

% if a == 1 || b == 1
%     if sig > 1
%         if sig < siglim
%             sig = sig;
%         else
%             sig = 9999;
%         end
%         %     a = 0.3313;
%         %     b = 0.5275;
%         %     c = 1;
%         %     L = (rc-rh) .* (11.23*sig -10.99);
%         L = (rc-rh) .* (-5.56*sig.^3 + 24.16*sig.^2 - 23.41*sig + 5.413);
%     else
%         L = 0;
%     end
% else
for i=1:a
    for j=1:b
        if sig > 1
            if sig < siglim(i,j)
                s = sig;
            else
                s = 9999;
            end
            %     a = 0.3313;
            %     b = 0.5275;
            %     c = 1;
            %     L = (rc-rh) .* (11.23*sig -10.99);
            L(i,j) = ((rc(i,j)-rh(i,j)) .* (-5.56*s.^3 + 24.16*s.^2 - 23.41*s + 5.413));

        else
            L(i,j) = 0;
        end
    end
end
% end

t = 1e-3; % 1 mm

%% Inlet
mass.inlet = (rhoPLA/1240) .* (3.90 * (rc - rh) .* (3*rc - rh + 0.001));

%% Blade Passage
mass.annul = ((2e-3)/(5e-3))^2 .* (rhoPLA/1240) .* (2.73*rh + 0.818*rc - 0.0062); % ((2e-3)/(5e-3))^2 .*

%% Blades
mass.blade = (rhoPLA/1240) .* (17.5*(rc.^2 - rh.^2));

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
mass.motor = 0.087 * (rh./20e-3);

%% Total
if chassis == 0
    mass.chassis = 0.529 + 0.466; % Dry weight + battery
else 
    mass.chassis = chassis;
end

mass.diff = mass.dhub + mass.case;
mass.fan = mass.inlet + mass.annul + mass.blade + mass.dhub + mass.case;
mass.total = mass.chassis + 4*(mass.inlet + mass.annul + mass.blade + mass.dhub + mass.case + mass.motor);
mass.CRDF = mass.chassis + 2*(mass.inlet + mass.annul + mass.blade + mass.dhub + mass.case) + 4*(mass.motor);
mass.prop = mass.chassis + 4*(mass.motor + 12e-3); % APC PROP = 12grams

diff.L = L;
diff.rce = rce;
diff.rhe = rhe;

mass.total(mass.case < 0) = NaN;
mass.CRDF(mass.case < 0) = NaN;
mass.dhub(mass.dhub > 100) = NaN;
mass.case(mass.case < 0) = NaN;

mass.total(mass.dhub < 0) = NaN;
mass.CRDF(mass.dhub < 0) = NaN;
mass.dhub(mass.dhub < 0) = NaN;
mass.case(mass.dhub < 0) = NaN;

diff.L(mass.case < 0) = NaN;
diff.rce(mass.case < 0) = NaN;
diff.rhe(mass.case < 0) = NaN;