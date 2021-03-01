%% Global Variables

% Mean-line Design
r_h             = 11.6e-3; % m Hub radius
r_c             = 35e-3; % m Casing radius

DF_r            = linspace(0.45,0.40,5)';
DF_s            = linspace(0.40,0.40,5)';
p               = 1.1; % psi distribution exponent (0 - forced, 2 - free)

m               = 100;

Vinf_range      = linspace(0,19.3558,m);

for index = 1:m
Vinf            = Vinf_range(index);    
    
    
%% CV analysis

% Gas & Wind Tunnel
rho             = 1.18;
Pinf            = 101048; % Gauge pressure only

% Aspect Ratio
AR_rotor        = 1;
AR_stator       = 0.75;

% Efficiency
eta.motor       = 0.9;

% Motor Properties
rpm             = 18e3;
omega           = rpm * 2 * pi / 60; % rad/s

% Geometry
r_m             = (r_h + r_c)/2;

% Areas
A0              = pi .* (r_c.^2);
A1              = pi .* (r_c.^2 - r_h.^2);
A2              = A1;
A3              = A2;

% Chord
Chord_r         = (r_c-r_h)/AR_rotor;
Chord_s         = (r_c-r_h)/AR_stator;

%% Radial Equilibrium
n               = 9999;

[V_x, phi, psi] = RadialEquilibriumFunction(r_h, r_c, rpm, phi_mid, psi_mid ,n, p );


%% Thrust and Power
U_m             = r_m * omega;

r_temp          = linspace(r_h,r_c,n)';
dr              = r_temp(2) - r_temp(1);
U               = r_temp .* omega;
dh_0            = psi .* U.^2;

% Reset values
mdot            = 0;
dH_0            = 0;

% Integrate
for i = 1:n
    mdot            = mdot + 2 * pi * rho * V_x(i) * r_temp(i) * dr;
    dH_0            = dH_0 + dh_0(i) * rho * V_x(i) * 2 * pi * r_temp(i) * dr;
   
end
mdot        = 0.089;

phi_m       = 0;
psi_m       = 0;

for i = 1:n
    % mass averaged phi and psi
    psi_m           = psi_m + (1/mdot) * dh_0(i)/U(i)^2 * rho * V_x(i) * 2 * pi * r_temp(i) * dr;
    phi_m           = phi_m + (1/mdot) * rho * V_x(i)^2 / U(i) * 2 * pi * r_temp(i) * dr;        
end

V0              = mdot / (rho * A0); % uniform inlet flow
P00             = Pinf + 0.5 * rho * Vinf^2;

dh_0_m          = dH_0 / mdot; % mass averaged stagnation enthalpy change

P01             = P00;

P4              = Pinf;

V1              = mdot / (rho * A1); % assume no swirl at inlet to rotor
V3              = mdot / (rho * A3); % assume no swirl at exit to stator
V4              = V3 * NR; % Set by continuity

P04             = P4 + 0.5 * rho * V4^2;
P03             = P04;
P02             = P03;
dP_0_m          = P02 - P01;

eta.flow(index) = dP_0_m / (rho * dh_0_m);

P3              = P03 - 0.5 * rho * V3^2;
P1              = P01 - 0.5 * rho * V1^2;
P0              = P00 - 0.5 * rho * V0^2;

% Exit Area
A4              = A3/NR;

%% Thrust and Power

Thrust(index)   = mdot * (V4 - Vinf);
Thrust2(index)  = mdot * (V4 - V0) + P4 * A4 + Pinf * (A0 - A4) - P0 * A0;

Power_KE        = 0.5 * mdot * (V4^2 - Vinf^2);
Power(index)    = dH_0 / eta.motor;

end
L               = 0.120;
D               = 0.078;
nu              = 1.53e-5;

Re_theory       = Vinf_range .* L / nu;
CF_theory       = Thrust ./ (0.5 * rho * U_m^2 * D^2);

Re_theory3       = Vinf_range .* L / nu;
CF_theory3       = Thrust2 ./ (0.5 * rho * U_m^2 * D^2);


%% Method 2

A2              = pi * ((56e-3)/2)^2;

V1              = Vinf_range;
P_in            = 60;
P               = P_in .* eta.flow;

for i  = 1:m
    coef            = [1/(rho*A2)-2*V1(i)^2 , 6*P(i)*V1(i), -(2*P(i))^2];
    temp            = roots(coef);
    F1(i)           = temp(1);
    F2(i)           = temp(2);
    V2_1(i)         = 2*P(i) / F1(i) - V1(i);
    V2_2(i)         = 2*P(i) / F2(i) - V1(i);
    dV_1(i)         = V2_1(i) - V1(i);
    dV_2(i)         = V2_2(i) - V1(i);
    
    if dV_1(i) >= 0
        F(i)        = F1(i);
        V2(i)       = V2_1(i);
    elseif dV_2(i) >= 0
        F(i)        = F2(i);
        V2(i)       = V2_2(i);
    end
end

Re_theory2      = L * V1 / nu;
CF_theory2      = F ./ (0.5 * rho * U_m^2 * D^2);

% figure
% hold on
% plot(Re_theory,CF_theory);
% plot(Re_theory2,CF_theory2);
% plot(Re_theory3,CF_theory3);
% legend('Method 1','Method 2')
% grid on
% hold off

