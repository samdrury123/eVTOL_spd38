clear all
close all

%% Global Variables

% Mean-line Design
r_h             = 11.6e-3; % m Hub radius
r_c             = 35e-3; % m Casing radius
NR              = 1.35;

psi_mid         = 0.39;
phi_mid         = 0.40;  
DF_r            = linspace(0.45,0.40,5)';
DF_s            = linspace(0.40,0.40,5)';
p               = 1.1; % psi distribution exponent (0 - forced, 2 - free)

Vinf            = 0;

m               = 100;
rpm_range       = linspace(0e3,18.850e3,m);

for index = 1:m
    
%% CV analysis

% Gas & Wind Tunnel
rho             = 1.18;
Pinf            = 101048; % Gauge pressure only

% Aspect Ratio
AR_rotor        = 1;
AR_stator       = 0.75;

% Efficiency
eta.motor       = 0.98;

% Motor Properties
rpm             = rpm_range(index);
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
Thrust2         = mdot * (V4 - V0) + P4 * A4 + Pinf * (A0 - A4) - P0 * A0;

Power_KE        = 0.5 * mdot * (V4^2 - Vinf^2);
Power(index)    = dH_0 / eta.motor;

end

figure
plot(rpm_range,Thrust)
figure
plot(rpm_range,Power)