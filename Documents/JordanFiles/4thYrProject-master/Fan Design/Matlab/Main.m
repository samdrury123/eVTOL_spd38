clear all
close all
clc

%% Global Variables
% Declare global variables
global flagplot ac bc cc acr bcr ccr shiftvary sectionlocation

% Set global flags
flag            = 2;                        % 1 - variable NR, 2 - fixed NR
flagplot        = 1;                        % 0 - no plots, 1 plots

% Test points
N               = 5;                        % number of points (3,5,7)

% Drone empty weight
payload = 0;

%% Get geometry and flow coefficients from Optimise Mass
[r_c, r_h, phi, psi, omega, thrust_p, sigma] = OptimiseMass(payload);

if flag == 2
    NR  = sigma;
end
%% Geometry
% Blade shapes
% r_h             = 25e-3;                    % Hub radius (m)
% r_c             = 55e-3;                   % Casing radius (m)
r_m             = (r_h+r_c)/2;              % Mean radius (m)
span            = r_c-r_h;                  % Blade span (m)
r.radius        = linspace(r_h,r_c,N)';     % Radial test points (m)
HTR             = r_h/r_c;                  % Hub to tip ratio
AR_rotor        = 1.5;                      % Aspect ratio of rotor
AR_stator       = 1.5;                      % Aspect ratio of stator
tip_gap_percent = 0.0;                        
r.I1            = linspace(-5,-5,N)';       % Rotor incidence angle
s.I2            = linspace(-5,-5,N)';       % Stator incidence angle
ac_r            = 0.5;%???????????????????? % Max Camber position to chord ratio

% Areas
A0              = pi.*(r_c.^2);
A1              = pi.*(r_c.^2-r_h.^2);
A2              = A1;
A3              = A2;

% Chord
Chord_r         = (r_c-r_h)/AR_rotor;
Chord_s         = (r_c-r_h)/AR_stator;

% Stator Shape
theta1 = 25*pi/180;
theta2 = 10*pi/180;
ac = (tan(theta1)+tan(theta2))/(2*r_h-2*r_c);
bc = tan(theta1)-2*ac*r_h;
cc = -ac*r_h^2-bc*r_h;

% Rotor Shape
rtheta1 = -25*pi/180;
rtheta2 = -60*pi/180;
rtheta2 = -20*pi/180;
acr = (tan(rtheta1)+tan(rtheta2))/(2*r_h-2*r_c);
bcr = tan(rtheta1)-2*acr*r_h;
ccr = -acr*r_h^2-bcr*r_h;

% File Name
temp            = {'BladeGeom'};
pathname        = [pwd, '/Geometry'];
file_name       = fullfile(pathname,temp);

%% Air Conditions in Test
rho             = 1.225;                    % Density (Kg/m^3)
Pinf            = 1.01325*1e5;              % Air gauge pressure only (Pa)

%% Motor Properties
% motor.rpmax     = 8000;                     % Motor operating speed(RPM)
% motor.rpm       = motor.rpmax;
% omega           = motor.rpm*2*pi/60;        % (rad/s)
U_m             = omega*r_m;

%% Assumed Parameters
% Performance parameters
% thrust_p        = thrust_p;                 % Required propulsor thrust = Mg/4
% V_2             = sqrt(thrust_p/(rho*A2));  % Required theoretical exit velocity
eta.motor       = 0.9;                      % Assumed motor efficiency 
eta.inlet       = 1;                        % Assumed stag pressure ratio
eta.stator      = 1;                        % Assumed stag pressure ratio
eta.nozzle      = 1;                        % Assumed stag pressure ratio
eta.flow        = 0.7;                      % Isentropic propulsor efficiency 
psi_mid         = psi;                      % Assumed stage loading coefficient
phi_mid         = phi;                      % Assumed flow coefficient
DF_r            = linspace(0.45,0.45,N)';   % Assumed stator diffusion factor
DF_s            = linspace(0.45,0.45,N)';   % Assumed rotor diffusion factor
p               = 1.75;                      % psi distribution exponent (0 - forced, 2 - free)
% External flow conditions
Vinf            = 0;                        % Forward velocity of propulsor

%% Mean-line Design
Mean_Line_Design;

%% Metal Blade Angles 
MetalAngles;

%% Print Blades
PrintGeometry;

figure; hold on; grid on; box on; axis equal;
mesh(squeeze(r_XYZ(:,1,:)),squeeze(r_XYZ(:,3,:)),squeeze(r_XYZ(:,2,:)));
mesh(squeeze(s_XYZ(:,1,:)),squeeze(s_XYZ(:,3,:)),squeeze(s_XYZ(:,2,:)));