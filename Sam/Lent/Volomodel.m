clear all

m = 900;        % Aircraft weight, kg
v = 25;         % Flight speed, m/s
th = 60;        % Time spent hovering,s
Cd = 0.4;       % Drag coefficient in forward flight COMPLETE GUESS
Af = 15;        % Frontal area, m^2 GUESS


bm = 300;       % Battery Mass, kg
be = 250;       % Batter energy density Wh/kg
SOCmin = 0.1;

n = 18;         % Number of props
d = 2.3;        % Prop diameter

ro = 1.225;     % Air density kg/m^3
g = 9.81;       % acceleration due to gravity, m/s^2

etaA=0.8;       % Aerodynamic efficiency
etaM=0.95;      % Motor efficiency
etaE=0.95;      % Electric efficiency
etaB=0.85;      % Battery efficiency

% Cruise flight model
theta = atand((0.5*ro*Af*Cd*v^2)/(m*g));                % Angle to balance weight and drag

% Propulsor analysis
A = n*pi*(d/2)^2;                                             % Disc area of propulsors, m^2
Thc = m*g /cosd(theta);                                     % Thrust in cruise, N
Pc = (Thc^(3/2) / (2*ro*A)^0.5) / (etaA*etaM*etaE*etaB);    % Power in cruise, W
Thh = m*g;                                                  % Thrust in hover, N
Ph = (Thh^(3/2) / (2*ro*A)^0.5) / (etaA*etaM*etaE*etaB);    % Power in hover, W

% Overall power and flight time
Eb = bm*be*(1-SOCmin)*3600;                                      % Energy available from battery Wh
tc = (Eb - Ph*(th))/Pc;                           % Flight time, mins
R = tc*v / 1000;                                         % Range, km

