%% Light-weight Design Code for Ducted Fan Blades
% Master
clear all;
close all;

global plotflag crdfflag
plotflag = 1;
crdfflag = 0;
contourflag = 0;

%% Set globals
rho = 1.225;

phi.m = 0.8; %0.8
psi.m = 0.25; %0.25

if crdfflag==1
    phi.m = 0.8;
    psi.m = phi.m*phi.m*0.125/(0.8*0.8);
end

psi_ts.m = 2*psi.m - phi.m^2;      % Total to Total stage loading = p02 - p01 / (0.5*rho*U^2)

rc = 60e-3;
rh = 20e-3;

R.AR = 2.5; R.phi = phi.m; R.psi = psi.m;
S.AR = 2.5; S.phi = phi.m; S.psi = psi.m;

%% Mass Model - Set Omega

n = 500;
meshv.omega = zeros(n,n);

switch contourflag
    case 1
        [meshv.phi, meshv.psi] = meshgrid(linspace(0,1,n),linspace(0,0.5,n));
        [meshv.propmass, meshv.diff] = MassModel(meshv.phi, meshv.psi, rc, rh, 0);
        meshv.omega = ThrustWeight(meshv.phi, meshv.psi, rc, rh, meshv.propmass);        
        meshv.rpm = meshv.omega*30/pi;
        
        meshv.psi(meshv.rpm > 20000) = NaN;
        meshv.phi(meshv.rpm > 20000) = NaN;
        meshv.rpm(meshv.rpm > 20000) = NaN;
        meshv.omega(meshv.rpm > 20000) = NaN;
        meshv.diff.L(meshv.rpm > 20000) = NaN;
        meshv.diff.rce(meshv.rpm > 20000) = NaN;
        meshv.diff.rhe(meshv.rpm > 20000) = NaN;
        meshv.propmass.total(meshv.rpm > 20000) = NaN;
        meshv.propmass.dhub(meshv.rpm > 20000) = NaN;
        meshv.propmass.case(meshv.rpm > 20000) = NaN;
        
        meshv.fom = sqrt(meshv.phi.^2 ./ meshv.psi);
        meshv.power = rho*pi.*meshv.phi.*meshv.psi.*meshv.omega.^3.*(rc^2 - rh^2)*((rc^2 + rh^2)/2)^(1.5);
    case 0
end

[propmass, diff] = MassModel(phi.m, psi.m, rc, rh, 0);
omega = ThrustWeight(phi.m, psi.m, rc, rh, propmass);
rpm = omega*30/pi;
% omega = RPM2RADS(rpm);

if crdfflag == 1
    omega = 2*omega;
    rpm = 2*rpm;
end

%% Continue 
filepath = [pwd, '/Geometry/'];

R.type = 'rotor';
S.type = 'stator';

DESIGNSECTIONS = 5;

%% VORTEX CONDITION
% 'free'; 'forced'; 'constangle';

p = 'custom';

if strcmp(p, 'custom')
    pp = 1.7;
else
    pp = 0;
end

%% Set PSI distribution

% 1 = MASS AVERAGED
% 2 = FRACTION OF SPAN
% 3 = PSI AVERAGED
% 4 = PSI*RADIUS AVERAGED

R_POS = 1;

switch R_POS
    case 1 % MASS AVERAGE
        rr = [rc rh];
        rm = rms(rr);
        %rm = rms(rr);
    case 2 % FRACTION OF SPAN
        midlineposition = 0.5;                 % Fraction of span
        rm = rh + (rc-rh)*midlineposition;
    case 3 % PSI AVERAGE
        rm = sqrt(2*(rc^2 * rh^2)/(rc^2 + rh^2));
    case 4 % PSI*RADIUS AVERAGE
        rm = 2*(1/rc + 1/rh)^(-1);
end

%% OPTINAL LOAD EXISTING META DATA
% load('JONNY_RIG.mat');
% omega = 5000*pi/30;
%% Calculate flow angles and velocities
[V, ang, phi, psi, psi_ts, radius, sections, reaction, pc] = VelocitiesAngles(DESIGNSECTIONS, omega, phi, psi, psi_ts, rc, rh, rm, p, pp);

%% Calculate delta, metal angles, and pitch to chord ratio
[V, ang, R, S, carter, delta] = Deviation('blade', V, ang, phi, psi, R, S, rc, rh, rm);

%% Energy 
[power, thrust, torque, mass, phi, psi, FOM, P] = Energy(V, phi, psi, radius, omega, rho, rc, rh);

%% Assemble blade section variables
[R,S] = Assemble(R,S,ang,sections,radius, rc, rh, rm);

%% Blade Shapes
R = Blades(R, p);
S = Blades(S, p);

%% Reynolds number
Re.m = V.m.rel1.*R.m.chord / (1.48e-5);
Re.sec = V.sec.rel1.*R.sec.chord / (1.48e-5);
Re.span = V.span.rel1.*R.span.chord / (1.48e-5);

%% Plot
if plotflag == 1
    x = PlotFlow(rc, rh, radius, sections, phi, psi, psi_ts, delta, R, S, carter, V, ang, omega, reaction, contourflag, meshv);
else
    x = (radius-rh)./(rc-rh);
end

sigma = SIG(phi.m, psi.m);