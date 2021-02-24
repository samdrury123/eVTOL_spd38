%% ---- FUNCTION: CONTROL VOLUME ANALYSIS FOR DUCTED FAN ---- %%
function [d,g,q] = CVanalysis_comp(d,g,q,L)

% REDO THIS
% Function using non-dimensional compressible relations to analyse ducted
% fans at the mean line

% Inputs: 
% - g: geometry of fan
% - d.Th: required thrust
% - d.u0: required forward speed
% - d.phi: flow coefficient of fan
% - d.sigma: area ratio of exit nozzle
% - q.atm: atmospheric conditions at flight altitude
% - L.L1: entropy change across blade 1
% - L.L2: entropy change across blade 2
% Outputs: 
% - d.u1, d.u4: inlet/exit velocities
% - d.Um: mean line blade speed
% - d.rpm: shaft speed
% - d.psi: non-dimensional pressure rise, total-total
% - g: updated geometry of fan
% - q: Mach numbers, temperatures and pressures (aero quantities)

%% Unpack
Tatm = q.atm.T; patm = q.atm.P;
gam = q.gam;
R = q.R;
cp = q.cp;
phi = d.phi; 
sigma = d.sigma;

% Calculate geometry
g.rcas = g.D/2; 
g.rhub = g.rcas*g.htrat;
g.rm = ((g.rcas^2+g.rhub^2)/2)^0.5;
g.A1 = pi*(g.rcas^2 - g.rhub^2);
g.A4 = g.A1 * sigma;

% Flight Mach Number and stagnation quantities
M0 = d.u0 / (gam*R*Tatm)^0.5;
p00 = patm ./ (1 + (gam-1)/2.*M0.^2).^-(gam/(gam-1));
T00 = Tatm ./ (1 + (gam-1)/2.*M0.^2).^-1;

% Set up Mach Number vector for lookups
Mach = 0.001:0.001:1;

% Non-dim impulse function for lookups
Fnd = (gam-1).^0.5 ./ gam .* (1+gam.*Mach.^2) ./ Mach .* (1 + (gam-1)/2.*Mach.^2).^-0.5;
% Non-dim mass flow functions for lookups
Mndp = gam / (gam-1)^0.5 .* Mach .* (1 + (gam-1)/2.*Mach.^2).^0.5;
Mndp0 = gam / (gam-1)^0.5 .* Mach .* (1 + (gam-1)/2.*Mach.^2).^(-0.5*(gam+1)/(gam-1));

% Assume exit pressure is atmospheric
p4 = patm;

% Guess T03 and initiate losses for iteration loops
T03 = 293; 
T04 = T03;


deltaT=1;
while abs(deltaT) > 0.001/100
    % Non-dimensional momentum equation with impulse function (see 3A3 EP2,Q1)
    Thcalc = g.A4.* (Fnd .* Mndp .* p4 - patm) - Mndp .* (cp*T04).^-0.5 .* p4 .* g.A4 .* d.u0;
    M4 = interp1(Thcalc,Mach,d.Th);
    Fnd4 = (gam-1).^0.5 ./ gam .* (1+gam.*M4.^2) ./ M4 .* (1 + (gam-1)/2.*M4.^2).^-0.5;
    
    % Exit conditions 
    p04 = p4 ./ (1 + (gam-1)/2.*M4.^2).^-(gam/(gam-1));
    T4 = T04 ./ (1 + (gam-1)/2.*M4.^2);
    u4 = M4.*(gam*R*T4)^0.5;
    
    % Inlet and Exit duct isentropic
    p03 = p04;
    p01 = p00;
    T01 = T00;

    % Calculate inlet Mach Number, M1
    Mndp4 = gam / (gam-1)^0.5 .* M4 .* (1 + (gam-1)/2.*M4.^2).^0.5;
    Mndp01 = Mndp4 .* p4./p01 .* (T01./T04).^0.5 .* sigma;
    M1 = interp1(Mndp0,Mach,Mndp01);
    T1 = T01 .* (1 + 0.5*(gam-1).*M1.^2).^-1;
    u1 = M1 .* (gam*R*T1)^0.5;
    
    % Blade speed from flow coefficient
    Um = u1/phi;
    rpm = Um / (2*pi*g.rm) * 60;

    % Stage loading
    psi = cp*(T04-T01) / Um^2;

    % Use loss to get T03
    T03calc = T01* exp((L.L1 + L.L2 + R*log(p03/p01))/cp);
    
    % Isentropic case
    T03s = T01 .* (p03/p01).^((gam-1)/gam);
    
    % Calculate change for while loop
    deltaT = (T03 - T03calc)/T03; %% Add relative term
    T03 = T03calc;
    T04 = T03;   
end

%% Repacking into structs
d.u1=u1; d.u4=u4; d.Um=Um; d.rpm=rpm; d.psi=psi;
q.M0=M0; q.M1=M1; q.M4=M4; q.Mndp4=Mndp4; q.Mndp01=Mndp01;
q.T00=T00; q.T01=T01; q.T03=T03; q.T03s=T03s; q.T04=T04;
q.T1=T1; q.T4=T4;
q.p00=p00; q.p01=p01; q.p03=p03; q.p04=p04; q.p4=p4;
