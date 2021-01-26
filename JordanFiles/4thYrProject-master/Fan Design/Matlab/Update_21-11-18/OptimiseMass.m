function [r_cCRIT, r_h, phi, psi, omega, OPThrust, sig] = OptimiseMass()
% Output r_c r_h phi and psi based on analyses of motor mass
% 

%% INPUTS YOU CAN CHANGE
global n

% ND Coefficients
phi         = 0.85;         % FLOW COEFFICIENT          
psi         = 0.25;         % STAGE LOADING

% Nominal Geometry

r_h         = 20e-3;        % HUB RADIUS
t           = 2e-3;         % CASING THICKNESS
r_c         = 56.6e-3;      % NOMINAL CASING RADIUS

turboL = 70e-3;
inletL = 20e-3;

% Motor Speed
rpm         = 7500;         % NOMINAL MOTOR SPEED

% Drone body mass
totalmass   = 1.5;          % DRONE EMPTY WEIGHT
motormass   = 0.06;

% Variables to plot
variables = 'phipsi';
mesh = 1;
chics = 1;
plott = 1;

%% Efficiencies
% Set, save in .mat

% Electrical
eta.motor       = 0.95;
eta.control     = 0.99;

% Aerodynamic
eta.blades      = 0.9;
eta.inlet       = 0.95;
eta.diffusor    = 0.95;

eta.propulsor   = eta.blades * eta.inlet * eta.diffusor;

eta.total       = eta.propulsor * eta.motor * eta.control;

%% Save
save('efficiencies', 'eta');

%% Show operating range for diffusion factor
phicontour = linspace(0, 1, 100);
psicontour = linspace(0, 0.5, 100);
[X,Y] = meshgrid(phicontour, psicontour);

desiredsig = 1.1;
sigmamesh = sqrt((X .^2) ./ (2*Y));
% figure(107); subplot(4,2,[2 4 6]); contourf(X,Y,sigmamesh, linspace(1,3,1000), 'LineColor', 'none'); xlabel("\phi"); ylabel("\psi"); title("Diffusion Factor (\sig)");
% hold on; plot([0, 1], [0.25,0.25]); plot([0.65, 0.65], [0,.5]); plot(phicontour, 0.5*phicontour.^2,'r'); colorbar;

%% Set fixed variables

% Test points
n           = 500;                          % Number of test r_c data points

% Flow conditions                         % ND
phirange = linspace(0.0, 1.0, n);
psirange = linspace(0.0, 0.5, n);

% Geometry
sig       = sqrt(phi^2/(2*psi));

r_m         = (r_c + r_h)/2;
r_cMAX      = 100e-3; 
tMAX        = 5e-3;
trange      = linspace(1e-3, tMAX, n);
r_hrange    = linspace(r_h, r_cMAX,n);
r_crange    = linspace(r_h, r_cMAX,n);% m
sigmarange  = linspace(1, 2, n);

% Blade speeds (mean)
rpmrange    = linspace(5000,10000,n);
omegarange  = rpmrange*pi/30;
omega       = rpm*pi/30;

% Densities
rhoPLA      = 1240;                         % kg/m^3
rhoAIR      = 1.225;                        % kg/m^3

% Drone mass and configuration
numberofprop= 4;

%% Diffuser and Inlet Volume for r_cCRITICAL
% Exit dimensions
r_ecrange = 0.5*(r_crange*(1+sig) + r_h*(1-sig));
r_ehrange = 0.5*(r_crange*(1-sig) + r_h*(1+sig));

r_icrange = r_crange + (r_crange-r_h).*0.25;
r_ihrange = r_h - (r_crange-r_h).*0.25;

L = (r_crange-r_h)/DF(phi,psi);

% Diffuser Casing
ch1 = L.*(r_ecrange + t) ./ (r_ecrange - r_crange);
ch2 = L.*(r_ecrange) ./ (r_ecrange - r_crange);
volcase = (pi/3)* (ch1.*(r_ecrange + t).^2 - (ch1-L).*(r_crange - t).^2 - ch2.*r_ecrange.^2 + (ch2-L).*r_crange.^2);

% Diffuser Hub 
hh1 = L.*(r_h) ./ (r_h - r_ehrange);
hh2 = L.*(r_h - t) ./ (r_h - r_ehrange);
volhub = (pi/3)* (hh1*r_h^2 - (hh1-L).*(r_ehrange).^2 - hh2*(r_h - t)^2 + (hh2-L).*(r_ehrange - t).^2);

% Inlet Casing
ich1 = inletL.*(r_icrange + t) ./ (r_icrange - r_crange);
ich2 = inletL.*(r_icrange) ./ (r_icrange - r_crange);
volicase = (pi/3)* (ich1.*(r_icrange + t).^2 - (ich1-inletL).*(r_crange - t).^2 - ich2.*r_icrange.^2 + (ich2-inletL).*r_crange.^2);

% Inlet Hub 
ihh1 = inletL.*(r_h) ./ (r_h - r_ihrange);
ihh2 = inletL.*(r_h - t) ./ (r_h - r_ihrange);
volihub = (pi/3)* (ihh1*r_h^2 - (ihh1-inletL).*(r_ihrange).^2 - ihh2*(r_h - t)^2 + (ihh2-inletL).*(r_ihrange - t).^2);

% Mass of propulsor, excluding blades (negligible)
casingmass = rhoPLA*pi*turboL*((r_crange + t).^2 - (r_crange).^2);
hubmass = rhoPLA*pi*turboL*(r_h^2 - (r_h - 10e-3)^2);
unitmass = casingmass + hubmass + motormass;

% Weight of diffuser
exitmass = rhoPLA*(volhub + volcase);

% Weight of Inlet
inletmass = rhoPLA*(volihub + volicase);

% Mass of propulsor
mass = inletmass + exitmass + unitmass;

%% r_cCRITICAL 
thrust = phi*sqrt(2*psi)*rhoAIR*pi*(r_crange.^2 - r_h^2).*(omega^2).*((r_crange + r_h)./2).^2;
power = rhoAIR*pi*phi*psi*(r_crange.^2 - r_h^2).*(omega^3).*((r_crange + r_h)./2).^3;

figureofmerit = (thrust./power) .* sqrt(thrust ./ (2*rhoAIR*pi*(r_crange.^2 - r_h^2)) );
    
% Thrust excess/deficit after mass taken off
deficit = thrust ./ ((mass + totalmass/numberofprop)*9.81);

% Find r_c for breakeven
x = linspace(r_h, r_cMAX, n);
y = ones(n,1);

for i = 1:n
    if deficit(i) > 1
        r_cCRIT = r_crange(i);
        OPThrust = thrust(i);
        OPPower = power(i);
        OPFigureOfMerit = figureofmerit(i);
        critidx = i;
        break
    elseif i == n
        disp("INSUFFICIENT THRUST TO LIFTOFF");
        r_cCRIT = r_cMAX;
    end
end 

disp("Thrust deficit positive when r_c =")
disp(r_cCRIT);

W.total = (mass(critidx) + totalmass/numberofprop)*9.81;
W.body = (totalmass/numberofprop)*9.81;
W.diff = exitmass(critidx)*9.81;
W.case = casingmass(critidx)*9.81;
W.hub = hubmass*9.81;
W.inlet = inletmass(critidx)*9.81;
W.motor = motormass*9.81;
W.numberofprop = numberofprop;

diffgeom.rc = r_ecrange(critidx);
diffgeom.rh = r_ehrange(critidx);
diffgeom.L = L(critidx);

diffgeom.range.rc = [r_ecrange(critidx) (r_ecrange(critidx)+t)];
diffgeom.range.rh = [(r_ehrange(critidx)-t) r_ehrange(critidx)];

inletgeom.rc = r_icrange(critidx);
inletgeom.rh = r_ihrange(critidx);
inletgeom.L = inletL;

inletgeom.range.rc = [r_icrange(critidx) (r_icrange(critidx)+t)];
inletgeom.range.rh = [(r_ihrange(critidx)-t) r_ihrange(critidx)];

geom.phi = phi;
geom.psi = psi;
geom.t = t;
geom.r_h = r_h;
geom.r_cMAX = r_cMAX;
geom.turboL = turboL;
geom.inletL = inletL;
geom.omega = omega;
geom.rhoPLA = rhoPLA;
geom.rhoAIR = rhoAIR;
geom.variables = variables;
geom.plott = plott;
geom.totalmass = totalmass;
geom.numberofprop = numberofprop;
geom.sig=sig;

save('meshvars', 'geom');
save('vars', 'phi', 'psi', 'r_h', 'r_cCRIT', 't', 'rpm', 'W', 'diffgeom', 'inletgeom');

%% Run MESH CALCULATIONS
if mesh == 1 
    OptimiseMassMESH;
end

%% Determine diffusion factor and length of diffuser

sig = phi/sqrt(2*psi);
L = (r_cCRIT - r_h)/DF(phi, psi);

Ax = pi*(r_cCRIT^2 - r_h^2);
Ae = Ax*sig;

r_ec = r_ecrange(critidx);
r_eh = r_ehrange(critidx);

r_ic = r_icrange(critidx);
r_ih = r_ihrange(critidx);

disp("Non-dimensional OP  [ phi psi ] ")
disp([ phi psi ]);

disp("Operating Point  [ Thrust Power FoM ] ")
disp([ OPThrust OPPower OPFigureOfMerit ]);

disp("Required diffuser  [ Length Sigma ] ")
disp([ L sig ]);

figure(107);
subplot(4,2,1); plot(r_crange, mass); xlabel("Casing Radius"); ylabel("Propulsor Mass / kg"); xlim([r_h r_cMAX]);
subplot(4,2,3); plot(r_crange, thrust); xlabel("Casing Radius"); ylabel("Thrust / N"); xlim([r_h r_cMAX]);
subplot(4,2,5); hold on; plot(r_crange, deficit); plot(x, y); xlabel("Casing Radius"); ylabel("T - Mg"); xlim([r_h r_cMAX]);


%% Print diffuser (SAME AS printDiff())

% Length of turbo machinery section 

figure(100); hold on;

% Inlet
% Inlet Internal
plot([0 -inletL], [r_h r_ih],'b'); plot([0 -inletL], [-r_h -r_ih],'b');         % Diffuser hub
plot([0 -inletL], [r_cCRIT r_ic],'b'); plot([0 -inletL], [-r_cCRIT -r_ic],'b'); % Diffuser casing

% Inlet External
plot([0 -inletL], [r_h-t r_ih-t],'b'); plot([0 -inletL], [-r_h+t -r_ih+t],'b');         % Diffuser hub
plot([0 -inletL], [r_cCRIT+t r_ic+t],'b'); plot([0 -inletL], [-r_cCRIT-t -r_ic-t],'b'); % Diffuser casing

% Turbomachinery
plot([0, turboL], [r_cCRIT r_cCRIT],'r'); plot([0, turboL], [-r_cCRIT -r_cCRIT],'r');           % Casing
plot([0, turboL], [r_cCRIT+t r_cCRIT+t],'r'); plot([0, turboL], [-r_cCRIT-t -r_cCRIT-t],'r');
plot([0, turboL], [r_h r_h],'r'); plot([0, turboL], [-r_h -r_h],'r');                           % Hub
plot([0, turboL], [r_h-t r_h-t],'r'); plot([0, turboL], [-r_h+t -r_h+t],'r');

% Diffuser Internal
plot([turboL (turboL+L)], [r_h r_eh],'b'); plot([turboL (turboL+L)], [-r_h -r_eh],'b');         % Diffuser hub
plot([turboL (turboL+L)], [r_cCRIT r_ec],'b'); plot([turboL (turboL+L)], [-r_cCRIT -r_ec],'b'); % Diffuser casing

% Diffuser External
plot([turboL (turboL+L)], [r_h-t r_eh-t],'b'); plot([turboL (turboL+L)], [-r_h+t -r_eh+t],'b');         % Diffuser hub
plot([turboL (turboL+L)], [r_cCRIT+t r_ec+t],'b'); plot([turboL (turboL+L)], [-r_cCRIT-t -r_ec-t],'b'); % Diffuser casing

%plot([turboL+L, turboL+L], [r_eh -r_eh],'b');                                                   % Hub cover
plot([0, turboL], [0 0],'r'); plot([turboL-10e-3, turboL], [10e-3 0],'r');plot([turboL-10e-3, turboL], [-10e-3 0],'r');

title('Target Diffuser Geometry');daspect([1 1 1]); xlim([-20e-3 0.3]); ylim([-0.1 0.1]);

fignamediff = ['DiffGeom_' num2str(phi) '-' num2str(psi) '_' num2str(sig) '.png'];
saveas(gcf,fignamediff);

figure(107); subplot(4,2,[2 4 6]); contourf(X,Y,sigmamesh, linspace(1,3,1000), 'LineColor', 'none'); xlabel("\phi"); ylabel("\psi"); title("Diffusion Factor (\sigma)");
hold on; plot([0, 1], [psi psi], 'k'); plot([phi phi], [0,.5], 'k'); plot(phicontour, 0.5*phicontour.^2,'r'); plot(phicontour, (0.5*phicontour.^2/((r_c + r_h)/(r_c - r_h))), 'r'); colorbar;

if chics == 1
    chicsatdesign;
end
end