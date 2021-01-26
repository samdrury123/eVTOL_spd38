% Required paylod for fan to be preferred method (lower absolute power use)

% Number of test points
n = 100;

% Densities
rhoAIR      = 1.225;                        % kg/m^3

% Payload
payload = 0;    % kilograms

massrange = linspace(1, 100, n);

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

motormass   = 0.06;

% Drone body mass
totalmass   = massrange(i);          % DRONE EMPTY WEIGHT

% Flow conditions                         % ND
phirange = linspace(0.0, 1.0, n);
psirange = linspace(0.0, 0.5, n);

% Geometry
sig       = sqrt(phi^2/(2*psi));

r_m         = (r_c + r_h)/2;
r_cMAX      = 1000e-3; 
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
deficit = thrust ./ ((mass + totalmass./numberofprop)*9.81);

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

if OUTPUT == 1
    disp("Thrust deficit positive when r_c =")
    disp(r_cCRIT);
end

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

    % Values
    Ax = pi*(r_cCRIT^2 - r_h^2);
    U = rpm*(pi/30)*(r_cCRIT + r_h)/2;
    sigma = sqrt((phi^2)/(2*psi));

    weighttotal = W.total + (payload/W.numberofprop)*9.81;     % Newtons
for i =1:n
    %% Thrust and Power
    T = phi*sqrt(2*psi)*rhoAIR*Ax*U^2;

    % Power wrt thrust
    fanpower(i) = (T^(3/2))/sqrt(4*sigma*rhoAIR*Ax);
    
    [proprpm, propweight, propp] = proprpm();
    
    proppower(i) = propp;
    
    if fanpower(i) < proppower(i)
        disp('Limit of Fan Superiority');
        disp(massrange(i));
        break
    end
end

