%% ---- FUNCTION: CONTROL VOLUME ANALYSIS FOR DUCTED FAN ---- %%
function [d,g,eta] = CVanalysis(d,g,eta)
% Function using continuity, SFME and SFEE to analyse ducted fan
% Inputs: 
% - atm: atmospheric conditions
% - g: geometry of fan
% - eta: aerodynamic and electric system efficiencies
% - d.T: required thrust
% - d.u0: required forward speed
% - d.phi: flow coefficient of fan, can be an array
% - d.sig: area ratio of exit nozzle, can be an array
% Outputs: 
% - P: Ideal, shaft and electrical power in kW
% - vel: upstream, inlet and exit velocity, and blade velocities
% - rpm: shaft speed
% - psi: non-dimensional pressure rise, total-total and total-to-static
% - load: disk loading and thrust to power

% Unpack input structures and calculate secondary inputs
T = d.T;
u0 = d.u0;
phi = d.phi;
sig = d.sig;
ro = d.atm.ro;
rcas = g.rcas;
rhub = g.rcas * g.htrat;
A1 = pi*(rcas^2 - rhub^2);
% Define rmid so that equal mass flow above and below line:
rmid = ((rcas^2+rhub^2)/2)^0.5;

% Inlet velocity from SFME and quadratic equation. Eq (2) in .ppt
u1 = (u0 + (u0^2 + 4*T./(sig*ro*A1)).^0.5) ./ (2./sig);
u0 = u0;
u4 = u1 ./ sig;

% Powers calculated from SFEE and efficiencies. Eq (3),(4),(5) in .ppt
Pideal = 0.5*ro*A1.*u1.*((u1.^2./sig.^2) - u0^2);
Pshaft = Pideal ./ eta.aero;
Pelec = Pshaft ./ eta.elec;

% RPM and blade speeds calculated from u1, flow coefficient and geometry for hover
% % Eqs (6) and (7) in .ppt
% u1hov = ((sig.*T)./(ro*A1)).^5;
% 
% Umidhov = u1hov ./ phi;
% psitthov = phi.^2 ./ (2.*sig.^2);
% % Calculate C=tan(alpha2r) to estimate shape of ch'ic for off-design behaviour. Eq (17) in .ppt
% C = (psitthov-1) ./ phi;
% 
% % Calculate U if there is forward flight speed. Eq (19) in .ppt
% Umid = (-C .* u1 + ((C .* u1).^2 - 2*(u0.^2 - u4.^2)).^0.5)/2;
Umid = u1 ./ phi; 
Utip = Umid * rcas / rmid;
d.rpm =  Umid / (2*pi*rmid) * 60;
% Calculate flow coefficient for cruise. Eq (20) in .ppt
% d.phicru = u1 ./ Umid;

% Stage loading coefficient from inlet and exit pressure boundary
% conditions. Eqs (15) and (11) in .ppt. 
psitt = (u1.^2 - sig.^2 .* u0.^2)./(2 .* sig.^2 .* Umid.^2);
psits = 2*psitt - phi.^2;

% Disk Loading and Thrust to Power
TtoA = T / A1; 
TtoP = T ./ Pshaft;  

% Figure of merit for hover. Eq (13) in .ppt
eta.mf = sqrt(2*sig);

% Propulsive efficiency for cruise. Eq (14) in .ppt
eta.prop = 2 * u0 ./(u0 + u4);

% Torque
d.Tq = Pideal ./ (d.rpm / 60 * 2*pi);

% Mach Numbers (based on atmospheric temperature)
a = sqrt(1.4*287*d.atm.T);
Mtip = Utip ./ a;
Mjet = u4 ./ a;

% Sort variables into output structures
g.rhub = rhub; g.rmid = rmid; g.A1 = A1;
d.u0 = u0; d.u1 = u1; d.u4 = u4; d.Umid = Umid; d.Utip = Utip;
d.Pideal = Pideal / 1000; d.Pshaft = Pshaft / 1000; d.Pelec = Pelec / 1000;
d.psitt = psitt; d.psits = psits;
d.TtoA = TtoA; d.TtoP = TtoP * 1000;
d.Mtip = Mtip;
d.Mjet = Mjet;
d.eta_prop = eta.prop;

end