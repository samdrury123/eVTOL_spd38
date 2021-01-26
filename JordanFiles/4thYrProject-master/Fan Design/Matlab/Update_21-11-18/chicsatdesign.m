% Output characteristics of the vehicle at specified design point 
% [] = chisatdesign(r_cCRIT, r_h, t, phi, psi, rpm, W)
etasave;
load('vars');
load('efficiencies');
n = 1000; 

% Densities
rhoPLA      = 1240;                         % kg/m^3
rhoAIR      = 1.225;                        % kg/m^3

% Electrical properties
kv          = 1100;

% Values
Ax = pi*(r_cCRIT^2 - r_h^2);
U = rpm*(pi/30)*(r_cCRIT + r_h)/2;
sigma = sqrt((phi^2)/(2*psi));

% Payload
payload = 0;    % kilograms
weighttotal = W.total + (payload/W.numberofprop)*9.81;     % Newtons

%% Thrust and Power
T = phi*sqrt(2*psi)*rhoAIR*Ax*U^2;

% Power wrt thrust
Pt = (T^(3/2))/sqrt(4*sigma*rhoAIR*Ax);

% Power wrt coefficeints
% P = phi*psi*rhoAIR*Ax*U^3;

% Electrical at OP
I = Pt*kv/(rpm*eta.total);
V = rpm/kv;

% Propulsor operating point
OP.current = I; OP.voltage = V; OP.elecpower = Pt/eta.total; OP.flowpower = Pt; 
OP.thrust = T; OP.phi = phi; OP.psi = psi; OP.mass = weighttotal/9.81; 
OP.rpm = rpm; OP.sigma = sigma; OP.rc = r_cCRIT; OP.rh = r_h; OP.diff = diffgeom;

% Total operating point
OPtotal.current = OP.current*W.numberofprop; OPtotal.voltage = V; 
OPtotal.elecpower = OP.elecpower*W.numberofprop; OPtotal.flowpower = W.numberofprop*Pt; 
OPtotal.thrust = W.numberofprop*T; OPtotal.phi = phi; OPtotal.psi = psi; 
OPtotal.mass = W.numberofprop*weighttotal/9.81; OPtotal.rpm = rpm; OPtotal.sigma = sigma;
OPtotal.rc = r_cCRIT; OPtotal.rh = r_h; OPtotal.diff = diffgeom;

OP
OPtotal
diffgeom

save('OP', 'OP', 'OPtotal', 'diffgeom');
% message = [
%     'phi = ' num2str(OP.phi) ... 
%     ';    psi = ' num2str(OP.psi) ... 
%     ';    RPM = ' num2str(OP.rpm) ... 
%     ';    Casing Radius = ' num2str(1000*OP.rc) 'mm'];
% f = msgbox(message,'Operating Point');
% 
% messagediff = [
%     'Exit casing r = ' num2str(1000*diffgeom.rc) 'mm'... 
%     ';    Exit hub r = ' num2str(1000*diffgeom.rh) 'mm'... 
%     ';    Length = ' num2str(1000*diffgeom.L) 'mm'... 
%     ';    Sigma = ' num2str(sigma)];
% f = msgbox(messagediff,'Diffuser Geometry');
%% Variables
% Speed
rpmrange = linspace(5000, 10000, n);
Urange = rpmrange*(pi/30)*(r_cCRIT + r_h)/2;

% Thrust
Trange = phi*sqrt(2*psi)*rhoAIR*Ax*Urange.^2;

% Power
Prange = (Trange.^(3/2))/sqrt(4*sigma*rhoAIR*Ax);

% ND Thrust Excess
Terange = Trange ./ weighttotal;

% Current draw
Irange = Prange.*kv./(rpmrange*eta.total);
%% Plot WEIGHTS
figure(101); b = bar([W.body/(9.81/1000) W.motor/(9.81/1000) 0 0 0 0; W.body/(9.81/1000) W.motor/(9.81/1000) W.diff/(9.81/1000) W.case/(9.81/1000) W.hub/(9.81/1000) W.inlet/(9.81/1000); W.body/(9.81/1000) 0 0 0 0 0;W.motor/(9.81/1000) 0 0 0 0 0; W.diff/(9.81/1000) 0 0 0 0 0; W.case/(9.81/1000) 0 0 0 0 0; W.hub/(9.81/1000) 0 0 0 0 0; W.inlet/(9.81/1000) 0 0 0 0 0], 'stacked'); legend('Body', 'Motor', 'Diffuser', 'Casing', 'Hub', 'Inlet'); xticklabels({'Prop Total', 'Fan Total', 'Body','Motor', 'Diffuser', 'Casing', 'Hub', 'Inlet'});grid on;

% xtips2 = b(1).XEndPoints(3:8); ytips2 = b(1).YEndPoints(3:8);
% labels2 = string(round(b(1).YData(3:8)));
% text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom');

b(1).FaceColor = 'flat';
b(1).CData(3,:) = [.5 .5 .5];b(1).CData(4,:) = [.5 .5 .5];b(1).CData(5,:) = [.5 .5 .5];b(1).CData(6,:) = [.5 .5 .5];b(1).CData(7,:) = [.5 .5 .5];b(1).CData(8,:) = [.5 .5 .5];
title('Component Masses'); ylabel('Mass / g');

%% Plot ch'ics
figure(102); hold on; plot(rpmrange, Trange, 'b'); plot([5000 rpm], [weighttotal weighttotal],'r'); plot([rpm rpm], [0 weighttotal],'r'); title('RPM vs Thrust'); xlabel('RPM'); ylabel('Thrust');
figure(103); hold on; plot(rpmrange, Prange, 'b'); plot([5000 rpm],[Pt Pt], 'r'); plot([rpm rpm], [0 Pt],'r');title('RPM vs Power'); xlabel('RPM'); ylabel('Power');
figure(104); hold on; plot(rpmrange, Terange, 'b'); plot([5000 rpm], [1 1],'r'); plot([rpm rpm], [0 1],'r');title('RPM vs Thrust Excess'); xlabel('RPM'); ylabel('Thrust Excess');
figure(105); hold on; plot(rpmrange, Irange, 'b'); plot([5000 rpm], [I I],'r'); plot([rpm rpm], [0 I],'r');title('RPM vs Motor Current'); xlabel('RPM'); ylabel('Motor Current');
figure(106); hold on; plot(Trange, Prange, 'b'); plot([weighttotal weighttotal],[0 Pt], 'r'); plot([0 weighttotal], [Pt Pt],'r');title('Thrust vs Power'); xlabel('Thrust'); ylabel('Power');