R1 = load('free_dfconst_45.mat');
R2 = load('free_dfvary_25-45.mat');
R3 = load('mixed_dfconst_45.mat');
R4 = load('mixed_dfvary_25-45.mat');

rm = sqrt(mean([(20e-3^2),(60e-3^2)]));

cm1 = (40e-3)/R1.R.AR;
cm2 = (40e-3)/R2.R.AR;
cm3 = (40e-3)/R3.R.AR;
cm4 = (40e-3)/R4.R.AR;

N1 = 2*pi*rm / (R1.R.pitchchord*cm1);
N2 = 2*pi*rm / (R2.R.pitchchord*cm2);
N3 = 2*pi*rm / (R3.R.pitchchord*cm3);
N4 = 2*pi*rm / (R4.R.pitchchord*cm4);

%% Plot
figure(5);
hold on; grid on; box on; 
ylabel('% of Span', 'FontSize', 12); xlabel('c / mm', 'FontSize', 12);
plot(1000*R1.R.span.chord,R1.x, 'k');
plot(1000*R2.R.span.chord,R2.x, '-.k');
plot(1000*R3.R.span.chord,R3.x, 'b');
plot(1000*R4.R.span.chord,R4.x, '-.b');
% legend('Free DF= 0.45','Free 0.25 < DF < 0.45','Custom DF = 0.45','Custom 0.25 < DF < 0.45', 'Location', 'east');
set(gcf, 'Position', [0 0 300 200]);
% title('Chord'); 

figure(51);
hold on; grid on; box on; 
yticklabels(['','','','','','']);xlabel('s/c', 'FontSize', 12);yticks([0 0.2 0.4 0.6 0.8 1]);
plot(R1.R.span.pitchchord,R1.x, 'k');
plot(R2.R.span.pitchchord,R2.x, '-.k');
plot(R3.R.span.pitchchord,R3.x, 'b');
plot(R4.R.span.pitchchord,R4.x, '-.b');
legend('Free Vortex DF= 0.45','Free Vortex 0.25 < DF < 0.45','Mixed Vortex DF = 0.45','Mixed Vortex 0.25 < DF < 0.45', 'Location', 'southeast');
set(gcf, 'Position', [0 0 300 200]);
% title('Pitch-Chord');

figure(6);
hold on; grid on; box on; 
ylabel('% of Span', 'FontSize', 12); xlabel('\delta / Degrees', 'FontSize', 12);
plot(R1.R.span.delta,R1.x, 'k');
plot(R2.R.span.delta,R2.x, '-.k');
plot(R3.R.span.delta,R3.x, 'b');
plot(R4.R.span.delta,R4.x, '-.b');
% legend('Free DF= 0.45','Free 0.25 < DF < 0.45','Custom DF = 0.45','Custom 0.25 < DF < 0.45', 'Location', 'east');
set(gcf, 'Position', [0 0 300 200]);
% title('Deviation'); 

figure(61);
hold on; grid on; box on;
xlabel('\chi_2 / Degrees', 'FontSize', 12); yticklabels(['','','','','','']); yticks([0 0.2 0.4 0.6 0.8 1]);
plot(R1.ang.span.chi2,R1.x, 'k');
plot(R2.ang.span.chi2,R2.x, '-.k');
plot(R3.ang.span.chi2,R3.x, 'b');
plot(R4.ang.span.chi2,R4.x, '-.b');
legend('Free Vortex DF= 0.45','Free Vortex 0.25 < DF < 0.45','Mixed Vortex DF = 0.45','Mixed Vortex 0.25 < DF < 0.45', 'Location', 'southeast');
set(gcf, 'Position', [0 0 300 200]);
% title('Rotor Exit Metal Angle'); 
%% APC Figure of merit system efficiency calculation
load('APC_data.mat');
load('EXP_DATA_PROP.mat', 'APC');
figure(12); hold on; grid on; box on;
xlabel('RPM'); ylabel('M_F'); 
ylim([0, 0.8]); xlim([2500 7500]);

plot(APC.rpmmean, APC.FOM,'ro'); 
plot([min(APC_data.RPM(2:6)) max(APC_data.RPM(2:6))], [mean(APC.FOM) mean(APC.FOM)], 'r'); 
plot(APC_data.RPM(2:6), APC_data.FOM(2:6), 'bo'); 
plot([min(APC_data.RPM(2:6)) max(APC_data.RPM(2:6))], [mean(APC_data.FOM(2:6)) mean(APC_data.FOM(2:6))], 'b');

legend('Experiment', 'Experiment mean', 'Calibration data', 'Calibration data mean', 'Location', 'southeast');
set(gcf, 'Position', [0 0 380 250]);

APC.eff = mean(APC.FOM) / mean(APC_data.FOM(2:6));
sys_eff = APC.eff;
%% Plot Component WEIGHTS
load('PropMass.mat');
onshape.total = 0.358*4 + propmass.chassis;
onshape.diff = 0.072;

figure(101); b = bar([propmass.chassis (propmass.motor + propmass.blade) 0 0 0; 
    propmass.chassis (propmass.motor + propmass.blade) (propmass.dhub+propmass.case) propmass.annul propmass.inlet; 
    propmass.chassis 0 0 0 0;
    propmass.motor 0 0 0 0; 
    (propmass.dhub+propmass.case) 0 0 0 0; 
    propmass.annul 0 0 0 0; 
    propmass.inlet 0 0 0 0], 'stacked'); 

legend('Flying Test Bed', 'Motor', 'Exit Duct', 'Blade Duct', 'Inlet'); 
xticklabels({'Prop Total', 'Fan Total', 'Body','Motor', 'Exit Duct', 'Blade Duct', 'Inlet'});
grid on;

% xtips2 = b(1).XEndPoints(3:8); ytips2 = b(1).YEndPoints(3:8);
% labels2 = string(round(b(1).YData(3:8)));
% text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom');

b(1).FaceColor = 'flat';
b(1).CData(3,:) = [.5 .5 .5];b(1).CData(4,:) = [.5 .5 .5];b(1).CData(5,:) = [.5 .5 .5];b(1).CData(6,:) = [.5 .5 .5];b(1).CData(7,:) = [.5 .5 .5];
title('Component Masses'); ylabel('Mass / kg');

%% Plot Test Bed WEIGHTS
load('V3.mat');
fom_fan_ish = V3.P20.FOM;
phi = 0.8;
psi = 0.25;
sig = SIG(phi,psi);
[mass, ~] = MassModel_SIG(sig,V3.rc,V3.rh,0);
% prop_hover_thrust = mass.prop*9.81/4;
% prop_hover_power = (prop_hover_thrust / mean(APC_data.FOM(2:6))) * sqrt(prop_hover_thrust/(2*1.225*pi*(APC_data.D/2)^2));
% fan_hover_thrust = ((fom_fan_ish/APC.eff) * prop_hover_power * sqrt(2*1.225*0.0092))^(2/3);
% teestmass = fan_hover_thrust*4/9.81;

% load('PropMass.mat');
onshape.total = 0.358*4 + mass.chassis;
onshape.diff = 0.072;

figure(101); b = bar([mass.chassis+4*mass.motor+4*12e-3 4*mass.blade 4*(mass.dhub+mass.case) 4*mass.annul 4*mass.inlet; 
    onshape.total 0 0 0 0;
    mass.prop 0 0 0 0], 'stacked'); 

legend('Flying Test Bed and Motor', 'Blades', 'Exit Duct', 'Blade Duct', 'Inlet', 'Location', 'northeast'); 
xticklabels({'Fan - Mass Model', 'Fan - Actual', 'Propeller'});
ylabel('Mass / kg'); ylim([0 3]);
grid on;

xtips1 = b(1).XEndPoints(1); ytips1 = b(5).YEndPoints(1);
xtips2 = b(5).XEndPoints(2); ytips2 = b(1).YEndPoints(2);
xtips3 = b(1).XEndPoints(3); ytips3 = b(1).YEndPoints(3);

labels1 = string(round(sum([b(1).YData(1) b(2).YData(1) b(3).YData(1) b(4).YData(1) b(5).YData(1)]), 2));
labels2 = string(round(sum([b(1).YData(2) b(2).YData(2)]),2));
labels3 = string(round(b(1).YData(3),2));

text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom');
text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom');
text(xtips3,ytips3,labels3,'HorizontalAlignment','center','VerticalAlignment','bottom');

b(1).FaceColor = 'flat';
b(1).CData(2,:) = [.8 .8 .8]; b(1).CData(3,:) = [.8 .8 .8];

set(gcf, 'Position', [0 0 600 200]);

%% Aerodynamic Efficiencies
load('V2.mat');
load('V3.mat');
load('APC_data.mat');
sys_eff = 0.4477;

maxMF.v2_1 = sqrt(2*V2.P12.sig);
maxMF.v2_2 = sqrt(2*V2.P20.sig);
maxMF.v3_1 = sqrt(2*V3.P20.sig);

aero_eff.v2_1 = (V2.P12.FOM/sys_eff)/maxMF.v2_1;
aero_eff.v2_2 = (V2.P20.FOM/sys_eff)/maxMF.v2_2;
aero_eff.v3_1 = (V3.P20.FOM/sys_eff)/maxMF.v3_1;

load('EXP_DATA_PROP.mat', 'APC');
figure(13); hold on; grid on; box on;
xlabel('RPM', 'FontSize', 12); ylabel('M_{F,s}', 'FontSize', 12); 
ylim([0, 1.6]); xlim([2500 7500]);

plot(V2.P12.rpmmean, V2.P12.FOM/sys_eff, 'xr');
plot(V2.P20.rpmmean, V2.P20.FOM/sys_eff, 'or');
plot(6500, V3.P20.FOM/sys_eff, '*r');
plot([3000 7000], [maxMF.v2_1 maxMF.v2_1], 'r');

plot(APC.rpmmean, APC.FOM./sys_eff, 'ko'); 
plot([3000 7000], [mean(APC.FOM(2:6))./sys_eff mean(APC.FOM(2:6))./sys_eff], 'k');
plot([3000 7000], [1 1], '--k');

legend('Fan A', 'Fan B', 'Fan C', 'Ducted fan design', 'Prop. exp.', 'Prop. mean', 'Prop. max', 'Location', 'southwest');
set(gcf, 'Position', [0 0 300 300]);
set(gcf, 'Position', [0 0 500 300]);

% Bar chart of M_F
figure(101); hold on; 

plot([0.5 3.5], [1.51 1.51], '-.b', 'LineWidth', 2);
plot([3.5 4.5], [1 1], '-.k', 'LineWidth', 2);

b = bar([V2.P12.FOM/sys_eff; V2.P20.FOM/sys_eff; V3.P20.FOM/sys_eff; mean(APC_data.FOM(2:6))]); 

% legend('Flying Test Bed and Motor', 'Blades', 'Exit Duct', 'Blade Duct', 'Inlet', 'Location', 'northeast'); 
xticks([1,2,3,4]);
xticklabels({'Fan A', 'Fan B', 'Fan C', 'Propeller'});
ylabel('$M_{F,s}$', 'Interpreter', 'Latex', 'FontSize', 15); ylim([0 1.75]);
grid on; box on;

xtips1 = b(1).XEndPoints(1); ytips1 = b(1).YEndPoints(1);
xtips2 = b(1).XEndPoints(2); ytips2 = b(1).YEndPoints(2);
xtips3 = b(1).XEndPoints(3); ytips3 = b(1).YEndPoints(3);
xtips4 = b(1).XEndPoints(4); ytips4 = b(1).YEndPoints(4);

labels1 = string(round(b(1).YData(1), 2));
labels2 = string(round(b(1).YData(2), 2));
labels3 = string(round(b(1).YData(3), 2));
labels4 = string(round(b(1).YData(4), 2));

text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom');
text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom');
text(xtips3,ytips3,labels3,'HorizontalAlignment','center','VerticalAlignment','bottom');
text(xtips4,ytips4,labels4,'HorizontalAlignment','center','VerticalAlignment','bottom');

b(1).FaceColor = 'flat';
b(1).CData(1,:) = [0 0 1];
b(1).CData(2,:) = [0 0 1];
b(1).CData(3,:) = [0 0 1];
b(1).CData(4,:) = [.8 .8 .8];

legend('Ducted Fan Design M_F', 'Propeller Maximum M_F', 'Location', 'northeast');

set(gcf, 'Position', [0 0 700 200]);
%% Shroud clearance flows
b1 = -59.34;
b2 = -55.05;
V2 = 37.6; % V2 in paper is V rel1 for the rotor reference frame
U = 29.17;
psi = 0.25;
cd = 0.5;

g = 2e-3;
h = 40e-3;

ml_mm = (g/h)*cd*sqrt(tand(b1)^2 - tand(b2)^2);
shroud_loss_coeff = 2*ml_mm*(1 - (tand(b2)/tand(b1))*sind(b1)^2) * (0.5*V2^2) / (psi*U^2);

% %% Reynolds number across span
% figure; grid on; box on;
% xlabel('Re', 'FontSize', 12);
% ylabel('Span', 'FontSize', 12);
% xlim([1.5e4 5.5e4]);
% 
% plot(Re.span, x, 'k');
% 
% set(gcf, 'Position', [0 0 250 250]);
% 
% % Casing flow coefficient
% phi_exp = V2.P20.Vx / (V2.P20.rpmmean * (pi / 30) * V2.rc);
% 
% % Having run 'Run'
% Vrel1 = sqrt((V.span.x.*((phi_exp/phi.sec(end))*(V2.P20.rpmmean/rpm))).^2 + (V2.P20.rpmmean*(pi/30).*radius).^2);
% Respan_exp = Vrel1.*R.span.chord / (1.48e-5);
% hold on;
% plot(Respan_exp, x, 'b');
% legend('Design', 'Experiment', 'Location', 'southeast');

%% Design Space
% close 24;
n = 500;
mo = imread('maffioli_overlay.png');
[M,N,~] = size(mo);
figure(24); hold on;
im = image('CData',mo,'XData',[0.3055 0.8693],'YData',[0.5186 0.1266]);
set(im, 'AlphaData', 0.5*ones(M,N));
box on; grid on; axis([0 1 0 0.5]);
phi=linspace(0,1,n);
psi=linspace(0,0.5,n);

plot(phi, (phi.^2)/2, 'k');
plot(phi, (phi.^2)/8, '-.k');

% Plot WF/WP in terms of phi psi rc rh
[phi, psi] = meshgrid(phi,psi);
area_fan = pi*(60e-3^2 - 20e-3^2);
area_prop = pi*(0.254/2)^2;
sig = SIG(phi, psi);
% [W,~] = MassModel_SIG(sig, 60e-3, 20e-3, 5);
% 
% WF_WP = ( (phi.^2 .* pi .* area_fan) ./ (2 .* psi .* area_prop) ).^(1/3);
% 
% cont = WF_WP ./ (W.total ./ W.prop);

% contourf(phi, psi, cont, [1 2]);
% [h,C1] = contour(phi, psi, WF_WP); clabel(h,C1);
legend('\sigma = 1', 'Symmetric diffuser limit', 'Location', 'northwest');
set(gcf, 'Position', [500 300 300 300]);
xlabel('$\phi_m$', 'Interpreter', 'Latex', 'FontSize',15);ylabel('$\psi_m$', 'Interpreter', 'Latex', 'FontSize',15);
c = colorbar;
c.Ticks = [0 1/7 2/7 3/7 4/7 5/7 6/7 1];
c.TickLabels = [0.8 0.83 0.87 0.9 0.91 0.92 0.93 0.94];
map = [0 0 1; 1 1 1; 1 0 0];
colormap(redblue(50));
c.Label.String = '\eta_a';

%% Plot LHS and RHS as functions of rc
n = 500;
sig_rhs = sqrt(0.8^2 ./ (2*0.25));
% sig_rhs = 1.0;
rc_rhs = linspace(15e-3, 200e-3, n);
% rc_rhs = 100e-3;

% rh_rhs = 20e-3; % Fixed hub
rh_rhs = linspace(15e-3, 200e-3, n); % Fixed htr

[rc_rhs, rh_rhs] = meshgrid(rc_rhs, rh_rhs);

area_rhs = pi.*(rc_rhs.^2 - rh_rhs.^2);
rc_rhs(area_rhs < 0) = NaN;
rh_rhs(area_rhs < 0) = NaN;
area_rhs(area_rhs < 0) = NaN;

HTR_rhs = rh_rhs ./ rc_rhs;
% 
% [mass_rc, ~] = MassModel_SIG(sig_rhs, rc_rhs, rh_rhs, 15);
% LHS = mass_rc.total ./ mass_rc.prop;
% RHS = ( (2.*sig_rhs.*area_rhs) ./ (pi*(0.254/2)^2) ).^(1/3);
% 
% [mass_rc_p1, ~] = MassModel_SIG(sig_rhs, rc_rhs, rh_rhs, 2);
% LHS_p1 = mass_rc_p1.total ./ mass_rc_p1.prop;
% RHS_p1 = ( (2.*sig_rhs.*area_rhs) ./ (pi*(0.254/2)^2) ).^(1/3);
% 
% [mass_rc_p2, ~] = MassModel_SIG(sig_rhs, rc_rhs, rh_rhs, 5);
% LHS_p2 = mass_rc_p2.total ./ mass_rc_p2.prop;
% RHS_p2 = ( (2.*sig_rhs.*area_rhs) ./ (pi*(0.254/2)^2) ).^(1/3);
% 
% [mass_rc_p3, ~] = MassModel_SIG(sig_rhs, rc_rhs, rh_rhs, 10);
% LHS_p3 = mass_rc_p3.total ./ mass_rc_p3.prop;
% RHS_p3 = ( (2.*sig_rhs.*area_rhs) ./ (pi*(0.254/2)^2) ).^(1/3);

figure(25); hold on; box on; grid on;
% plot([0.5 2.5], [0.5 2.5], 'k');
% plot(RHS,LHS,'r');
% plot(RHS_p1,LHS_p1,'-.b');
% plot(RHS_p2,LHS_p2,'--b');
% plot(RHS_p3,LHS_p3,':b');

% contourf(1000*rc_rhs, HTR_rhs, (RHS - LHS));


xlabel('$A_{x,fan} / A_{x,prop.}$', 'Interpreter', 'Latex', 'FontSize', 12);
ylabel('$HTR$', 'Interpreter', 'Latex', 'FontSize', 12);

% xlabel('$\big[ 2\sigma \frac{A_F}{A_P} \big] ^{\frac{1}{3}}$', 'Interpreter', 'Latex', 'FontSize', 12);
% ylabel('$\frac{W_F}{W_P}$', 'Interpreter', 'Latex', 'FontSize', 12);
% legend('LHS = RHS', 'FTB payload', '2kg payload', '5kg payload', '10kg payload', 'Location', 'northwest');
set(gcf, 'Position', [500 300 600 300]);

[mass_rc_fw, diff_rc_fw] = MassModel_SIG(sig_rhs, rc_rhs, rh_rhs, 10);
LHS = mass_rc_fw.total./ mass_rc_fw.prop;
RHS = ( (((1.36/0.67).^2) * area_rhs) ./ (1*pi*(0.254/2)^2) ).^(1/3);%( (2*sig_rhs*area_rhs) ./ (1*pi*(0.254/2)^2) ).^(1/3);

omega = sqrt(2*sig_rhs.*(mass_rc_fw.total*9.81/4) ./ (1.225*pi*0.8^2*rc_rhs.^4.*(1-HTR_rhs.^4)) );

po = max(max(RHS - LHS));
LHS_min = LHS((RHS - LHS) == po);
rc_max = rc_rhs((RHS - LHS) == po);
htr_max = HTR_rhs((RHS - LHS) == po);
contourf((area_rhs./((pi*(0.254/2)^2))), HTR_rhs, (RHS - LHS));

% plot(1000*rc_max, htr_max, 'xk');
c = colorbar;
c.Label.String = '\Sigma';
colormap(redblue(50));
plot([0 0], [0 0], 'k');
plot([1 1], [0 1], '-.k','LineWidth',2);
[C,h] = contour((area_rhs./((pi*(0.254/2)^2))), HTR_rhs, omega*30/pi, [2000 3000 4000 6000  10000 20000 50000] ); clabel(C,h);
caxis([-1 1]);
legend('\Sigma', 'RPM', 'Prop. area', 'Location', 'northeast');

%% Plot LHS and RHS as functions of rc
n = 500;

sig_rhs = sqrt(0.8^2 ./ (2*0.25));
% sig_rhs = 1.0;
rc_rhs = linspace(15e-3, 200e-3, n);
% rc_rhs = 100e-3;

payload = linspace(.5,15,n);
htr = 0.2;

[rc_rhs, payload] = meshgrid(rc_rhs, payload);

rh_rhs = rc_rhs .* htr;
area_rhs = pi.*rc_rhs.^2 * (1 - htr^2);

rc_rhs(area_rhs < 0) = NaN;
rh_rhs(area_rhs < 0) = NaN;
area_rhs(area_rhs < 0) = NaN;

figure(25); hold on; box on; grid on;
xlabel('$\Lambda = A_{x,fan} / A_{x,prop.}$', 'Interpreter', 'Latex', 'FontSize', 12);
xlabel('$A_{x,fan} / A_{x,prop.}$', 'Interpreter', 'Latex', 'FontSize', 15);
% ylabel('$\frac{W_P}{g} / kg$', 'Interpreter', 'Latex', 'FontSize', 12);
% ylabel('$W / kg$', 'Interpreter', 'Latex', 'FontSize', 15);
set(gcf, 'Position', [500 300 300 300]);

for a =1:n
    for b =1:n
        [mass_p, ~] = MassModel_SIG(sig_rhs, rc_rhs(a,b), rh_rhs(a,b), payload(a,b));
        LHS(a,b) = mass_p.CRDF./ mass_p.prop;
        RHS(a,b) = ( ((1.36/0.67)^2 * area_rhs(a,b)) ./ (1*pi*(0.254/2)^2) ).^(1/3);%( (2*sig_rhs*area_rhs) ./ (1*pi*(0.254/2)^2) ).^(1/3);
        omega(a,b) = sqrt(2*sig_rhs.*(mass_p.total*9.81/4) ./ (1.225*pi*0.8^2*rc_rhs(a,b).^4.*(1-0.2.^4)) );
    end
end

po = max(max(RHS - LHS));
LHS_min = LHS((RHS - LHS) == po);
rc_max = rc_rhs((RHS - LHS) == po);
htr_max = HTR_rhs((RHS - LHS) == po);

contourf((area_rhs./((pi*(0.254/2)^2))), payload, (RHS - LHS), [-10 -5 -3 -2 -1 -0.5 -0.25 0 0.25 0.5 1], 'LineColor','None');
c = colorbar;
c.Label.String = '\Sigma';
colormap(redblue(50));
plot([0 0], [0 0], 'k');
% plot([min(min((area_rhs./((pi*(0.254/2)^2))))) max(max((area_rhs./((pi*(0.254/2)^2)))))], [1 1], '-.k','LineWidth',2);
% [C,h] = contour((area_rhs./((pi*(0.254/2)^2))), payload, omega*30/pi, [1500 2000 3000 4000 6000  10000 20000 50000] ); clabel(C,h);
[C,h] = contour((area_rhs./((pi*(0.254/2)^2))), payload, (RHS - LHS), [-10 -5 -3 -2 -1 -0.5 -0.25 0 0.25 0.5 1] ); clabel(C,h);
caxis([-1 1]);
% legend('\Sigma', 'RPM', 'Location', 'southeast');
legend('\Sigma', 'Location', 'southeast');
ylim([0.5,15]);
xlim([0.014 2.381]);

%% Operating speed
rho = 1.225;
phi = 0.8;
psi = 0.25;
sig = SIG(phi,psi);
rc = 60e-3;
rh = 20e-3;
[mass,diff] = MassModel_SIG(sig,rc,rh,0);

omega = sqrt(2*sig*(mass.total*9.81/8) ./ (rho*pi*(phi^2)*(rc^4 - rh^4)));
rpm = omega*30/pi;

power = pi*rho*phi^3*omega^3*(rc^2 - rh^2)*((rc^2 + rh^2)/2)^(3/2)/(2*sig^2);
torque = power/omega;

LHS = mass.total ./ mass.prop;
LHS_CRDF = mass.CRDF ./ mass.prop;

fan_fom = V3.P20.FOM / sys_eff;
prop_fom = mean(APC_data.FOM(2:6));

RHS = ( (fan_fom^2.*area_rhs) ./ (prop_fom*pi*(0.254/2)^2) ).^(1/3);