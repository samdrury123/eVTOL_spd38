%% Plot Graphs (LONG INTAKE)
% Thrust vs RPM
% Power vs RPM
% FOM vs Re
% Pressure vs Axial Position
% Pressure Rise vs Flow Coefficient

%% Import
intersectrpmpwrAPC;

load('EXP_DATA.mat');
load('EXP_DATA_PROP.mat');

%% Calculations

S08 = calcs(S08);S08S = calcs(S08S);S08L = calcs(S08L);
S10 = calcs(S10);S10S = calcs(S10S);S10L = calcs(S10L);
S12 = calcs(S12);S12S = calcs(S12S);S12L = calcs(S12L);
APC = calcs(APC);BASE = calcs(BASE);

%% Thrust vs RPM
figure(2); subplot(2,3,1); hold on; title('Thrust vs RPM'); xlabel('RPM'); 
ylabel('Thrust / N'); ylim([0 max([max(S08L.T), max(S10L.T), max(S12L.T), max(APC.T), max(S12.T)])]);
% plot(S08.rpmmean, S08.T, '-or');
% plot(S10.rpmmean, S10.T, '-ok');
plot(S12S.rpmmean, S12S.T, 'xb');
plot(S08L.rpmmean, S08L.T, '-xr');
plot(S10L.rpmmean, S10L.T, '-xk');
plot(S12L.rpmmean, S12L.T, '-xb');
plot(APC.rpmmean, APC.T, 'xg');
legend('\sigma = 0.8', '\sigma = 1.0', '\sigma = 1.2', '\sigma = 1.2 (LONG)', 'APC Propellor', 'Location', 'southeast');
% legend('SHORT Sigma = 0.8', 'SHORT Sigma = 1.0', 'SHORT Sigma = 1.2', 'LONG Sigma = 0.8', 'LONG Sigma = 1.0', 'LONG Sigma = 1.2', 'Location', 'southeast');

%% Power vs RPM
figure(2); subplot(2,3,4); hold on; title('Power vs RPM'); xlabel('RPM'); 
ylabel('Power / W'); ylim([0 max([max(S08L.P), max(S10L.P), max(S12L.P), max(APC.P), max(S12.P)])]);
% plot(S08.rpmmean, S08.P, '-or');
% plot(S10.rpmmean, S10.P, '-ok');

plot(S08L.rpmmean, S08L.P, '-xr');
plot(S10L.rpmmean, S10L.P, '-xk');
plot(S12L.rpmmean, S12L.P, '-xb');
plot(S12S.rpmmean, S12S.P, 'xb');
plot(APC.rpmmean, APC.P, 'xg');
legend('\sigma = 0.8', '\sigma = 1.0', '\sigma = 1.2', '\sigma = 1.2 (LONG)', 'APC Propellor', 'Location', 'southeast');
% legend('SHORT Sigma = 0.8', 'SHORT Sigma = 1.0', 'SHORT Sigma = 1.2', 'LONG Sigma = 0.8', 'LONG Sigma = 1.0', 'LONG Sigma = 1.2', 'Location', 'southeast');

%% FOM vs Re
figure(2); subplot(2,3,[2 5]); hold on; title('FoM vs Re'); xlabel('Re'); 
ylabel('FoM / M_f'); ylim([0 1]);
% plot(S08.rpmmean, S08.FOM, '-or');
% plot(S10.rpmmean, S10.FOM, '-ok');
plot(S08L.Re, S08L.FOM, '-xr');
plot(S10L.Re, S10L.FOM, '-xk');
plot(S12L.Re, S12L.FOM, '-xb');
plot(S12S.Re, S12S.FOM, 'xb');
plot([min(S08L.Re) max(S12.Re)], [mean(S12S.FOM) mean(S12S.FOM)], 'b');
plot([min(S08L.Re) max(S12.Re)], [mean(BASE.FOM) mean(BASE.FOM)], '-.g');
plot([min(S08L.Re) max(S12.Re)], [mean(APC.FOM) mean(APC.FOM)], '-g');
% plot(APC.Re, APC.FOM, 'xg');
legend('\sigma = 0.8', '\sigma = 1.0', '\sigma = 1.2', '\sigma = 1.2 (LONG)', '\sigma = 1.2 (MEAN LONG)', 'Baseline Propellor', 'APC Propellor', 'Location', 'northeast');
% legend('SHORT Sigma = 0.8', 'SHORT Sigma = 1.0', 'SHORT Sigma = 1.2', 'LONG Sigma = 0.8', 'LONG Sigma = 1.0', 'LONG Sigma = 1.2', 'Baseline Propellor', 'APC Propellor', 'Location', 'east');

%% Thrust vs Power
figure(2); subplot(2,3,[3 6]); hold on; title('Thrust vs Power'); xlabel('Power / W'); 
ylabel('Thrust / N'); 
ylim([0 max([max(S08L.T), max(S10L.T), max(S12L.T), max(APC.T), max(BASE.T), max(S12S.T)])]);
xlim([0 max([max(S08L.P), max(S10L.P), max(S12L.P), max(APC.P), max(BASE.P), max(S12S.P)])]);
% plot(S08.rpmmean, S08.FOM, '-or');
% plot(S10.rpmmean, S10.FOM, '-ok');

plot(S08L.P, S08L.T, '-xr');
plot(S10L.P, S10L.T, '-xk');
plot(S12L.P, S12L.T, '-xb');
plot(S12S.P, S12S.T, 'xb');
plot(BASE.P, BASE.T, '-.g');
plot(APC.P, APC.T, '-g');
% plot(APC.Re, APC.FOM, 'xg');
legend('\sigma = 0.8', '\sigma = 1.0', '\sigma = 1.2', '\sigma = 1.2 (LONG)', 'Baseline Propellor', 'APC Propellor', 'Location', 'southeast');
% legend('SHORT Sigma = 0.8', 'SHORT Sigma = 1.0', 'SHORT Sigma = 1.2', 'LONG Sigma = 0.8', 'LONG Sigma = 1.0', 'LONG Sigma = 1.2', 'Baseline Propellor', 'APC Propellor', 'Location', 'east');

save('EXP_DATA.mat','S08S','S10S','S12S','S08L','S10L','S12L','S08','S10','S12','R10','R17','R20');
save('EXP_DATA_PROP.mat', 'BASE','APC');