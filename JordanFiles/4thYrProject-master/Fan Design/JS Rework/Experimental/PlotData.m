%% v1.0
% Graphs
% - T_ND vs sig
% - FOM vs sig

load('V1');

% Vortex Parameter
figure(1); 
subplot(2,1,1); hold on; grid on;
title('Vortex Parameter effect on Non-Dimensional Thrust vs Diffusion Factor')
xlabel('\sigma'); ylabel('T_{ND}');
plot(V1.P10.sig, V1.P10.T_ND, 'bo');
plot(V1.P17.sig, V1.P17.T_ND, 'ro');
plot(V1.P20.sig, V1.P20.T_ND, 'ko');
legend('Vortex Parameter = 1.0', 'Vortex Parameter = 1.75', 'Vortex Parameter = 2.0');

subplot(2,1,2); hold on; grid on;
title('Vortex Parameter effect on Figure of Merit vs Diffusion Factor')
xlabel('\sigma'); ylabel('M_f');
plot(V1.P10.sig, V1.P10.FOM, 'bo');
plot(V1.P17.sig, V1.P17.FOM, 'ro');
plot(V1.P20.sig, V1.P20.FOM, 'ko');
legend('Vortex Parameter = 1.0', 'Vortex Parameter = 1.75', 'Vortex Parameter = 2.0');

% Diffusion Factor
figure(2); 
subplot(2,1,1); hold on; grid on;
title('Diffusion Factor effect on Non-Dimensional Thrust vs RPM')
xlabel('\sigma'); ylabel('T_{ND}');
plot(V1.S08.sig, V1.S08.T_ND, 'bo');
plot(V1.S10.sig, V1.S10.T_ND, 'ro');
plot(V1.S12.sig, V1.S12.T_ND, 'ko');
legend('Diffusion Factor = 0.8000', 'Diffusion Factor = 1.0000', 'Diffusion Factor = 1.2649');

subplot(2,1,2); hold on; grid on;
title('Diffusion Factor effect on Figure of Merit vs Diffusion Factor')
xlabel('\sigma'); ylabel('M_f');
plot(V1.S08.sig, V1.S08.FOM, 'bo');
plot(V1.S10.sig, V1.S10.FOM, 'ro');
plot(V1.S12.sig, V1.S12.FOM, 'ko');
legend('Diffusion Factor = 0.8000', 'Diffusion Factor = 1.0000', 'Diffusion Factor = 1.2649');

%% v2.0
% Graphs
% - T_ND vs Vortex Parameter
% - FOM vs Vortex Parameter

load('V2');