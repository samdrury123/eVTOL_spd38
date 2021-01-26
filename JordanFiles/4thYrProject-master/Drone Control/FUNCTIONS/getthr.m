%% LOAD IN RPM DATA FIRST (Raw ADC output)
clear all;
close all;

importthr;

voltage = 16.56;
current = 200/voltage;
%% FIDDLE

%% PLOT
P = voltage*current;
T = mean(THR(25:end-25))/1000;
A = 0.0092;
rho = 1.225;

FOM = (T/P)*sqrt(T/(2*A*rho))

figure(2); hold on; plot([0 max(timethr)], [1000*T 1000*T]); plot(timethr, THR);title('Thrust'); xlabel('Time'); ylabel('mN');

