%% DATA COLLECTION INSTRUCTION
% COLLECTION OF DATA
% - Run Widget > Analyse:
%       select  SYS_STATUS.battery_current
%               SYS_STATUS.battery_voltage
% - Check MANUAL KILL SWITH ENGAGED (G)
% - Begin fan by typing from the Mavlink console (QGROUNDCONTROL)
%       "pwm test -p [1000 - 2000] -c 3"
% - Run 'EXP_COLLECT.m' and follow instructions. Be sure to close each
%   message box AFTER completing the instruction.
% - Start fan by disengaging MANUAL KILL SWITCH
% - Click 'START LOGGING' on QGC and save as 'power' in ./Logs
% 
% STOP RECORDING
% - 'STOP LOGGING' on QGC
% - RPi will automatically stop recording data and save file. EXP_COLLECT
% will scp the file from the RPi to the local machine

clear all; %#ok<CLALL>
close all;

%% FILE NAME
file = '201103_test';

waitfor(msgbox('Run STATIC_TEST.py on RPI','Info','warn'));
disp('Initialising DSA, PLEASE WAIT');
%getpressure;
disp('Processing Data...')
disp('Collecting RaspberryPi data (Thrust/RPM)...');
rpi = scp_simple_get('raspberrypi.local','pi','raspberry','STATIC_TEST.txt','.\','~/drone/');

%% TEST DEFINITION
% 08, 10, 12
% 08S, 10S, 12S
% 08L, 10L, 12L
SIGMANAME = 'N10';
if strcmp(SIGMANAME(1), 'R') || strcmp(SIGMANAME(1), 'N')
    sigma = 1.0; % [0.8, 1.0 1.2649]
    intake = 'LONG';
    exit = 'SHORT';
end
pwm = 1400;

LCfit = load('LC_FIT.mat','LC_cal');
thrustfit = LCfit.LC_cal;

%% PROCESS DATA
[timepwr, current, voltage] = importpower();
[timerpm, RPM, thr] = importrpm('TEST');
[pks, ind] = findpeaks(RPM);

thr = 9.81*thrustfit(thr)/1000;

time2 = linspace(0,max(timerpm),length(timerpm));

ind = ind(pks > 0.8);
n = length(ind);
timerpm = time2(ind);

current = current/100;
voltage = voltage/1000;
systemresist = voltage/current;

%% CALCULATE RPM
% load('RPMTEST.mat'); % RPM RAW TEST DATA

if var(RPM) == 0
    rpmsensor = 0;
else
    rpmsensor = 1;
end

switch rpmsensor
    case 0
    case 1
        averagewindow = 1;

        lim = averagewindow;
        timerpm = timerpm(lim:end-lim);
        
        separation = zeros([1, n-1]);
        
        for j = 1:n-1
            separation(j) = (time2(ind(j+1)) - time2(ind(j)));
        end

        timerpm(separation<0.002) = [];
        separation(separation<0.002) = [];

        f = length(separation);

        averagewindow = 2;

        for i = 1+averagewindow:f-averagewindow+1
            separation(i) = (sum(separation(i-averagewindow:i+averagewindow-1))/(2*averagewindow));
        end

        frequency = 1./separation;
        rpmmeasure = zeros([1, f]);

        for i = 1+averagewindow:f-averagewindow+1
            rpmmeasure(i) =  0.5 * 60 * (sum(frequency(i-averagewindow:i+averagewindow-1))/(2*averagewindow));
        end

        averagewindow = 2;

        k = length(rpmmeasure);

        liml = averagewindow;
        limh = averagewindow-1;

        for l = 1+averagewindow:k-averagewindow+1
            rpmmeasure(l) =  (sum(rpmmeasure(l-liml:l+limh-1)) / (2*averagewindow - 1));
        end
        
        rpmmean = mean(rpmmeasure(rpmmeasure > 0.8*max(rpmmeasure)));
end
%% CALCULATE OP and PLOT RAW DATA

power = voltage.*current;
P = mean(power(power > 0.8*max(power)));
T = mean(thr(thr > 0.8*max(thr)));
A = 0.0092;
rho = 1.225;

FOM = (T/P)*sqrt(T/(2*A*rho));
mass = 1000*T/9.81;

%% Get pressure data
load('EXP_PRESSURES.mat');
radius = [60 50 40 30] / 1000;

DSA.Pa = 1.0125e5;
DSA.time.P1 = e.P(:,9); DSA.time.P1b = e.P(:,14); DSA.time.P1c = e.P(:,15); DSA.time.P1d = e.P(:,16);
DSA.time.span.P1 = [DSA.time.P1 DSA.time.P1b DSA.time.P1c DSA.time.P1d];
DSA.time.P12 = e.P(:,10);
DSA.time.P2 = e.P(:,11);
DSA.time.PDI = e.P(:,12);
DSA.time.PDE = e.P(:,13);

DSA.P1 = mean(DSA.time.P1); DSA.P1b = mean(DSA.time.P1b); DSA.P1c = mean(DSA.time.P1c); DSA.P1d = mean(DSA.time.P1d);
DSA.span.P1 = [DSA.P1 DSA.P1b DSA.P1c DSA.P1d];
DSA.P12 = mean(DSA.time.P12);
DSA.P2 = mean(DSA.time.P2);
DSA.PDI = mean(DSA.time.PDI);
DSA.PDE = mean(DSA.time.PDE);

DSA.time.Vx = sqrt(2*(-DSA.time.P1)/rho);
DSA.Vx = mean(DSA.time.Vx);

switch rpmsensor
    case 1
        DSA.span.U = radius .* rpmmean*pi/30;
        DSA.U = (60e-3) * rpmmean*pi/30;

        % Flow Coefficients
        DSA.time.psi_ts = (DSA.time.P2 ./ (0.5*rho*DSA.U^2));
        DSA.time.phi = DSA.time.Vx/DSA.U;
        DSA.time.psi = 0.5*(DSA.time.psi_ts + DSA.time.phi.^2);
        DSA.time.CP1 = DSA.time.P1/(0.5*rho*DSA.U^2); DSA.time.CP12 = DSA.time.P12/(0.5*rho*DSA.U^2); DSA.time.CP2 = DSA.time.P2/(0.5*rho*DSA.U^2);
        DSA.time.CdelP = (DSA.time.P1 - DSA.time.P2)/(0.5*rho*DSA.U^2);
        DSA.time.DdelV = sqrt(2*(DSA.time.PDE - DSA.time.PDI)/rho);

        DSA.Vx = mean(DSA.time.Vx);
        DSA.phi = mean(DSA.time.phi);
        DSA.psi_ts = mean(DSA.time.psi_ts);
        DSA.psi = mean(DSA.time.psi);
        DSA.CP1 = DSA.P1/(0.5*rho*DSA.U^2); DSA.CP12 = DSA.P12/(0.5*rho*DSA.U^2); DSA.CP2 = DSA.P2/(0.5*rho*DSA.U^2);
        DSA.CdelP = (DSA.P1 - DSA.P2)/(0.5*rho*DSA.U^2);
        DSA.DdelV = sqrt(2*(DSA.PDE - DSA.PDI)/rho);

        % Span
        DSA.span.Vx = sqrt(2*(-DSA.span.P1)/rho);
        DSA.span.phi = DSA.span.Vx./DSA.span.U;
    case 0
end

%% Plot 
span = (radius - 20e-3) ./ (60e-3 - 20e-3);

figure(1); 
switch rpmsensor
    case 1
        disp(' '); disp(file);
        X = ['FOM = ' num2str(FOM) ';  RPM = ' num2str(rpmmean) ';  Thrust (g) = ' num2str(mass) ';'];
        disp(X); clear X;
        
        subplot(2,3,3); hold on; title('Spanwise Flow Coefficients'); ylabel('% Span'); ylim([0 1]);
        plot(DSA.phi, span, 'r'); plot(DSA.psi, span, 'b'); plot(DSA.psi_ts, span, 'k');
        legend('\phi','\psi','\psi_{TS}');
        subplot(2,3,6); hold on; title('Inlet Velocity Profile'); ylabel('% Span'); ylim([0 1]);
        plot(DSA.span.Vx, span);

        subplot(2,3,1); hold on; plot(timerpm, rpmmeasure); plot([0 max(timerpm)], [rpmmean rpmmean], 'k'); title('RPM'); xlabel('Time'); ylabel('RPM');
        subplot(2,3,4); hold on; plot(time2, 1000*thr/9.81);plot([0 max(time2)], [1000*T/9.81 1000*T/9.81], 'k'); ylim([0 1.2*1000*T/9.81]);title('Thrust (in grams)'); xlabel('Time'); ylabel('grams');
        subplot(2,3,2); hold on; plot(timepwr/1000, current); plot(timepwr/1000, voltage); ylim([0 1.2*max(voltage)]); title('Electrical Input'); xlabel('Time'); ylabel('Amps or Volts'); legend('Current', 'Voltage');
        subplot(2,3,5); hold on; plot(timepwr/1000, power); plot([0 max(timepwr)/1000], [P P], 'k'); ylim([0 1.2*P]); title('Power'); xlabel('Time'); ylabel('Watts');
    case 0
        disp(' '); disp(file);
        X = ['FOM = ' num2str(FOM) ';  Thrust (g) = ' num2str(mass) ';'];
        disp(X); clear X;
        
        subplot(2,2,[1 3]); hold on; plot(time2, 1000*thr/9.81);plot([0 max(time2)], [1000*T/9.81 1000*T/9.81], 'k'); ylim([0 1.2*1000*T/9.81]);title('Thrust (in grams)'); xlabel('Time'); ylabel('grams');
        subplot(2,2,2); hold on; plot(timepwr/1000, current); plot(timepwr/1000, voltage); ylim([0 1.2*max(voltage)]); title('Electrical Input'); xlabel('Time'); ylabel('Amps or Volts'); legend('Current', 'Voltage');
        subplot(2,2,4); hold on; plot(timepwr/1000, power); plot([0 max(timepwr)/1000], [P P], 'k'); ylim([0 1.2*P]); title('Power'); xlabel('Time'); ylabel('Watts');
end

savefig(file);
%% CALCULATE META, SAVE AND EXIT
if mean(voltage) > 12.5
    nomV = 14.8;
else
    nomV = 12;
end

save('EXP_META', 'pwm','SIGMANAME','nomV','intake','exit','sigma');

switch rpmsensor
    case 1
        save('EXP_OP', 'FOM', 'P', 'T', 'mass', 'rpmmean', 'DSA');
        save(file, 'FOM', 'P', 'T', 'mass', 'rpmmean', 'DSA');
        exppath = ['C:\Users\Sam Drury\Documents\CAMBRIDGE\IIB\eVTOL_Project\JordanFiles\4thYrProject-master\Drone Control\SAM\' file];
        save(exppath, 'FOM', 'P', 'T', 'mass', 'rpmmean', 'DSA');
    case 0
        save('EXP_OP', 'FOM', 'P', 'T', 'mass', 'DSA');
        save(file, 'FOM', 'P', 'T', 'mass', 'DSA');
        exppath = ['C:\Users\Sam Drury\Documents\CAMBRIDGE\IIB\eVTOL_Project\JordanFiles\4thYrProject-master\Drone Control\SAM\' file];
        save(exppath, 'FOM', 'P', 'T', 'mass', 'DSA');
end

clear all; %#ok<CLALL>

[P12, P20, P20_NS, P20_NS_S, CALC_P12, CALC_P20] = loaddata;
% EXP_UPDATE;