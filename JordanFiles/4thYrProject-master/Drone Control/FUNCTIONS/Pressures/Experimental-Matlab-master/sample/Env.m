% Gibbons Experimental Environment
% Written by James Taylor 2014

% Clear workspace and close instruments
% instrreset; close all; clear all; clc; delete(instrfindall);

%% Air properties

exp_air

%% Directories

% Set working directory 
directory = 'C:\Users\jvt24\Documents\Gibbons\Exp\';
cfd_directory = 'C:\Users\jvt24\Documents\Gibbons\TURBOSTREAM\';

%% Instrument information

% Use DSA-343 (Gibbons DSA)
dsa.ip = '129.169.103.86'; 

% Set default DSA logging variables
dsa.port = 23; dsa.nchan = 16;
dsa.period = 150; dsa.avg = 100;

% Vaisala PTH sensor
ptu.loc = 'LSL';

% Use temperature logger on USB
tc08.library = {'C:\Users\jvt24\Apps\UsbTC08sdk_r10_5_0_28\x64\usbtc08.dll' ...
    'C:\Users\jvt24\Apps\UsbTC08sdk_r10_5_0_28\usbtc08.h'};
tc08.type = 'K'; tc08.units = 2; tc08.period = 100; tc08.nchan = 8;

% Use NI PXIe 1073
pxie.sesh = 'ni';
pxie.dev = {'Slot2' 'Slot2' 'Slot2' 'Slot2' 'Slot2' 'Slot2' 'Slot2' 'Slot2'...
    'Slot3' 'Slot3' 'Slot3' 'Slot3' 'Slot3' 'Slot3' 'Slot3' 'Slot3'};
pxie.ai = [0:7 0:7]; pxie.nchan = 16; pxie.rate = 200e3;

% Accelerometer
acc.sens = 0.968; % V/g
acc.chan = 1; %AI 0

% Use Parker Paragon stepper motor driver
park.port = 'COM1'; park.rate = 9600; park.term = 'CR'; park.nchan = 3;
park.current = [2.4 0.6 0.6]; % Amps / phase
park.resistance = [0.72 6.2 4.2]; % Ohms / phase
park.inductance = [1 8.8 2.5]; % milli-henrys / phase

% Constants for each motor - steps per degree or metre
park.d = [160000 2000 200];
park.c = [-53589 -2000 -18000];

% Homing config for each motor
park.home.e = ['-' '+' '+'];
park.home.t = [1 1 0];
park.home.V = [-1 1 1];

% Limit configuration for each traverse position
park.limit.T3 = [-5556 2868 ; -1000 20000 ; -20000 20000];
park.limit.T4 = [-5576 2920 ; -8000 9000 ; -20000 20000];
park.limit.T5 = [-5624 2920 ; -2800 25500 ; -20000 20000];

%% Probe and rig measurements

% Five hole probes for paper with Sam
probe.FHP_HL.r_off = -10.07e-3; probe.FHP_HL.Iota_off = 0;
probe.FHP_HL.time = 0; probe.FHP_HL.D = 3e-3; probe.FHP_HL.arrange = 'plus';

probe.FHP_HM.r_off = -10.07e-3; probe.FHP_HM.Iota_off = 0;
probe.FHP_HM.time = 0.4; probe.FHP_HM.D = 2e-3; probe.FHP_HM.arrange = 'plus';

probe.FHP_HS.r_off = -10.07e-3; probe.FHP_HS.Iota_off = 0;
probe.FHP_HS.time = 3.5; probe.FHP_HS.D = 1e-3; probe.FHP_HS.arrange = 'plus';

% Offsets and time constants for probes
probe.FHP_G.r_off = -10.07e-3; probe.FHP_G.Iota_off = 0;
probe.FHP_G.time = 0.3; probe.FHP_G.D = 1.05e-3; probe.FHP_G.arrange = 'cross';

probe.ToP_A.r_off = -5.78e-3; probe.ToP_A.time = 2; 
probe.ToP_A.D = 2.3e-3;

probe.PoP_A.r_off = -5.57e-3; probe.PoP_A.time = 15; 
probe.PoP_A.L = 25e-3; probe.PoP_A.D = 0.303e-3;

probe.HW_A.r_off = -0.81e-3; probe.HW_A.Iota_off = 80; 
probe.HW_A.time = 1.29; probe.HW_A.D = 0;

probe.HW_B.r_off = -1.51e-3; probe.HW_B.Iota_off = 80; 
probe.HW_B.time = 1.29; probe.HW_B.D = 0;

% Rig measurements
r_hub = 0.310; r_cas = 0.3647; r_mid = 0.5 * (r_hub + r_cas);
x_2 = 0.1; x_3 = 0.1866; x_4 = 0.2466; x_5 = 0.3226; x_6 = 0.39;
x.T4 = x_4; x.T5 = x_5;
