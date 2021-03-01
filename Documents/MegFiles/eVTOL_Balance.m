% Balance eVTOL rig by logging once per rev and accelerometer

clear
close all
clc

% Load variables
eVTOL

% Channel AI0 is once per rev
% Channel AI1 is accelerometer

% A is 0 deg
% B is 45 deg
% C is 90 deg
% D is 135 deg
% E is 180 deg
% F is 225 deg
% G is 270 deg
% H is 315 deg

% A, C, E and G are used regulary for balancing blade sets
% B, D, F and H are used to balance the disc and are not often moved

% Log imbalance at 1000 rpm with filter set to 50. 
% 0.1 peak-peak imbalance is okay. Less than 0.05 is very good.

% Select channels
N.V_accel = 2;
N.V_shaft = 1;

% Select balancing positions
% th_bal = [0 90 180 270];
% th_bal = [45 135 225 315]; 
th_bal = 0:45:315;
% th_bal = [0 80 180 260]; % Yellow compressor positions - note labelled wrong way round

% NI-DAQ rates
pxie.rate = 40000; pxie.time = 2;

% Trial mass (g) and speed
th_trial = 0; mass_trial = 1.7;
rpm = 1650;

% Calculated factors for the Gibbons at 1500 rpm
% r_calib = 4.3598; ph0._trial = -151.7673;
% r_calib = 24.3403; ph_trial = 29.9896;

% Run balancing script
Exp_Balance

