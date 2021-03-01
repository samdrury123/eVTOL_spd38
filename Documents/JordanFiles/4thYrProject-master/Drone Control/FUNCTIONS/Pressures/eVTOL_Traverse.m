%   eVTOL_Traverse                                
%                               
%   A script to run a traverses on the eVTOL rig
%
%   Written by Megan McEveley & James Taylor               
%   August 2019

% Clear workspace
close all; clear; clc; instrreset; delete(instrfindall);

%Load eVTOL data
eVTOL; dc = 0; settle_fac = 2; V_des = 30;

% Select probe
probe_name = 'FHP_M'; calib_name = 'FHP_M_30_fine_reduced.mat';

% Constants for current from 
park.port = 'COM1'; park.rate = 9600; park.term = 'CR'; park.nchan = 4;
park.current = [0.6 1.5 0.6 1.5]; park.resistance = [4.2 2 4.2 2]; park.inductance = [2.5 1 2.5 1];
park.V = [10 2 12 12]; park.A = [100 20 200 200];
park.d = [160000 2000 800000 800000];
park.c = [0 0 0 0];
park.limit = [0 240000  ;0 240000; -100000 240000 ; -200000 240000];

% Midspan radius
r_mid = 0.0875;

% Target RPM
rpm_target = 3000
N.V_shaft = 1; 

%% Five hole probe traverse settings 

% Record probe channel numbers
varnames = fieldnames(probe.(probe_name).N); 
for v = 1:length(varnames); N.(varnames{v}) = probe.(probe_name).N.(varnames{v}); end

% % Fine double passage at gauze exit
% st = 1; bl = [1 2]; dir_type = 'radial'; Iota = 25; ni = 1;
% nk = 19; t_min = -7.5; t_max = 0; t_clust = 'linear'; t_pitch = 7.5;
% 
% % Radial points distribution for half span traverse
% r_hub = 0.31; r_cas = 0.36; r_mid = 0.5 * (r_hub + r_cas);
% nj = 15; r_min = r_hub + 0.4e-3 + probe.(probe_name).D/2; r_max = r_mid; r_fac = [0.15 2.7];


%% Select operating point configuration

% Set tolerances in operating point to 0.25%
err.nrt = 0.0025; err.mrtp = 0.0025; 

% Channels for rig operating point
N_rig.P = N.P_Rint; N_rig.To = N.To_Rint; N_rig.Po = N.Po_Rint; N_rig.Pa = N.Pa; N_rig.V = N.V_shaft;
N_rig.P_FHP = N.P_FHP;
N_rig.PTH = 'Vaisala';

% Traverse plotting channels
N_plot.Po = N.Po_Sout; N_plot.P = N.P_Sout; N_plot.T = N.To_Rin;

% Set desired operating points
des.mrtp = 0.2685; des.nrt = 0.559;                              

% Frequency to check operating point at
n_check = 5;

% Set geometry name
geom_name = 'datum';


%% Run the traverse program

% Dry run to set up grid or traverse gear
% run_type = 'grid';
% run_type = 'off';

% Turn on profiler
profile clear; profile on;

% Call the script to initialise, execute, monitor and save the traverse
  Exp_Traverse

% View profiler
  profile viewer;



