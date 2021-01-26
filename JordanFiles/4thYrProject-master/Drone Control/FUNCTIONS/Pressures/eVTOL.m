%   eVTOL                                
%                               
%   Load eVTOL environment variables for CFD and experiment
%
%   Written by Megan McEveley              
%   July 2018


%% Set directories and job variables

% Set working directory on Windows computer
directory = 'C:\Users\mm2191\Desktop\testlog\190107_NHP_HWtests\';

% Root directory depending on host
% [~,host] = system('hostname'); host = strtrim(host);
% if strcmp(host,'rjm-225.eng.cam.ac.uk') == 1
%     dr.root = '/home4/mm2191/eVTOL/';
% end

% Sub directories saved in struct
% dr.ts = [dr.root 'TURBOSTREAM' filesep];
% dr.exp = [dr.root 'Exp' filesep];
% dr.geom = [dr.root 'Geometry' filesep];
% dr.pv = [dr.root 'PARAVIEW' filesep];
% dr.ag = [dr.root 'AUTOGRID' filesep];
% dr.fig = [dr.root 'Documents' filesep];


% Variables for Wilkes 1
% dr.hpc = '/scratch/jvt24/iGibbons/';
% job.host = 'login.hpc.cam.ac.uk';
% job.user = 'jvt24';
% job.directory = dr.ts;
% job.hpcdir = dr.hpc;
% job.version = '3.2.0';
% job.ppn = 2;
% job.mem = 63900;
% job.arch = 'tesla';

% Variables for Wilkes 2
% dr.hpc = '/rds/user/mm2191/hpc-work/eVTOL/';
% job.host = 'login-gpu.hpc.cam.ac.uk';
% job.user = 'mm2191';
% job.directory = dr.ts;
% job.hpcdir = dr.hpc;
% job.version = '3.6.2';
% job.ppn = 4;
% job.mem = 96000;
% job.arch = 'pascal';


%% Fluid properties

% Air properties for experiments
exp_air

% Gas constants
cp = 1003.308; ga = 1.40047; mu = 1.8000e-05;
cv = cp / ga; R = cp - cv;

% Global variables for post processing CFD
global Pref Tref
Pref = 101325.8;
Tref = 288.15;

%% Geometry
r_cas = 0.12;
r_hub = 0.055;
r_mid = (r_cas+r_hub)/2;
Apass = pi*(r_cas^2 - r_hub^2);

% Torque calibration
% cal_grd = 10;



%% Instrument information

% Use DSA-343 (eVTOL DSA)
dsa.ip = '129.169.103.101'; 

% Set default DSA logging variables
dsa.port = 23; dsa.nchan = 16;
dsa.period = 150; dsa.avg = 150;

% Vaisala PTH sensor
ptu.loc = 'LSL';

% Use temperature logger on USB
tc08.library = {'C:\Users\mm2191\Apps\UsbTC08sdk_r10_5_0_28\x64\usbtc08.dll' ...
    'C:\Users\mm2191\Apps\UsbTC08sdk_r10_5_0_28\usbtc08.h'};
tc08.type = 'K'; tc08.units = 2; tc08.period = 100; tc08.nchan = 8;

% Use NI PXIe 1073 for voltage readings
pxie.sesh = 'ni'; pxie.dev = repmat({'Dev1'},[1 16]);
pxie.type = 'V'; pxie.coupling = repmat({'DC'},[1 16]); 
pxie.ai = [0:7 0:7]; pxie.nchan = 16; pxie.rate = 20e3;


%% Instrument channel numbers

% Record total number of channels for all instruments
tc08.chan = 1:2; pxie.chan = 1:2; rtd.chan = 1:12;

% Pressure channels
N.P_FHP = 1:5;
N.Po_bell = 6;
N.P_bell = 7;
% N.Po_Sint = 7;
N.Po_Rin = 8;
N.Po_Rint = 9;
N.P_Rin = 10;
N.P_Rint = 11;
N.Po_Sou = 12;
N.Po_Sout = 13;
N.P_Sou = 14;
N.P_Sout = 15;
N.Pa = 16;

% Thermocouple channels
N.T_disc = 1 + length(rtd.chan);
N.T_casing = 2 + length(rtd.chan);

% RTD channels
% N.Ta = 1;
N.To_Rin = 2:5;
N.To_Rint = [1 6];
N.To_Sou = 8:11;
N.To_Sout = [7 12];

% Voltage channels
N.V_shaft = 1;
% N.V_Ha = 2;
% N.V_Pa = 3;


%% Probe values and measurements

% Five hole probes
% probe.FHP_H.r_off = 0; probe.FHP_H.Iota_off = -1.41;
% probe.FHP_H.time = 1.2; probe.FHP_H.D = 1.05e-3; probe.FHP_H.arrange = 'plus';
% N.P_cn = 1; N.P_ur = 2; N.P_dr = 3; N.P_dl = 4; N.P_ul = 5; probe.FHP_H.N = N;
% 
% probe.FHP_I.r_off = 5.913e-3; probe.FHP_I.Iota_off = -1.11;
% probe.FHP_I.time = 0.8; probe.FHP_I.D = 1.05e-3; probe.FHP_I.arrange = 'cross';
% N.P_cn = 1; N.P_ur = 2; N.P_dr = 3; N.P_dl = 4; N.P_ul = 5; probe.FHP_I.N = N;

probe.FHP_M.r_off = 0; probe.FHP_M.Iota_off = 0;
probe.FHP_M.time = 1.2; probe.FHP_M.D = 1.3e-3; probe.FHP_M.arrange = 'plus';
N.P_cn = 3; N.P_up = 4; N.P_dn = 5; N.P_lf = 7; N.P_rt = 6; probe.FHP_M.N = N;

% % probe.FHP_M.r_off = 0; probe.FHP_M.Iota_off = 0;
% % probe.FHP_M.time = 1.2; probe.FHP_M.D = 1.3e-3; probe.FHP_M.arrange = 'plus';
% N.P_cn = 1; N.P_up = 2; N.P_dn = 3; N.P_lf = 5; N.P_rt = 4; probe.FHP_M.N = N;

% probe.FHP_FERG.Iota_off = 0; probe.FHP_FERG.time = 0.5; probe.FHP_FERG.arrange = 'plus';
% N.P_cn = 1; N.P_up = 2; N.P_dn = 4; N.P_lf = 3; N.P_rt = 5; probe.FHP_FERG.N = N;

% Stagnation temperature probes
probe.To_A.r_off = -5.78e-3; probe.To_A.time = 2; 
probe.To_A.D = 2.3e-3;

% Stagnation pressure probes
probe.Po_A.r_off = 0; probe.Po_A.Iota_off = 0; probe.Po_A.time = 0; 
probe.Po_A.D = 1.05e-3; probe.Po_A.N.Po = 1;

% Static pressure probes
probe.P_A.r_off = 0; probe.P_A.Iota_off = 0; probe.P_A.time = 0; 
probe.P_A.D = 1.05e-3; probe.P_A.N.P = 1;

% Hot-wires
probe.HW_C.r_off = 15.545e-3; probe.HW_C.Iota_off = 75.81; 
% probe.HW_C.r_off = 6.97e-3; probe.HW_C.Iota_off = 75.81; 
probe.HW_C.time = 1.25; probe.HW_C.D = 0; probe.HW_C.N.V_hw = 2; probe.HW_C.N.V_hw_av = 3;
probe.HW_C.lpf = 1; probe.HW_C.hpf = 30e3; probe.HW_C.gain = 20;
probe.HW_C.ohr = 0.8; probe.HW_C.alpha = 0.0036;
pxie.V_thresh = 4;


%% Rig Measurments
% 
% % SMURF scale factor
% scale = 0.47861; 
% 
% % Meridional radii
% r_hub = 0.310; r_cas = 0.3647; r_mid = 0.5 * (r_hub + r_cas);
% 
% % Axial measurment locations
% x_1 = -0.025; x_2 = 0.114; x_3 = 0.184; x_4 = 0.25125; x_5 = 0.324; x_6 = 0.369;
% 
% % Blade locations
% x_IGV = 0.09251 + 0.5 * 7.98e-3 + 0.045;
% x_Rotor = 0.213657 + 0.5 * 8.025e-3;
% x_Stator = 0.273710 + 0.5 * 7.98e-3;
% x_mid_IGV = -0.026234480 - 0.045; x_mid_Rotor = -0.01828903; x_mid_Stator = -0.01625015;
% 
% % Nose bullet axial gap size
% x_gap = [0.188 0.189]; dr_gap = 11.5e-3;
% 
% % Meridional lines
% xr_hub = [linspace(-0.1,0.6,100)' ones(100,1) * r_hub];
% xr_cas = [linspace(-0.1,0.6,100)' ones(100,1) * r_cas];
% save([dr.geom 'xr_hub.mat'],'xr_hub') 
% save([dr.geom 'xr_cas.mat'],'xr_cas') 


%% Operating points for CFD

% % Nominal shaft speed
% rpm_des = 3500; U_mid = r_mid * 2 * pi * rpm_des / 60;
% 
% % Flow conditions relative to IGV inlet
% % m_p1 = 9.1084; p_p1 = 104200; % For R_inlet_Rig_260_R8_scaled.dat inlet condition
% m_p1 = 9.1084; p_p1 = 104600; % For GZ27_measured_inlet.dat inlet condition
% 
% 
% %% CFD mesh parameters
% 
% % Boundary conditions
% I.pout = p_p1; I.m = m_p1; I.rpm = [rpm_des 0];
% 
% % Iteration counts
% I.nstep = 100000; I.poisson_nstep = 10000; I.k = 100;
% 
% % Preallocate cell arrays
% B = cell(2,1); F = cell(2,1); D = cell(2,1);
% 
% % Define fillet parameters
% s_cl = [0 1]; rad = [1 1] * 2e-3;
% 
% % Rotor hub fillet
% F{1}.hub.s1 = s_cl; F{1}.hub.s2 = s_cl; F{1}.hub.rad1 = rad; F{1}.hub.rad2 = rad; 
% F{1}.flag = [1 0]; F{1}.tol = 1e-6;    
% 
% % Stator casing fillet
% F{2}.cas.s1 = s_cl; F{2}.cas.s2 = s_cl; F{2}.cas.rad1 = rad; F{2}.cas.rad2 = rad; 
% F{2}.flag = [0 1]; F{2}.tol = 1e-6;
% 
% % Define gap parameters
% D{1}.cas = 0.8e-3; D{1}.flag = [0 1];
% D{2}.hub = 0.6e-3; D{2}.flag = [1 0];
% 
% % Add interplatform gap between nose bullet and disc
% D{1}.shroud = [3 0]; D{1}.dx_mp = 2e-3;  D{1}.dr_shroud = dr_gap; D{1}.x_shroud = x_gap;
% 
% % Global mesh parameters for a fine iGibbons grid
% M.dx_row = 3.3e-3;                                                  % Axial constant
% M.n_si1 = [23 0]; M.n_si2 = [17 0];                                 % Shroud z-constant lines
% M.n_zi = [0 17];  M.n_zo = [17 37];                                 % z-constant lines
% M.n_j = 121; M.n_jc = 30;                                           % Flow path and constant
% M.n_g = 25; M.n_f = 25; M.n_fb = 13; M.n_gb = 13;                   % Gap, fillet and butterfly
% M.n_bu = 69; M.n_bl = 101; M.n_b = 25; M.r_b = 1.2;                 % Blade blocks
% M.n_ai = 21; M.n_aiu = 21; M.n_aid = 21;                            % Azimuthal inlet block
% M.n_ao = 21; M.n_aou = 21; M.n_aod = 21;                            % Azimuthal outlet block
% M.n_si = 21; M.n_so = 9;                                            % Streamwise blocks
% M.e_j = 3e-6; M.e_b = 3e-6; M.e_f = 2e-4; M.e_s = 5e-6;             % First point spacing
% M.c_sh = 10;                                                        % Axial re-clustering
% M.i_j = 2;                                                          % Spanwise interpolation
% M_Y1 = M;
% 
% % Refined mesh to check grid independance
% M.dx_row = 3.5e-3;                                                  % Axial constant
% M.n_si1 = [27 0]; M.n_si2 = [21 0];                                 % Shroud z-constant lines
% M.n_zi = [0 21];  M.n_zo = [21 41];                                 % z-constant lines
% M.n_j = 125; M.n_jc = 30;                                           % Flow path and constant
% M.n_g = 29; M.n_f = 29; M.n_fb = 17; M.n_gb = 17;                   % Gap, fillet and butterfly
% M.n_bu = 77; M.n_bl = 125; M.n_b = 29; M.r_b = 1.15;                % Blade blocks
% M.n_ai = 25; M.n_aiu = 25; M.n_aid = 25;                            % Azimuthal inlet block
% M.n_ao = 25; M.n_aou = 25; M.n_aod = 25;                            % Azimuthal outlet block
% M.n_si = 25; M.n_so = 13;                                           % Streamwise blocks
% M.e_j = 5e-6; M.e_b = 5e-6; M.e_f = 2e-4; M.e_s = 8e-6;             % First point spacing
% M.c_sh = 10;                                                        % Axial re-clustering
% M.i_j = 2;                                                          % Spanwise interpolation
% M_YF = M;
% 
% % Global mesh parameters for a coarse iGibbons grid
% M.dx_row = 3.3e-3;                                                  % Axial constant
% M.n_si1 = [23 0]; M.n_si2 = [17 0];                                 % Shroud z-constant lines
% M.n_zi = [0 17]; M.n_zo = [17 37];                                  % z-constant lines
% M.n_j = 101; M.n_jc = 30;                                           % Flow path and constant
% M.n_g = 25; M.n_f = 25; M.n_fb = 13; M.n_gb = 13;                   % Gap, fillet and butterfly
% M.n_bu = 65; M.n_bl = 97; M.n_b = 25; M.r_b = 1.15;                 % Blade blocks
% M.n_ai = 21; M.n_aiu = 17; M.n_aid = 17;                            % Azimuthal inlet block
% M.n_ao = 21; M.n_aou = 17; M.n_aod = 17;                            % Azimuthal outlet block
% M.n_si = 17; M.n_so = 9;                                            % Streamwise blocks
% M.e_j = 1e-5; M.e_b = 1e-5; M.e_f = 2e-4; M.e_s = 2e-5;             % First point spacing
% M.c_sh = 10;                                                        % Axial re-clustering
% M.i_j = 2;                                                          % Spanwise interpolation
% M.skew = 12;                                                        % Skewness limit
% M_Y7 = M;
% 
% % Mesh boundary parameters
% x_out = 0.45; x_in = 0.17;
% 
% % Fix rotor-stator mixing plane location
% x_mp = [0.253 0.252 0.247];
% 
% % Coordinates for running chics
% coords{1}.bids = 2; coords{1}.xr = [x_3 0 ; x_3 1];
% coords{2}.bids = 11; coords{2}.xr = [x_5 0 ; x_5 1];

