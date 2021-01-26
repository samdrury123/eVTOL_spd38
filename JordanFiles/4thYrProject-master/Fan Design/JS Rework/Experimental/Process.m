%% Process Data
% - Form non-dimensionals
% - Process

rho = 1.225;

%% v1.0
% Variables
% - T_ND
% - FOM
% - Power Efficiency (mf/mf_ideal)
% - phi

% Load data
load('V1');

% Calculate flow area
V1.rc = 57.7e-3; 
V1.rh = 20e-3;
V1.area = calc_area(V1.rc, V1.rh);

V1.P10.rc = V1.rc;
V1.P17.rc = V1.rc;
V1.P20.rc = V1.rc;

V1.S08.rc = V1.rc;
V1.S10.rc = V1.rc;
V1.S12.rc = V1.rc;

% Phi mean-line
% V1.P10 = calc_phi(V1.P10);
% V1.P17 = calc_phi(V1.P17);
% V1.P20 = calc_phi(V1.P20);
% 
% V1.S08 = calc_phi(V1.S08);
% V1.S10 = calc_phi(V1.S10);
% V1.S12 = calc_phi(V1.S12);

% Non-dim Thrust
% V1.P10 = calc_T_ND(V1.P10);
% V1.P17 = calc_T_ND(V1.P17);
% V1.P20 = calc_T_ND(V1.P20);
% 
% V1.S08 = calc_T_ND(V1.S08);
% V1.S10 = calc_T_ND(V1.S10);
% V1.S12 = calc_T_ND(V1.S12);

% Figure of Merit
V1.P10.FOM = calc_FOM(V1.P10.T, V1.P10.P, V1.area);
V1.P17.FOM = calc_FOM(V1.P17.T, V1.P17.P, V1.area);
V1.P20.FOM = calc_FOM(V1.P20.T, V1.P20.P, V1.area);

V1.S08.FOM = calc_FOM(V1.S08.T, V1.S08.P, V1.area);
V1.S10.FOM = calc_FOM(V1.S10.T, V1.S10.P, V1.area);
V1.S12.FOM = calc_FOM(V1.S12.T, V1.S12.P, V1.area);

% Power efficiency
V1.P10 = calc_powereta(V1.P10);
V1.P17 = calc_powereta(V1.P17);
V1.P20 = calc_powereta(V1.P20);

V1.S08 = calc_powereta(V1.S08);
V1.S10 = calc_powereta(V1.S10);
V1.S12 = calc_powereta(V1.S12);

save('V1','V1');
%% v2.0
% Variables
% - T_ND
% - FOM
% - Power Efficiency
% - 

% Load data
load('V2');

V2.rc = 60e-3;
V2.rh = 20e-3;
V2.area = calc_area(V2.rc, V2.rh);

V2.P12.rc = V2.rc;
V2.P20.rc = V2.rc;

% Phi mean-line
V2.P12 = calc_phi(V2.P12);
V2.P20 = calc_phi(V2.P20);

% Non-dim Thrust
V2.P12 = calc_T_ND(V2.P12, V2.area);
V2.P20 = calc_T_ND(V2.P20, V2.area);

% Figure of Merit
V2.P12.FOM = calc_FOM(V2.P12.T, V2.P12.P, V2.area);
V2.P20.FOM = calc_FOM(V2.P20.T, V2.P20.P, V2.area);

% Power Efficiency
V2.P12 = calc_powereta(V2.P12);
V2.P20 = calc_powereta(V2.P20);

% Static Pressure Rise
V2.P12 = calc_pressurerise(V2.P12);
V2.P20 = calc_pressurerise(V2.P20);

% Inlet Vx
V2.P12 = calc_inletVx(V2.P12);
V2.P20 = calc_inletVx(V2.P20);



save('V2','V2');
%% v3.0
% Variables
% - T_ND - NO
% - FOM
% - Power Efficiency 

load('V3');

V3.P20.rc = V3.rc;

% Non-dim Thrust
% V3.P20.T_ND = calc_T_ND(V3.P20.T, V3.P20.rpmmean, V3.area, V3.rc, V3.rh);

% Figure of Merit
V3.P20.FOM = calc_FOM(V3.P20.T, V3.P20.P, V3.area);

% Power Efficiency
V3.P20 = calc_powereta(V3.P20);

% Static Pressure Rise
V3.P20 = calc_pressurerise(V3.P20);

% Inlet Vx
V3.P20 = calc_inletVx(V3.P20);

save('V3','V3');
%% Prop
prop_radius = 10 * 0.0254 / 2;
% Baseline
base_c = 0.0234;
% APC
apc_c = 0.026;
area = pi*(prop_radius^2);
data_apc = load('APC.mat');
data_base = load('Base.mat');