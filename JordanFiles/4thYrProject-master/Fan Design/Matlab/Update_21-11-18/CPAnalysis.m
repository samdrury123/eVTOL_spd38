clear all
close all

%% Global Variables
flag            = 2; % 2 - fixed NR
NR              = 1.411;
eta.flow        = 0.45;


% Mean-line Design
r_h             = 11.6e-3; % m Hub radius
r_c             = 35e-3; % m Casing radius

psi_mid         = 0.35;
phi_mid         = 0.6;  
DF_r            = linspace(0.45,0.40,5)';
DF_s            = linspace(0.40,0.40,5)';
p               = 1.1; % psi distribution exponent (0 - forced, 2 - free)

Vinf            = 0;

N               = 5;

% Metal Angles
tip_gap_percent = 3;

r.I1            = [-5 -5 -5 -5 -5]';  
s.I2            = [-5 -5 -5 -5 -5]' ;


%r.I1            = [0 0 0 0 0]';  
%delta1          = [0 0 0 0 0]';
%s.I2            = [0 0 0 0 0]' ;

% File Name
temp            = {'BladeGeom'};
pathname        = 'C:\Users\Fritz Kirchner\Google Drive\Project\MatLab Scripts\Blade Sizing\BladeDesign\Geometry';
file_name       = fullfile(pathname,temp);

%% Mean-line Design

Mean_Line_Design

delta1          = delta';

%% Metal Blade Angles 

MetalAngles

%% Print Blades

PrintGeometry

%% Display Results
CP_exp      = [-0.377148737, -0.43,	0.18397499,	0.220769992,	0.025572524];

close all





