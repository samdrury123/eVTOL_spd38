
%% Global Variables
global flagplot

flag            = 1; % 2 - fixed NR
flagplot        = 0;% 0 - no plots 1 - plots

if flag == 2
    NR  = 1.41;
end

% Mean-line Design
r_h             = 11.6e-3; % m Hub radius
r_c             = 35e-3; % m Casing radius

eta.flow        = 0.5;


psi_mid         = 0.35;
phi_mid         = 0.6;  
DF_r            = linspace(0.45,0.40,N)';
DF_s            = linspace(0.40,0.40,N)';

Vinf            = 20;

% Metal Angles
tip_gap_percent = 3;

r.I1            = linspace(-5,-5,N)';  
s.I2            = linspace(-5,-5,N)' ;


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

Output.phi      = phi;
Output.psi      = psi;
Output.R        = R;
Output.DF_r     = DF_r;
Output.DF_s     = DF_s;
Output.rBlades  = r.blades;
Output.sBlades  = s.blades;
Output.thrust   = Thrust;
Output.power    = Power;
Output.camber   = camber(3);
