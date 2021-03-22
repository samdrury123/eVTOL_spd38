%   eVTOL
%                               
%   Load eVTOL environment variables for CFD and experiment
%
%   Written by James Taylor & Sam Grimshaw               
%   May 2020


%% Set directories and job variables

% Root directory depending on host
[~,host] = system('hostname'); host = strtrim(host);
if strcmp(host,'rjm-125.eng.cam.ac.uk') == 1
    dr.root = '/home/jvt24/eVTOL/';
    job.user = 'jvt24'; 
    dr.hpc = '/rds/project/amjd2/rds-amjd2-compressors/eVTOL/';
elseif strcmp(host,'DESKTOP-4OIBD51') == 1
    dr.root = 'C:\Users\jvt24\Documents\eVTOL\';
    job.user = 'jvt24';
    dr.hpc = '/rds/project/amjd2/rds-amjd2-compressors/eVTOL/';
elseif strcmp(host,'WL-JT-002') ==1
    dr.root = '/home/jc2071/eVTOL/';
    job.user = 'jc2071';
    dr.hpc = '/rds/hpc-work/eVTOL/';
elseif strcmp(host,'Jonathans-MacBook-myname-7.local') ==1
    dr.root = '/Users/MCJON/OneDrive - University of Cambridge/eVTOL/';
    job.user = 'jc2071';
    dr.hpc = '/rds/hpc-work/eVTOL/';
elseif strcmp(host,'LAPTOP-S3CMJ7V6') ==1
    dr.root = 'C:\Users\Sam Drury\Documents\CAMBRIDGE\IIB\eVTOL\';
    job.user = 'spd38';
    dr.hpc = '/rds/hpc-work/eVTOL/';
elseif strcmp(host,'LAPTOP-S3CMJ7V6') ==1
    dr.root = '/home/spd38/eVTOL/';
    job.user = 'spd38';
    dr.hpc = '/rds/hpc-work/eVTOL/';
end

% Sub directories saved in struct
dr.tb = [dr.root 'TBLOCK' filesep];
dr.ts = [dr.root 'TURBOSTREAM' filesep];
dr.exp = [dr.root 'Exp' filesep];
dr.geom = [dr.root 'Geometry' filesep];
dr.pv = [dr.root 'PARAVIEW' filesep];
dr.ag = [dr.root 'AUTOGRID' filesep];
dr.fig = [dr.root 'Documents' filesep];
dr.cad = [dr.root 'CAD' filesep];

% Variables for Wilkes 2
job.host = 'login-gpu.hpc.cam.ac.uk';
job.directory = dr.ts;
job.hpcdir = dr.hpc;
job.version = '3.6.3';
job.ppn = 4;
job.mem = 96000;
job.arch = 'pascal';

%% Fluid properties

% Air properties for experiments
% exp_air

% Gas constants
cp = 1003.308; ga = 1.40047; mu = 1.8000e-05;
cv = cp / ga; R = cp - cv;

% Global variables for post processing CFD
global Pref Tref
Pref = 101325.8;
Tref = 288.15;


%% CFD mesh parameters

% Boundary conditions
I.pout = Pref; I.inlet = [dr.ts 'inlet_2percent.dat']; I.guess = [1.1 50 0 0 2e5 0.0002];

% Iteration counts
I.nstep = 100000; I.poisson_nstep = 10000; 

% Preallocate cell arrays
B = cell(2,1); F = cell(2,1); D = cell(2,1);

% Define fillet parameters
s_cl = [0 1]; rad = [1 1] * 8e-3;

for n = 1:2
    F{n}.hub.s1 = s_cl; F{n}.hub.s2 = s_cl; F{n}.hub.rad1 = rad; F{n}.hub.rad2 = rad; 

    F{n}.cas.s1 = s_cl; F{n}.cas.s2 = s_cl; F{n}.cas.rad1 = rad; F{n}.cas.rad2 = rad; 

    F{n}.flag = [2 2]; F{n}.tol = 1e-6;
end


% Define gap parameters

D{1}.flag = [0 0]; D{2}.flag = [0 0];

% Global mesh parameters for ducted fan
M.dx_row = 25e-3;                                                   % Axial constant
M.n_zi = [21 9]; M.n_zo = [9 33];                                   % Z-constant lines
M.n_j = 89; M.n_jc = 30;                                            % Flow path and constant
M.n_g = 25; M.n_f = 21; M.n_fb = 13; M.n_gb = 13;                   % Gap, fillet and butterfly
M.n_bu = 65; M.n_bl = 65; M.n_b = 25; M.r_b = 1.25;                 % Blade blocks
M.n_ai = 21; M.n_aiu = 17; M.n_aid = 17;                            % Azimuthal inlet block
M.n_ao = 21; M.n_aou = 17; M.n_aod = 17;                            % Azimuthal outlet block
M.n_si = 17; M.n_so = 17;                                           % Streamwise blocks
M.e_j = 1.3e-5; M.e_b = 1e-5; M.e_f = 2e-3; M.e_s = 4e-5;           % First point spacing
M.i_j = 2;                                                          % Spanwise interpolation
M.skew = 12;                                                        % Skewness limit


