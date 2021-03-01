%%%% BLUE BEAR DUCTED FAN DESIGN v1 %%%% 
% Dec 2020 for Blue Bear project
%   - Updated with altitude and new operating points

% Inputs
% - 1.1: Define phi and sigma for one design or set phi=0 to analyse array of designs
% - 1.2: Set hover thrust, cruise thrust and cruise speed
% - 1.3: Set casing radius and hub-to-tip ratio
% - 1.4: Set atmospheric conditions
% - 1.5: Set aerodynamic efficiency - shaft to fluid, and electrical efficiency - battery to shaft
% - 1.6: Set vortex design, blade counts, aspect ratio, Lieblein diffusion factor and how deviation is accounted for. 

% Outputs
% - Most useful information are in structures NRF and CRF (conventional and contra-rotating)
% - These contain 'd' which are related to hover and cruise design operating points, 'a'
%   which are aerodynamic outputs like flow velocity and angle, and 'g' which
%   has geometry including blade metal angles (chi), chord (c) and blade number (N)
%   chi are arrays with blade metal angles and have dimensions:
%   - radius
%   - 1: LE, 2: TE
%   - 1: Blade 1 2: Blade 2
%   c is chord in m with dimensions: radius and which blade

% How it works:
% - Control volume function is called to analyse hover and cruise. Design maps are plotted if an array 
%   of designs is present 
% - Velocity triangle function called to find flow velocities and angles for conventional fan (NRF) and 
%   counter rotating fan (CRF)
% - "DuctedFanDesign_sdg_140520.ppt" has more details on equations used and comparison of two designs
   
% Still to do (11/5):
% - Velocity triangle function needs accompanying slides,  more comments, 
% - Varying alpha1 at inlet
% - Make design maps more flexible, i.e. explore casing diameter, thrust

% Clear workspace and load eVTOL envrionment
clear; clc; 
close all;
eVTOL;

% Run 27 design cube
nn = 3;
phi = repmat(reshape([0.55 0.7 0.85],[nn 1 1]),[1 nn nn]);
sig = repmat(reshape([0.8 1 1.2],[1 nn 1]),[nn 1 nn]);
r_cas = repmat(reshape([0.5 0.7 0.9],[1 1 nn]),[nn nn 1]);
for o = 14%1:numel(phi)
    
% Make new figure
h = figure();

%% -- 1. INPUTS -- %%

% ::: 1.1 Design Name ::: %
d.name = ['P' num2str(phi(o)*100) 'S' num2str(sig(o)*100) 'R' num2str(r_cas(o)*100)];

% ::: 1.1 Design ::: %
d.phi = phi(o); % 0.285; % Flow coefficient, set to zero for array of designs 
%                       velocity triangles do not work for this option (yet)
d.sig = sig(o); %0.5; % Area ratio of exit nozzle, sig = A4/A1. Eq (1) in .ppt

if d.phi==0
   phi1d = (0.2:0.1:1)';
   sig1d = (0.6:0.1:1.4)';
   d.phi = repmat(phi1d',[size(sig1d),1]);
   d.sig = repmat(sig1d,[1,size(phi1d)]);
end

% Split d into take-off and cruise
cruise=d; dtkoff=d;

% ::: 1.2 Performance requirements ::: %
cruise.T = 2000; % Thrust required in N
cruise.u0 = 51.4; % 50 kts
cruise.alt = 1524;
[cruise.atm.T,cruise.atm.p,cruise.atm.ro,~] = Altitude(cruise.alt); % At sea level

dtkoff.T = 5500; %327; % Thrust for forward flight in N
dtkoff.u0 = 24.7; % Forward flight speed in m/s
dtkoff.alt = 0;
[dtkoff.atm.T,dtkoff.atm.p,dtkoff.atm.ro,~] = Altitude(dtkoff.alt); % At 5000ft

% ::: 1.3 Geometry ::: %
g.rcas = r_cas(o); % casing radius in m
g.htrat = 0.4; % hub to tip ratio

% ::: 1.4 Atmosphere ::: %
% Remove and replaced with Altitude function K

% ::: 1.5 Efficiency and loading ::: %
eta.aero = 0.88; % efficiency of aerodynamic components, i.e. shaft to fluid
eta.elec = 0.85; % efficiency of electrical components, i.e. battery to shaft

% ::: 1.6 3D Blade Design ::: %
dev = 1; % d is whether to design with constant deviation across span 
%        - d = 0: chord is set to given constant diffusion factor across span
%        - d = 1:  chord is set to give constant deviation across span based on midspan value

CRF.AR = 3;
CRF.DF = 0.4;
CRF.n = 0.3; % exponent for vortex design of contra-rotating fan, n=1 is free, n=-1 is forced
CRF.blades = [0 0]; % Manually set number of blades, as above for NRF
% CRF.cman = [1.4,1;1.2,1]; % Scale calculated chord distribution [hub B1,cas B1;hub B2, casB2]. 
  
%% ---- 2. CONTROL VOLUME ANALYSIS ---- %%
% ::: 2.1 Take Off ::: %
[cruise,g,eta_TO] = CVanalysis(cruise,g,eta);
CRF.cruise = cruise;

% ::: 2.2 Cruise ::: %
[dtkoff,g,eta_cru] = CVanalysis(dtkoff,g,eta);
CRF.dtkoff = dtkoff;

%% ---- 3. Velocity Triangles ---- %%
if size(d.phi) == 1   
   % ::: 3.2 Counter-rotating compressor with rotor-rotor ::: %
   k=0.5; % k = 0.5 for counter-rotating 50-50 split between Rotor A and B 
   [CRF.a,CRF.g,CRF.d] = veltriangles(cruise,g,k,CRF,dev,1);
   
else
   disp('Velocity triangles cannot be calculated for multiple designs (yet)');
end

%% ---- 4. Use Mafiomoli Sections to Plot Blades ---- %%
if size(d.phi) == 1 
   CRF.bl1 = Make_Blade(CRF.g.chi(:,1,1),CRF.g.chi(:,2,1),CRF.g.r,CRF.g.c(:,1),1,h);    
   CRF.bl2 = Make_Blade(CRF.g.chi(:,1,2),CRF.g.chi(:,2,2),CRF.g.r,CRF.g.c(:,2),2,h);    
else
   disp('Blade geometry cannot be calculated for multiple designs (yet)');
end

%% ---- 5. Save blade designs to file ---- %%

% % Save design, rename structures to "d"
d = CRF; save([dr.geom 'CRF_' CRF.d.name '.mat'],'d')

% Add title to plots
subtitle(CRF.d.name);
gtg
% Close even numbered figures
% if mod(o,2) == 0; close(h); end;

% Clear old variables
clear d; clear CRF; 

end

% Save designs
save([dr.geom 'design_cube.mat'],'phi','sig','r_cas');


