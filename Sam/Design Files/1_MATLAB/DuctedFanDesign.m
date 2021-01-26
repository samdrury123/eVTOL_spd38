    %%%% DUCTED FAN DESIGN %%%% 
% Sept 2020 for Sam Drury %

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
% - Varying alpha1 at inlet?
% - Make design maps more flexible, i.e. explore casing diameter, thrust

% Clear workspace and load eVTOL environment
clear; clc; close;
% close all;
% eVTOL;


%% -- 1. INPUTS -- %%

% CRDF flag
CRDF_flag = false;

% ::: 1.1 Design Name ::: %
d.name = '1';

% ::: 1.1 Design ::: %
d.phi = 0.8; % Flow coefficient, set to zero for array of designs 
%                       velocity triangles do not work for this option (yet)
d.sig = 1.1; %0.5; % Area ratio of exit nozzle, sig = A4/A1. Eq (1) in .ppt. Not used if d.phi=0;
d.rotor_lean = 20; % Lean angle in degrees
d.rotor_sweep = 4; % Sweep angle in degrees
d.stator_lean = 20;
d.stator_sweep = 4;

if d.phi==0
   phi1d = (0.2:0.05:1.0)'; %%CHANGED TO 0.9 from 0.8
   sig1d = (0.4:0.05:1.2)';
   d.phi = repmat(phi1d',[size(sig1d),1]);
   d.sig = repmat(sig1d,[1,size(phi1d)]);
end

% Split d into hov and cru
dhov=d; dcru=d;

% ::: 1.2 Performance requirements ::: %
dhov.T = 10; % Thrust required in N (2.43kg/4, add a little)
dhov.u0 = 0;
dcru.T = 3; % Thrust for forward flight in N (arbitrary)
dcru.u0 = 3; % Forward flight speed in m/s

% ::: 1.3 Geometry ::: %
g.rcas = 0.06; % casing radius in m
g.htrat = 1/3; % hub to tip ratio

% ::: 1.4 Atmosphere ::: %
atm.ro = 1.225; % density of air in kg/m^3
atm.p = 101325; % atmospheric pressure in N/m^2
atm.T = 273+15; % atmospheric temperature in K

% ::: 1.5 Efficiency and loading ::: %
eta.aero = 0.81; % efficiency of aerodynamic components, i.e. shaft to fluid
eta.elec = 0.85; % efficiency of electrical components, i.e. battery to shaft

% ::: 1.6 3D Blade Design ::: %
dev = 1; % d is whether to design with constant deviation across span 
%        - d = 0: chord is set to given constant diffusion factor across span
%        - d = 1:  chord is set to give constant deviation across span based on midspan value
NRF.AR = 1.8; % Aspect ratio of blades, span to meanspan-chord
NRF.DF = 0.35; % Lieblein Diffusion Factor
NRF.n = 0.5; % exponent for vortex design of conventional fan, n=1 is free, n=-1 is forced
NRF.blades = [0 0]; % Manually set number of blades:
                    % - If set to 0, meanspan chord set by AR, pitch-to-chord by DF and blade count by meanspan pitch
                    % - If set to a number, pitch is calculated from blade count and chord from pitch-to-chord from DF 
% NRF.cman = %[0.8,1.5;2,1]; % Scale calculated chord distribution [hub B1,cas B1;hub B2, casB2]. 


CRF.AR = 2.7;
CRF.DF = 0.38;
CRF.n = 0.3; % exponent for vortex design of contra-rotating fan, n=1 is free, n=-1 is forced
CRF.blades = [13 11]; % Manually set number of blades, as above for NRF
CRF.cman = [1.4,1;1.2,1]; % Scale calculated chord distribution [hub B1,cas B1;hub B2, casB2]. 
  
%% ---- 2. CONTROL VOLUME ANALYSIS ---- %%
% ::: 2.1 Hover ::: %
[dhov,g,eta_hov] = CVanalysis(dhov,g,eta,atm);

% ::: 2.2 Cruise ::: %
[dcru,g,eta_cru] = CVanalysis(dcru,g,eta,atm);
NRF.dcru = dcru;

if CRDF_flag == true
    CRF.dcru = dcru;
end

% ::: 2.3 Plot design performance maps ::: %
if size(d.phi) > 1
    % Hover
    figure(1); set(gcf, 'color', 'w'); grid off; box on; set(gcf,'Position',[20 200 1000 400]);
    subplot(2,4,1); contourf(d.phi,d.sig,dhov.Pideal,20); colorbar
    xlabel('phi'), ylabel('sigma'); title('Hover Ideal Power (kW)')
    text(0.25,0.45,'Halve for each CRF row','color','w')
        
    subplot(2,4,2); contourf(d.phi,d.sig,dhov.rpm,20); colorbar
    xlabel('phi'), ylabel('sigma'); title('Hover Shaft speed (rpm)')
    text(0.25,0.45,'Halve for each CRF row','color','w')
    
    subplot(2,4,3); contourf(d.phi,d.sig,dhov.Tq,20); colorbar
    xlabel('phi'), ylabel('sigma'); title('Hover Torque (Nm)')
    
    subplot(2,4,4); contourf(d.phi,d.sig,dhov.psitt,20); colorbar
    xlabel('phi'), ylabel('sigma'); title('Hover Stage Loading Coefficient')

    % Forward flight
    subplot(2,4,5); contourf(d.phi,d.sig,dcru.Pideal,20); colorbar
    xlabel('phi'), ylabel('sigma'); title('Cruise Ideal Power (kW)')
    text(0.25,0.45,'Halve for each CRF row','color','w')
    
    subplot(2,4,6); contourf(d.phi,d.sig,dcru.rpm,20); colorbar
    xlabel('phi'), ylabel('sigma'); title('Cruise Shaft speed (rpm)')
    text(0.25,0.45,'Halve for each CRF row','color','w')
    
    subplot(2,4,7); contourf(d.phi,d.sig,dcru.phicru,20); colorbar
     xlabel('phi'), ylabel('sigma'); title('Cruise Flow Coefficient')   
    
%      subplot(2,3,6); contourf(d.phi,d.sig,eta_cru.prop,20); colorbar
%      xlabel('phi'), ylabel('sigma'); title('Propulsive Efficiency')   
     subplot(2,4,8); contourf(d.phi,d.sig,dcru.psitt,20); colorbar
     xlabel('phi'), ylabel('sigma'); title('Cruise Stage Loading Coefficient')
     
 % SAM
    figure(2); set(gcf, 'color', 'w'); grid off; box on; set(gcf,'Position',[20 200 400 400]);
    contourf(d.phi,d.sig,dhov.Pideal,20); colorbar; hold on;
    xlabel('\phi'), ylabel('\sigma'); title('Hover Ideal Power (kW)')
    phi_des = 0.8; sig_des = 1.1;
    plot(phi_des, sig_des, 'ro', 'MarkerSize', 5, 'LineWidth', 5);
    hold off;
    %text(0.25,0.45,'Halve for each CRF row','color','w')

end

%% ---- 3. Velocity Triangles ---- %%
if size(d.phi) == 1   
   % ::: 3.1 Conventional compressor with rotor-stator ::: %
   k = 1; % k is Blade A proportion of rpm, k = 1 for conventional    
   [NRF.a,NRF.g,NRF.d] = veltriangles(dhov,g,k,NRF,dev,1,1);
   
   % ::: 3.2 Counter-rotating compressor with rotor-rotor ::: %
   if CRDF_flag == true
       k=0.5; % k = 0.5 for counter-rotating 50-50 split between Rotor A and B 
       [CRF.a,CRF.g,CRF.d] = veltriangles(dhov,g,k,CRF,dev,1,2);
   end
   
else
   disp('Velocity triangles cannot be calculated for multiple designs (yet)');
end

%% ---- 4. Use Mafiomoli Sections to Plot Blades ---- %%
if size(d.phi) == 1
   NRF.bl1 = Make_Blade(NRF.g.chi(:,1,1),NRF.g.chi(:,2,1),NRF.g.r,NRF.g.c(:,1),1,1,d);    
   NRF.bl2 = Make_Blade(NRF.g.chi(:,1,2),NRF.g.chi(:,2,2),NRF.g.r,NRF.g.c(:,2),2,1,d);
   txt = {'Blades:', NRF.g.N'};
   %text(0.01,0,txt)
    
   if CRDF_flag == true
       CRF.bl1 = Make_Blade(CRF.g.chi(:,1,1),CRF.g.chi(:,2,1),CRF.g.r,CRF.g.c(:,1),1,2,d);    
       CRF.bl2 = Make_Blade(CRF.g.chi(:,1,2),CRF.g.chi(:,2,2),CRF.g.r,CRF.g.c(:,2),2,2,d);
       txt = {'Blades:', CRF.g.N'};
       text(0.01,0,txt)
   end
else
   disp('Blade geometry cannot be calculated for multiple designs (yet)');
end

%% ---- 5. Save blade designs to file ---- %%

% % Save both types of design, rename structures to "d"
% d = NRF; save([dr.geom 'NRF_' NRF.d.name '.mat'],'d')
% save([dr.geom 'NRFibl_' NRF.d.name '.ibl'],'d')
% d = CRF; save([dr.geom 'CRF_' CRF.d.name '.mat'],'d')

%% Rotor Stator Interaction for Noise - WORK IN PROGRESS!
% nR = repmat([1:5],5,1);
% nS = repmat([1:5]',1,5);
% % Cut on frequency at tip
% wcuton = (atm.T*287*1.4)^0.5 / NRF.g.rcas; 
% R = 11;
% S = 15;
% wpattern = NRF.d.rpm1 .* nR .* R ./ (nR .* R - nS .* S) / 60 * 2 * pi / wcuton;
% mpattern = nR .* R - nS .* S;
%NRF.g.r = linspace(0.2,0.6,31)'

