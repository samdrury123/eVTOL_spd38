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
% eVTOL;

% Run 27 design cube
nn = 3;
phi = repmat(reshape([0.5 0.7 0.9],[nn 1 1]),[1 nn nn]);
sig = repmat(reshape([0.8 1 1.2],[1 nn 1]),[nn 1 nn]);
r_cas = repmat(reshape([0.5 0.7 0.9],[1 1 nn]),[nn nn 1]);


%% -- 1. INPUTS -- %%

% ::: 1.1 Design Name ::: %
d.name = '1';

% ::: 1.1 Design ::: %
d.phi = 0.5; % 0.285; % Flow coefficient, set to zero for array of designs 
%                       velocity triangles do not work for this option (yet)
d.sig = 1.3; %0.5; % Area ratio of exit nozzle, sig = A4/A1. Eq (1) in .ppt

if d.phi==0
   phi1d = (0.2:0.1:1)';
   sig1d = (0.6:0.1:1.4)';
   d.phi = repmat(phi1d',[size(sig1d),1]);
   d.sig = repmat(sig1d,[1,size(phi1d)]);
end

% Split d into take-off and cruise
dTO=d; dcru=d;

% ::: 1.2 Performance requirements ::: %
dTO.T = 2000; % Thrust required in N
dTO.u0 = 51.4; % 50 kts
dTO.alt = 1524;
[dTO.atm.T,dTO.atm.p,dTO.atm.ro,~] = Altitude(dTO.alt); % At sea level

dcru.T = 5500; %327; % Thrust for forward flight in N
dcru.u0 = 24.7; % Forward flight speed in m/s
dcru.alt = 0;
[dcru.atm.T,dcru.atm.p,dcru.atm.ro,~] = Altitude(dcru.alt); % At 5000ft

% ::: 1.3 Geometry ::: %
g.rcas = 0.9; % casing radius in m
g.htrat = 0.3; % hub to tip ratio

% ::: 1.4 Atmosphere ::: %
% Remove and replaced with Altitude function K

% ::: 1.5 Efficiency and loading ::: %
eta.aero = 0.88; % efficiency of aerodynamic components, i.e. shaft to fluid
eta.elec = 0.85; % efficiency of electrical components, i.e. battery to shaft

% ::: 1.6 3D Blade Design ::: %
dev = 1; % d is whether to design with constant deviation across span 
%        - d = 0: chord is set to given constant diffusion factor across span
%        - d = 1:  chord is set to give constant deviation across span based on midspan value
% NRF.AR = 2.5; % Aspect ratio of blades, span to meanspan-chord
% NRF.DF = 0.35; % Lieblein Diffusion Factor
% NRF.n = 0.5; % exponent for vortex design of conventional fan, n=1 is free, n=-1 is forced
% NRF.blades = [0 0]; % Manually set number of blades:
                    % - If set to 0, meanspan chord set by AR, pitch-to-chord by DF and blade count by meanspan pitch
                    % - If set to a number, pitch is calculated from blade count and chord from pitch-to-chord from DF 
% NRF.cman = %[0.8,1.5;2,1]; % Scale calculated chord distribution [hub B1,cas B1;hub B2, casB2]. 


CRF.AR = 3;
CRF.DF = 0.35;
CRF.n = 0.5; % exponent for vortex design of contra-rotating fan, n=1 is free, n=-1 is forced
CRF.blades = [0 0]; % Manually set number of blades, as above for NRF
% CRF.cman = [1.4,1;1.2,1]; % Scale calculated chord distribution [hub B1,cas B1;hub B2, casB2]. 
  
%% ---- 2. CONTROL VOLUME ANALYSIS ---- %%
% ::: 2.1 Take Off ::: %
[dTO,g,eta_TO] = CVanalysis(dTO,g,eta);

% ::: 2.2 Cruise ::: %
[dcru,g,eta_cru] = CVanalysis(dcru,g,eta);
NRF.dcru = dcru;
CRF.dcru = dcru;

% ::: 2.3 Plot design performance maps ::: %
if size(d.phi) > 1
    % Hover
    figure; set(gcf, 'color', 'w'); grid off; box on; set(gcf,'Position',[20 200 1000 400]);
    subplot(2,4,1); contourf(d.phi,d.sig,dTO.Pideal,20); colorbar
    xlabel('phi'), ylabel('sigma'); title('Hover Ideal Power (kW)')
    text(0.25,0.45,'Halve for each CRF row','color','w')
        
    subplot(2,4,2); contourf(d.phi,d.sig,dTO.rpm,20); colorbar
    xlabel('phi'), ylabel('sigma'); title('Hover Shaft speed (rpm)')
    text(0.25,0.45,'Halve for each CRF row','color','w')
    
    subplot(2,4,3); contourf(d.phi,d.sig,dTO.Tq,20); colorbar
    xlabel('phi'), ylabel('sigma'); title('Hover Torque (Nm)')
    
    subplot(2,4,4); contourf(d.phi,d.sig,dTO.psitt,20); colorbar
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
    
%     subplot(2,3,6); contourf(d.phi,d.sig,eta_cru.prop,20); colorbar
%     xlabel('phi'), ylabel('sigma'); title('Propulsive Efficiency')   
     subplot(2,4,8); contourf(d.phi,d.sig,dcru.psitt,20); colorbar
     xlabel('phi'), ylabel('sigma'); title('Cruise Stage Loading Coefficient')   


end

%% ---- 3. Velocity Triangles ---- %%
if size(d.phi) == 1   
   % ::: 3.1 Conventional compressor with rotor-stator ::: %
%    k = 1; % k is Blade A proportion of rpm, k = 1 for conventional    
%    [NRF.a,NRF.g,NRF.d] = veltriangles(dTO,g,k,NRF,dev,1);
   
   % ::: 3.2 Counter-rotating compressor with rotor-rotor ::: %
   k=0.5; % k = 0.5 for counter-rotating 50-50 split between Rotor A and B 
   [CRF.a,CRF.g,CRF.d] = veltriangles(dTO,g,k,CRF,dev,1);
   
else
   disp('Velocity triangles cannot be calculated for multiple designs (yet)');
end

%% ---- 4. Use Mafiomoli Sections to Plot Blades ---- %%
if size(d.phi) == 1 
%    NRF.bl1 = Make_Blade(NRF.g.chi(:,1,1),NRF.g.chi(:,2,1),NRF.g.r,NRF.g.c(:,1),1,1);    
%    NRF.bl2 = Make_Blade(NRF.g.chi(:,1,2),NRF.g.chi(:,2,2),NRF.g.r,NRF.g.c(:,2),2,1);    

   CRF.bl1 = Make_Blade(CRF.g.chi(:,1,1),CRF.g.chi(:,2,1),CRF.g.r,CRF.g.c(:,1),1,1);    
   CRF.bl2 = Make_Blade(CRF.g.chi(:,1,2),CRF.g.chi(:,2,2),CRF.g.r,CRF.g.c(:,2),2,1);    
else
   disp('Blade geometry cannot be calculated for multiple designs (yet)');
end

%% ---- 5. Save blade designs to file ---- %%

% % Save both types of design, rename structures to "d"
% d = NRF; save([dr.geom 'NRF_' NRF.d.name '.mat'],'d')
% d = CRF; save([dr.geom 'CRF_' CRF.d.name '.mat'],'d')

%% ---- 6. Write out inputs and outputs for spreadsheet ---- %%
% if size(d.phi) == 1 
%    display(num2str(dTO.T))
%    display(num2str(dTO.u0))
%    display(num2str(dTO.alt))
%    display(num2str(g.rcas*2))
%    display(num2str(g.htrat))
%    display(num2str(dTO.sig))
%    display(num2str(dTO.phi))
%    display(' ')
%    display(num2str(dTO.Pideal/2))
%    display(num2str(dTO.rpm/2))
%    display(num2str(dTO.Tq))
%    display(num2str(dTO.phi))
%    display(num2str(dTO.psitt))
%    display(num2str(dTO.Mtip/2))
%    display(num2str(dTO.Mtip))
%    display(num2str(dTO.Mjet))
%    display(num2str(eta_TO.prop))
%    display(num2str(eta_TO.mf))
%    display(' ')
%    display(num2str(dcru.T))
%    display(num2str(dcru.u0))
%    display(num2str(dcru.alt))
%    display(num2str(g.rcas*2))
%    display(num2str(g.htrat))
%    display(num2str(dcru.sig))
%    display(num2str(dcru.phi))
%    display(' ')
%    display(num2str(dcru.Pideal/2))
%    display(num2str(dcru.rpm/2))
%    display(num2str(dcru.Tq))
%    display(num2str(dcru.phicru))
%    display(num2str(dcru.psitt))
%    display(num2str(dcru.Mtip/2))
%    display(num2str(dcru.Mtip))
%    display(num2str(dcru.Mjet))
%    display(num2str(eta_cru.prop))
%    display(num2str(eta_cru.mf))
%    
% end

%% Rotor Stator Interaction - WORK IN PROGRESS!
% nR = repmat([1:5],5,1);
% nS = repmat([1:5]',1,5);
% % Cut on frequency at tip
% wcuton = (atm.T*287*1.4)^0.5 / NRF.g.rcas; 
% R = 11;
% S = 15;
% wpattern = NRF.d.rpm1 .* nR .* R ./ (nR .* R - nS .* S) / 60 * 2 * pi / wcuton;
% mpattern = nR .* R - nS .* S;



