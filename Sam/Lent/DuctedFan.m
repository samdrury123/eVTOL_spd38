clear; close all;

%% SPD Radial Eq 08/02
bluebear = 0;
% CRDF_flag = 0; % Change this maybe?

%% Still to add:
% Vortex exponent
% Radial equilibrium
% Couple these with the blade profile losses
% Check what needs to be sent to CAD
% Improve variable naming for velocities - make location an index
% Should gap be percentage of span?
% Exit swirl
% ERROR HANDLING
% Allow NRF and CRF structs
% Add efficiencies - two motors has greater inefficiencies 

%% To be added in this branch:
% DONE Vortex exponent
% DONE Radial Equilibrium
% DONELook at removing CRDF flag
% DONE Change tmax, tTe to be percentages of chord
% Improve variable naming for velocities - make location an index
% Add loss charts
% Rename rloss, sloss to b1, b2 for contra comaptability
% ERROR HANDLING
% Merge VelTriangles and BladeLoss
% Merge a and q structs?

%% Input variables

if bluebear == 1
    d.name = 'Blue Bear';
    d.Th = 1500; % Thrust
    d.u0 = 2*25.7; % Flight speed
    d.alt = 0; % Altitude in m
    d.Cdnom = 0.002; % Nominal dissipation coefficient for BL loss    
    d.Cc = 0.6; % Contraction coefficient for shroud loss calculation
    d.DF = 0.35; % Diffusion factor for Lieblein
    
    q.atm = Altitude(d.alt); % At sea level
    
    g.D = 1.4; % Casing diameter
    g.htrat = 0.3; % Hub-to-tip ratio
    g.AR = 4; % Blade aspect ratio
    g.n = 0.5; % Exponent for vortex design of conventional fan, n=1 is free, n=-1 is forced
    g.tTE = 0.025; % Trailing edge thickness, % of chord (ie 2.5%)
    g.tmax = 0.163; % Max blade thickness, % of chord
    g.gap = 0.001; % Shroud clearance in m
else % Whittle eVTOL
    d.name = 'The Rt.Hon Whittle eVTOL';
    d.Th = 10; % Thrust
    d.u0 = 0; % Flight speed
    d.alt = 0; % Altitude in m
    d.Cdnom = 0.0002; % Nominal dissipation coefficient for BL loss (Dickens p154)
    d.Cc = 0.6; % Contraction coefficient for shroud loss calculation - Yoon p5 quotes Cd = 0.7
    d.DF = 0.35; % Diffusion factor for Lieblein
    
    q.atm = Altitude(d.alt); % At sea level [Tatm,patm,roatm,~]
    
    g.D = 0.06*2; % Casing diameter
    g.htrat = 1/3; % Hub-to-tip ratio
    g.AR = 1.8; % Blade aspect ratio
    g.n = 0.5; % Exponent for vortex design of conventional fan, n=1 is free, n=-1 is forced
    g.tTE = 0.025; % Trailing edge thickness, % of chord (ie 2.5%)
    g.tmax = 0.15; % Max blade thickness, % of chord
    g.gap = 0.0005; % Shroud clearance in m %% CHANGE THIS TO % OF CHORD??
end

% Gas Constants
gam = 1.4; 
R = 287;
cp = gam*R / (gam-1);
q.gam=gam; q.R=R; q.cp=cp;

sigmalist = [0.8 0.9 1 1.1 1.2]; 
sigmalist = [0.7 0.8 0.9 1 1.1 1.2]; 
% sigmalist = [0.9 1 1.1];
philist = [0.5 0.6 0.7 0.8 0.9];
% philist = [0.5 0.6 0.7 0.8];
% sigmalist=1.2;
% philist=.8;

% Initialise structs, plotting arrays
design([ size(philist,2) size(sigmalist,2) ]) = struct();
% NRFdesign = design;
% CRFdesign = design;

% Filling with NaNs ensures duff data doesn't skew contour plots
init = NaN.*ones(size(philist,2), size(sigmalist,2)); 
N.phis = init;
N.sigmas = init;
N.etas = init;  % Fan efficiency 
N.Rs = init; % Reaction
N.rpms = init; % rpm
N.Nr = init; % Blade count row 1
N.Ns = init; % Blade count row 2
N.Frs = init; % Propulsive efficiency
N.psis = init;
N.FOMs = init; % Figure of Merit
C=N;

%% Begin loop

for pp = 1:size(philist,2)
    for ss = 1:size(sigmalist,2)     
tic

d.phi = philist(pp);
d.sigma = sigmalist(ss);

%% Iteation loop for NRF losses

L.Lr = 0;     % Entropy change in rotor
L.Ls = 0;     % Entropy change in stator
L.deltaLr = 1;
L.deltaLs = 1;

while abs(L.deltaLr) > 0.001/100 && abs(L.deltaLs) > 0.001/100
% Iteration loop to find T03 and T04, again improve so not absolute value
[d,g,q] = CVanalysis_comp(d,g,q,L);
% Meanline velocity triangles
k=1; %(NRF), k = % of blade speed in first row
[d,g,a,q] = VelTriangles(d,g,q,L,k);
% Evaluate analytical velocity profiles and loss for each blade row, plus calculate shroud clearance and endwall losses
[g,L] = BladeLoss(d,g,a,q,L); % Return Re?
end
NRF.d=d; NRF.g=g; NRF.a=a; NRF.q=q; NRF.L=L;

% Finished loss loop so everything defined. The rest is gathering together
% outputs
% Hacking together arrays to plot
if ~isnan(g.Nb(:))
    N.phis(pp,ss) = d.phi;
    N.sigmas(pp,ss) = d.sigma;
    N.etas(pp,ss) = d.eta;  % Fan efficiency 
    N.Rs(pp,ss) = d.Reaction; % Reaction
    N.rpms(pp,ss) = d.rpm1; % rpm
    N.Nr(pp,ss) = g.Nb(1); % Blade count row 1
    N.Ns(pp,ss) = g.Nb(2); % Blade count row 2
    N.Frs(pp,ss) = d.Fr; % Propulsive efficiency
    N.psis(pp,ss) = d.psi;
    N.FOMs(pp,ss) = d.Mf;
end
%% Iteation loop for CRF losses

L.Lr = 0;     % Entropy change in rotor
L.Ls = 0;     % Entropy change in stator
L.deltaLr = 1;
L.deltaLs = 1;
% Iteation loop for CRF losses
while abs(L.deltaLr) > 0.001/100 && abs(L.deltaLs) > 0.001/100
% Iteration loop to find T03 and T04, again improve so not absolute value
[d,g,q] = CVanalysis_comp(d,g,q,L);
% Meanline velocity triangles
k=0.5; %(CRF), k = % of blade speed in first row
[d,g,a,q] = VelTriangles(d,g,q,L,k);
% Evaluate analytical velocity profiles and loss for each blade row, plus calculate shroud clearance and endwall losses
[g,L] = BladeLoss(d,g,a,q,L); % Return Re?
end
CRF.d=d; CRF.g=g; CRF.a=a; CRF.q=q; CRF.L=L;

% Finished loss loop so everything defined. The rest is gathering together
% outputs
% Hacking together arrays to plot
if ~isnan(g.Nb(:))
    C.phis(pp,ss) = d.phi;
    C.sigmas(pp,ss) = d.sigma;
    C.etas(pp,ss) = d.eta;  % Fan efficiency 
    C.Rs(pp,ss) = d.Reaction; % Reaction
    C.rpms(pp,ss) = d.rpm1; % rpm
    C.Nr(pp,ss) = g.Nb(1); % Blade count row 1
    C.Ns(pp,ss) = g.Nb(2); % Blade count row 2
    C.Frs(pp,ss) = d.Fr; % Propulsive efficiency
    C.psis(pp,ss) = d.psi;
    C.FOMs(pp,ss) = d.Mf;
end

% Store outputs in structures
design(pp,ss).phisig = [d.phi, d.sigma];
design(pp,ss).NRF = NRF;
design(pp,ss).CRF = CRF;

% % Inlet streamtube capture 
% % Currently assumes isentropic inlet duct
% Mndp00 = gam / (gam-1)^0.5 .* M0 .* (1 + (gam-1)/2.*M0.^2).^(-0.5*(gam+1)/(gam-1));
% A0 = A1 * Mndp01 / Mndp00;
% CAR = A0/A1;
% rnrat = 1.2;
% Ain = pi*(rc*rnrat)^2;
% MFR = A0/Ain;

% Advance ratio and loading coefficient
% J = u0 / (rpm/60 * 2*rc);
% CT = Th / (roatm * (rpm/60)^2 * (2*rc)^4);

% Summarise
display(['phi = ' num2str(d.phi) '  sigma = ' num2str(d.sigma) '   Time = ' num2str(toc,3) ]);
    end
end

%% NRF Plotting
Nb_lim = [0 20]; % Sets limit of colormap for blade numbers
phis=N.phis; sigmas=N.sigmas; psis=N.psis;

figure(1);set(gcf, 'color', 'w'); grid off; box on;
set(gcf,'Position',[20 50 1300 600]);
subplot(2,4,1)
contourf(phis,sigmas,N.etas,25, 'edgecolor','none'); colorbar
title("Fan efficiency \eta_a"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,2)
contourf(phis,sigmas,N.Frs,25, 'edgecolor','none'); colorbar
title("Propulsive efficiency \eta_p"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,3)
contourf(phis,sigmas,N.etas.*N.Frs,25, 'edgecolor','none'); colorbar
title("Overall efficiency \eta_{ov}"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,4)
contourf(phis,sigmas,N.FOMs,25, 'edgecolor','none'); colorbar
title("Figure of Merit"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,5)
contourf(phis,sigmas,N.rpms,25, 'edgecolor','none'); colorbar
title("RPM"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,6)
contourf(phis,sigmas,N.Nr,25, 'edgecolor','none'); colorbar; caxis(Nb_lim);
title("Rotor blades"); xlabel("\phi"); ylabel("\sigma");
% text(1.1,1.3, ['CRDF = ' num2str(CRDF_flag)],'color','b')
subplot(2,4,7)
contourf(phis,sigmas,N.Ns,25, 'edgecolor','none'); colorbar; caxis(Nb_lim);
title("Stator blades"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,8)
contourf(phis,sigmas,N.psis,25, 'edgecolor','none'); colorbar
title("Psi \psi"); xlabel("\phi"); ylabel("\sigma");


figure(2);set(gcf, 'color', 'w'); grid off; box on;
set(gcf,'Position',[20 50 1000 600]);
subplot(2,3,1)
contourf(phis,psis,N.etas,25, 'edgecolor','none'); colorbar
title("Fan efficiency \eta_a"); xlabel("\phi"); ylabel("\psi");
subplot(2,3,2)
contourf(phis,psis,N.Frs,25, 'edgecolor','none'); colorbar
title("Propulsive efficiency \eta_p"); xlabel("\phi"); ylabel("\psi");
if d.u0 == 0
    subplot(2,3,3)
    contourf(phis,psis,N.FOMs,25, 'edgecolor','none'); colorbar
    title("Figure of Merit"); xlabel("\phi"); ylabel("\psi");
else
    subplot(2,3,3)
    contourf(phis,psis,N.etas.*N.Frs,25, 'edgecolor','none'); colorbar
    title("Overall efficiency \eta_{ov}"); xlabel("\phi"); ylabel("\psi");
end
subplot(2,3,4)
contourf(phis,psis,N.rpms,25, 'edgecolor','none'); colorbar
title("RPM"); xlabel("\phi"); ylabel("\psi");
% text(1.1,0.8, ['CRDF = ' num2str(CRDF_flag)],'color','b')
subplot(2,3,5)
contourf(phis,psis,N.Nr,25, 'edgecolor','none'); colorbar; caxis(Nb_lim);
title("Rotor blades"); xlabel("\phi"); ylabel("\psi");
subplot(2,3,6)
contourf(phis,psis,N.Ns,25, 'edgecolor','none'); colorbar; caxis(Nb_lim);
title("Stator blades"); xlabel("\phi"); ylabel("\psi");

%% CRF plotting

phis=C.phis; sigmas=C.sigmas; psis=C.psis;

figure(3);set(gcf, 'color', 'w'); grid off; box on;
set(gcf,'Position',[20 50 1300 600]);
subplot(2,4,1)
contourf(phis,sigmas,C.etas,25, 'edgecolor','none'); colorbar
title("Fan efficiency \eta_a"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,2)
contourf(phis,sigmas,C.Frs,25, 'edgecolor','none'); colorbar
title("Propulsive efficiency \eta_p"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,3)
contourf(phis,sigmas,C.etas.*C.Frs,25, 'edgecolor','none'); colorbar
title("Overall efficiency \eta_{ov}"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,4)
contourf(phis,sigmas,C.FOMs,25, 'edgecolor','none'); colorbar
title("Figure of Merit"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,5)
contourf(phis,sigmas,C.rpms,25, 'edgecolor','none'); colorbar
title("RPM"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,6)
contourf(phis,sigmas,C.Nr,25, 'edgecolor','none'); colorbar; caxis(Nb_lim);
title("Rotor blades"); xlabel("\phi"); ylabel("\sigma");
% text(1.1,1.3, ['CRDF = ' num2str(CRDF_flag)],'color','b')
subplot(2,4,7)
contourf(phis,sigmas,C.Ns,25, 'edgecolor','none'); colorbar; caxis(Nb_lim);
title("Stator blades"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,8)
contourf(phis,sigmas,C.psis,25, 'edgecolor','none'); colorbar
title("Psi \psi"); xlabel("\phi"); ylabel("\sigma");


figure(4);set(gcf, 'color', 'w'); grid off; box on;
set(gcf,'Position',[20 50 1000 600]);
subplot(2,3,1)
contourf(phis,psis,C.etas,25, 'edgecolor','none'); colorbar
title("Fan efficiency \eta_a"); xlabel("\phi"); ylabel("\psi");
subplot(2,3,2)
contourf(phis,psis,C.Frs,25, 'edgecolor','none'); colorbar
title("Propulsive efficiency \eta_p"); xlabel("\phi"); ylabel("\psi");
if d.u0 == 0
    subplot(2,3,3)
    contourf(phis,psis,C.FOMs,25, 'edgecolor','none'); colorbar
    title("Figure of Merit"); xlabel("\phi"); ylabel("\psi");
else
    subplot(2,3,3)
    contourf(phis,psis,C.etas.*C.Frs,25, 'edgecolor','none'); colorbar
    title("Overall efficiency \eta_{ov}"); xlabel("\phi"); ylabel("\psi");
end
subplot(2,3,4)
contourf(phis,psis,C.rpms,25, 'edgecolor','none'); colorbar
title("RPM"); xlabel("\phi"); ylabel("\psi");
% text(1.1,0.8, ['CRDF = ' num2str(CRDF_flag)],'color','b')
subplot(2,3,5)
contourf(phis,psis,C.Nr,25, 'edgecolor','none'); colorbar; caxis(Nb_lim);
title("Rotor blades"); xlabel("\phi"); ylabel("\psi");
subplot(2,3,6)
contourf(phis,psis,C.Ns,25, 'edgecolor','none'); colorbar; caxis(Nb_lim);
title("Stator blades"); xlabel("\phi"); ylabel("\psi");

