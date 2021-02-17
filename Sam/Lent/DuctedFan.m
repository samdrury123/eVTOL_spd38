clear; close all;

%% SPD Radial Eq 08/02
bluebear = 1;
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
% Make trailing edge thickness absolute

%% To be added in this branch:
% DONE Vortex exponent
% DONE Radial Equilibrium
% DONELook at removing CRDF flag
% DONE Change tmax, tTe to be percentages of chord
% DONE Improve variable naming for velocities - make location an index
% DONE Add loss charts
% DONE Rename rloss, sloss to b1, b2 
% DONE Rename Nr, Ns to N1, N2
% ERROR HANDLING
% DONE Merge VelTriangles and BladeLoss
% (not needed) Merge a and q structs?
% DONE Tidy up names for functions
% DONE Vary T,p,rho spanwise

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
    g.tTE = 2/1000; % Trailing edge thickness in m
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
    g.tTE = 1/1000; % Trailing edge thickness in m
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
% philist = [0.7 0.8];
% sigmalist=1.2;
% philist=.8;

% Initialise structs, plotting arrays
design([ size(philist,2) size(sigmalist,2) ]) = struct();
% NRFdesign = design;
% CRFdesign = design;

% Arrays for plotting - filling with NaNs ensures duff data doesn't skew contour plots
init = NaN.*ones(size(philist,2), size(sigmalist,2)); 
N.phis = init;
N.sigmas = init;
N.etas = init;  % Fan efficiency 
N.Rs = init; % Reaction
N.rpms = init; % rpm
N.N1 = init; % Blade count row 1
N.N2 = init; % Blade count row 2
N.Frs = init; % Propulsive efficiency
N.psis = init;
N.FOMs = init; % Figure of Merit
N.Ltot = init;
N.L_eta = init;
N.profL = init;
N.baseL = init;
N.tipL = init;
N.endwallL = init;
N.DH1 = init; N.DH2 = init;
N.name = 'NRF';
C=N; C.name = 'CRF';

%% Begin loop

for pp = 1:size(philist,2)
    for ss = 1:size(sigmalist,2)     
tic

d.phi = philist(pp);
d.sigma = sigmalist(ss);

%% Iteation loop for NRF losses

L.L1 = 0;     % Entropy change in rotor
L.L2 = 0;     % Entropy change in stator
L.deltaL1 = 1;
L.deltaL2 = 1;

while abs(L.deltaL1) > 0.001/100 && abs(L.deltaL2) > 0.001/100
% Iteration loop to find T03 and T04, again improve so not absolute value
[d,g,q] = CVanalysis_comp(d,g,q,L);
% Meanline velocity triangles
k=1; %(NRF), k = % of blade speed in first row
dev = (1-bluebear); % Deviation constant across the span (1) or DF constant (0)
[d,g,a,q,L] = VelTriangles(d,g,q,L,k,dev);
% Evaluate analytical velocity profiles and loss for each blade row, plus calculate shroud clearance and endwall losses
% [g,L] = BladeLoss(d,g,a,q,L,k); % Return Re?
end
NRF.d=d; NRF.g=g; NRF.a=a; NRF.q=q; NRF.L=L;

% Finished loss loop so everything defined. The rest is gathering together
% outputs
% Hacking together arrays to plot
if ~isnan(g.Nb(:))
    N.phis(pp,ss) = d.phi;
    N.sigmas(pp,ss) = d.sigma;
    N.etas(pp,ss) = d.eta;  % Fan efficiency 
%     N.Rs(pp,ss) = d.Reaction; % Reaction
    N.rpms(pp,ss) = d.rpm1; % rpm
    N.N1(pp,ss) = g.Nb(1); % Blade count row 1
    N.N2(pp,ss) = g.Nb(2); % Blade count row 2
    N.Frs(pp,ss) = d.Fr; % Propulsive efficiency
    N.psis(pp,ss) = d.psi;
    N.FOMs(pp,ss) = d.Mf;
    Ltot = L.L1 + L.L2; 
    N.Ltot(pp,ss) = Ltot;
    N.L_eta(pp,ss) = L.eta_s;
    N.profL(pp,ss) = (L.Loss1.prof + L.Loss2.prof)/Ltot;
    N.baseL(pp,ss) = (L.Loss1.base + L.Loss2.base)/Ltot;
    N.tipL(pp,ss) = (L.Loss1.tip + L.Loss2.tip)/Ltot;
    N.endwallL(pp,ss) = (L.Loss1.endwall + L.Loss2.endwall)/Ltot;
    N.DH1(pp,ss) = d.DH(1);
    N.DH2(pp,ss) = d.DH(2);
end
%% Iteation loop for CRF losses

L.L1 = 0;     % Entropy change in rotor
L.L2 = 0;     % Entropy change in stator
L.deltaL1 = 1;
L.deltaL2 = 1;
% Iteation loop for CRF losses
while abs(L.deltaL1) > 0.001/100 && abs(L.deltaL2) > 0.001/100
% Iteration loop to find T03 and T04, again improve so not absolute value
[d,g,q] = CVanalysis_comp(d,g,q,L);
% Meanline velocity triangles
k=0.5; %(CRF), k = % of blade speed in first row
dev = (1-bluebear); % Deviation constant across the span (1) or DF constant (0)
[d,g,a,q,L] = VelTriangles(d,g,q,L,k,dev);
% Evaluate analytical velocity profiles and loss for each blade row, plus calculate shroud clearance and endwall losses
% [g,L] = BladeLoss(d,g,a,q,L,k); % Return Re?
end
CRF.d=d; CRF.g=g; CRF.a=a; CRF.q=q; CRF.L=L;


% Finished loss loop so everything defined. The rest is gathering together
% outputs
% Hacking together arrays to plot

if ~isnan(g.Nb(:))
    C.phis(pp,ss) = d.phi;
    C.sigmas(pp,ss) = d.sigma;
    C.etas(pp,ss) = d.eta;  % Fan efficiency 
%     C.Rs(pp,ss) = d.Reaction; % Reaction
    C.rpms(pp,ss) = d.rpm1; % rpm
    C.N1(pp,ss) = g.Nb(1); % Blade count row 1
    C.N2(pp,ss) = g.Nb(2); % Blade count row 2
    C.Frs(pp,ss) = d.Fr; % Propulsive efficiency
    C.psis(pp,ss) = d.psi;
    C.FOMs(pp,ss) = d.Mf;
    Ltot = L.L1 + L.L2; 
    C.Ltot(pp,ss) = Ltot;
    C.L_eta(pp,ss) = L.eta_s;
    C.profL(pp,ss) = (L.Loss1.prof + L.Loss2.prof)/Ltot;
    C.baseL(pp,ss) = (L.Loss1.base + L.Loss2.base)/Ltot;
    C.tipL(pp,ss) = (L.Loss1.tip + L.Loss2.tip)/Ltot;
    C.endwallL(pp,ss) = (L.Loss1.endwall + L.Loss2.endwall)/Ltot;
    C.DH1(pp,ss) = d.DH(1);
    C.DH2(pp,ss) = d.DH(2);
end

%% Store outputs in structures
design(pp,ss).phisig = [d.phi, d.sigma];
design(pp,ss).NRF = NRF; design(pp,ss).NRF.phisig = [d.phi, d.sigma];
design(pp,ss).CRF = CRF; design(pp,ss).CRF.phisig = [d.phi, d.sigma];

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
% lots of small or one big
% thrusts for same op point
% compare different methods of op point matching
% MISMATCH OF cruise and hover op points - fundamental challenge

% Summarise
display(['phi = ' num2str(d.phi) '  sigma = ' num2str(d.sigma) '   Time = ' num2str(toc,3) ]);
    end
end

%% Plotting
% Limits used for plotting, ensures effective comparisons
lims.Nb = [0 25]; % Sets limit of colormap for blade numbers
lims.L = [0 0.5]; % Sets loss limits between 0 and 1
lims.Ltot = [min(min(N.Ltot(:)), min(C.Ltot(:))) max(max(N.Ltot(:)), max(C.Ltot(:)))];

% Need to improve de Haller plotting
PlotCharts(N,lims);
% PlotCharts(C,lims);
% PlotVels(design(24),1)
% PlotVels(blah)
% velocity triangles, profiles (hub tip mean), chord
