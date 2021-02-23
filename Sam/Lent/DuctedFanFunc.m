function f = DuctedFanFunc(x)
bluebear=1;
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
    g.gap = 1/1000; % Shroud clearance in m
else % Whittle eVTOL
    d.name = 'The Rt.Hon Whittle eVTOL';
    d.Th = 10; % Thrust
    d.u0 = 0; % Flight speed
    d.alt = 0; % Altitude in m
    d.Cdnom = 0.002; % Nominal dissipation coefficient for BL loss (Dickens p154)
    d.Cc = 0.6; % Contraction coefficient for shroud loss calculation - Yoon p5 quotes Cd = 0.7
    d.DF = 0.35; % Diffusion factor for Lieblein
    
    q.atm = Altitude(d.alt); % At sea level [Tatm,patm,roatm,~]
    
    g.D = 0.06*2; % Casing diameter
    g.htrat = 1/3; % Hub-to-tip ratio
    g.AR = 1.8; % Blade aspect ratio
    g.n = 0.5; % Exponent for vortex design of conventional fan, n=1 is free, n=-1 is forced
    g.tTE = 1/1000; % Trailing edge thickness in m
    g.tmax = 0.15; % Max blade thickness, % of chord
    g.gap = 0.5/1000; % Shroud clearance in m %% CHANGE THIS TO % OF CHORD??
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
% sigmalist=[.9 1 1.1];
% philist=.9;
% These only really work for BB
% sigmalist = [0.5 0.6 0.7 0.8 0.9 1 1.1 1.2];
% philist = [0.3 0.4 0.5 0.6 0.7 0.8 0.9];

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
d.phi = x(1);
d.sigma = x(2);

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
% L.deltaL1
% L.deltaL2
end
NRF.d=d; NRF.g=g; NRF.a=a; NRF.q=q; NRF.L=L;
f = 1-d.Fr*d.eta;

% Finished loss loop so everything defined. The rest is gathering together
% outputs
% Hacking together arrays to plot
% if ~isnan(g.Nb(:))
%     N.phis(pp,ss) = d.phi;
%     N.sigmas(pp,ss) = d.sigma;
%     N.etas(pp,ss) = d.eta;  % Fan efficiency 
% %     N.Rs(pp,ss) = d.Reaction; % Reaction
%     N.rpms(pp,ss) = d.rpm1; % rpm
%     N.N1(pp,ss) = g.Nb(1); % Blade count row 1
%     N.N2(pp,ss) = g.Nb(2); % Blade count row 2
%     N.Frs(pp,ss) = d.Fr; % Propulsive efficiency
%     N.psis(pp,ss) = d.psi;
%     N.FOMs(pp,ss) = d.Mf;
%     Ltot = L.L1 + L.L2; 
%     N.Ltot(pp,ss) = Ltot;
%     N.L_eta(pp,ss) = L.eta_s;
%     N.profL(pp,ss) = (L.Loss1.prof + L.Loss2.prof)/Ltot;
%     N.baseL(pp,ss) = (L.Loss1.base + L.Loss2.base)/Ltot;
%     N.tipL(pp,ss) = (L.Loss1.tip + L.Loss2.tip)/Ltot;
%     N.endwallL(pp,ss) = (L.Loss1.endwall + L.Loss2.endwall)/Ltot;
%     N.DH1(pp,ss) = d.DH(1);
%     N.DH2(pp,ss) = d.DH(2);
% end
% %% Iteation loop for CRF losses
% 
% L.L1 = 0;     % Entropy change in rotor
% L.L2 = 0;     % Entropy change in stator
% L.deltaL1 = 1;
% L.deltaL2 = 1;
% % Iteation loop for CRF losses
% while abs(L.deltaL1) > 0.001/100 && abs(L.deltaL2) > 0.001/100
% % Iteration loop to find T03 and T04, again improve so not absolute value
% [d,g,q] = CVanalysis_comp(d,g,q,L);
% % Meanline velocity triangles
% k=0.5; %(CRF), k = % of blade speed in first row
% dev = (1-bluebear); % Deviation constant across the span (1) or DF constant (0)
% [d,g,a,q,L] = VelTriangles(d,g,q,L,k,dev);
% % Evaluate analytical velocity profiles and loss for each blade row, plus calculate shroud clearance and endwall losses
% % [g,L] = BladeLoss(d,g,a,q,L,k); % Return Re?
% end
% CRF.d=d; CRF.g=g; CRF.a=a; CRF.q=q; CRF.L=L;