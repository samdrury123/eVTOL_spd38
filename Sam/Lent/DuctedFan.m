clear; close all;

%% SPD Radial Eq 08/02
bluebear = 0;
CRDF_flag = 0; % Change this maybe?

%% Still to add:
% Vortex exponent
% Radial equilibrium
% Couple these with the blade profile losses
% Check what needs to be sent to CAD
% Improve variable naming for velocities - make location an index

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
    g.tTE = 0.003; % Trailing edge thickness, this should be a % of chord
    g.tmax = 0.02; % Max blade thickness, this should be a % of chord
    g.gap = 0.001; % Shroud clearance in m
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
    g.tTE = 0.0005; % Trailing edge thickness, in m
    g.tmax = 0.003; % Max blade thickness, in m
    g.gap = 0.0005; % Shroud clearance in m
end

% Gas Constants
gam = 1.4; 
R = 287;
cp = gam*R / (gam-1);
q.gam=gam; q.R=R; q.cp=cp;

sigmalist = [0.8 0.9 1 1.1 1.2]; 
sigmalist = [0.7 0.8 0.9 1 1.1 1.2]; 
% sigmalist = [0.9 1 1.1];
% philist = [0.5 0.6 0.7 0.8 0.9];
philist = [0.5 0.6 0.7 0.8];
% sigmalist=1.2;
% philist=.8;

design([ size(philist,2) size(sigmalist,2) ]) = struct(); % Initialise

%% Begin loop

for pp = 1:size(philist,2)
    for ss = 1:size(sigmalist,2)
    tic

d.phi = philist(pp);
d.sigma = sigmalist(ss);

L.Lr = 0;     % Entropy change in rotor
L.Ls = 0;     % Entropy change in stator
L.deltaLr = 1;
L.deltaLs = 1;

% Iteation loop for losses, improve while criterion (i.e. % change)
while abs(L.deltaLr) > 0.001/100 && abs(L.deltaLs) > 0.01/100
%% Iteration loop to find T03 and T04, again improve so not absolute value
[d,g,q] = CVanalysis_comp(d,g,q,L, CRDF_flag);
%% Meanline velocity triangles
[a,q] = VelTriangles(d,q,L,CRDF_flag);
%% Evaluate analytical velocity profiles and loss for each blade row, plus calculate shroud clearance and endwall losses
[g,L] = BladeLoss(d,g,a,q,L); % Return Re?
end

%% Summarise
display(['phi = ' num2str(d.phi) '  sigma = ' num2str(d.sigma) '   Time = ' num2str(toc,3) ]);


% Finished loss loop so everything defined. The rest is gathering together
% outputs

% Power
mdot = q.Mndp01 * g.A1 * q.p01 / (cp*q.T01)^0.5;
P = mdot * cp*(q.T03-q.T01);

% Reaction
Reaction = (q.T2-q.T1)/(q.T3-q.T1);

% Tip Mach Number
Mtip2 = (d.rpm/60*2*pi*g.rc) / (gam*R*q.T1)^0.5;

% Propulsive efficiency and hover figure of merit
Fr = 2*d.u0 / (d.u0+d.u4);
Mf = d.Th / P * sqrt(d.Th / (2*q.atm.ro*g.A1));

% Fan Efficiency
eta = (q.T03s - q.T01) / (q.T03 - q.T01);

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


% Hacking together arrays to plot
phis(pp,ss) = d.phi;
sigmas(pp,ss) = d.sigma;
etas(pp,ss) = eta;  % Fan efficiency 
Rs(pp,ss) = Reaction; % Reaction
rpms(pp,ss) = d.rpm; % rpm
Nr(pp,ss) = g.Nb(1); % Blade count row 1
Ns(pp,ss) = g.Nb(2); % Blade count row 2
Frs(pp,ss) = Fr; % Propulsive efficiency
psis(pp,ss) = d.psi;
FOMs(pp,ss) = Mf;

% Store outputs in structures
design(pp,ss).phisig = [d.phi, d.sigma];
design(pp,ss).d = d;
design(pp,ss).g = g;
design(pp,ss).a = a;
design(pp,ss).q = q;
design(pp,ss).L = L;


    
    end
end
%% Plotting

figure(1);set(gcf, 'color', 'w'); grid off; box on;
set(gcf,'Position',[20 50 1300 600]);
subplot(2,4,1)
contourf(phis,sigmas,etas,25, 'edgecolor','none'); colorbar
title("Fan efficiency \eta_a"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,2)
contourf(phis,sigmas,Frs,25, 'edgecolor','none'); colorbar
title("Propulsive efficiency \eta_p"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,3)
contourf(phis,sigmas,etas.*Frs,25, 'edgecolor','none'); colorbar
title("Overall efficiency \eta_{ov}"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,4)
contourf(phis,sigmas,FOMs,25, 'edgecolor','none'); colorbar
title("Figure of Merit"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,5)
contourf(phis,sigmas,rpms,25, 'edgecolor','none'); colorbar
title("RPM"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,6)
contourf(phis,sigmas,Nr,25, 'edgecolor','none'); colorbar
title("Rotor blades"); xlabel("\phi"); ylabel("\sigma");
text(1.1,1.3, ['CRDF = ' num2str(CRDF_flag)],'color','b')
subplot(2,4,7)
contourf(phis,sigmas,Ns,25, 'edgecolor','none'); colorbar
title("Stator blades"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,8)
contourf(phis,sigmas,psis, 'edgecolor','none'); colorbar
title("Psi \psi"); xlabel("\phi"); ylabel("\sigma");


figure(2);set(gcf, 'color', 'w'); grid off; box on;
set(gcf,'Position',[20 50 1000 600]);
subplot(2,3,1)
contourf(phis,psis,etas,25, 'edgecolor','none'); colorbar
title("Fan efficiency \eta_a"); xlabel("\phi"); ylabel("\psi");
subplot(2,3,2)
contourf(phis,psis,Frs,25, 'edgecolor','none'); colorbar
title("Propulsive efficiency \eta_p"); xlabel("\phi"); ylabel("\psi");
if d.u0 == 0
    subplot(2,3,3)
    contourf(phis,psis,FOMs,25, 'edgecolor','none'); colorbar
    title("Figure of Merit"); xlabel("\phi"); ylabel("\psi");
else
    subplot(2,3,3)
    contourf(phis,psis,etas.*Frs,25, 'edgecolor','none'); colorbar
    title("Overall efficiency \eta_{ov}"); xlabel("\phi"); ylabel("\psi");
end
subplot(2,3,4)
contourf(phis,psis,rpms,25, 'edgecolor','none'); colorbar
title("RPM"); xlabel("\phi"); ylabel("\psi");
text(1.1,0.8, ['CRDF = ' num2str(CRDF_flag)],'color','b')
subplot(2,3,5)
contourf(phis,psis,Nr,25, 'edgecolor','none'); colorbar
title("Rotor blades"); xlabel("\phi"); ylabel("\psi");
subplot(2,3,6)
contourf(phis,psis,Ns,25, 'edgecolor','none'); colorbar
title("Stator blades"); xlabel("\phi"); ylabel("\psi");





