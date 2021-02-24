clear; close all;

%% SPD Radial Eq 08/02
bluebear = 0;

%% Still to add:
% Check what needs to be sent to CAD
% Should gap be percentage of span?
% Exit swirl
% ERROR HANDLING
% Add efficiencies - two motors has greater inefficiencies
% On PlotVels, change x acis to [-0.1 1.1] and plot vLE as well as V?

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


%% Begin loop
% 
% for pp = 1:size(philist,2)
%     for ss = 1:size(sigmalist,2)     
% tic
options = optimset('PlotFcns',@optimplotfval);
[x,fval] = fminsearch(@DuctedFanFunc, [0.6;0.7],options)


%% Iteation loop for NRF losses

% L.L1 = 0;     % Entropy change in rotor
% L.L2 = 0;     % Entropy change in stator
% L.deltaL1 = 1;
% L.deltaL2 = 1;
% 
% while abs(L.deltaL1) > 0.001/100 && abs(L.deltaL2) > 0.001/100
% % Iteration loop to find T03 and T04, again improve so not absolute value
% [d,g,q] = CVanalysis_comp(d,g,q,L);
% % Meanline velocity triangles
% k=1; %(NRF), k = % of blade speed in first row
% dev = (1-bluebear); % Deviation constant across the span (1) or DF constant (0)
% [d,g,a,q,L] = VelTriangles(d,g,q,L,k,dev);
% % Evaluate analytical velocity profiles and loss for each blade row, plus calculate shroud clearance and endwall losses
% % [g,L] = BladeLoss(d,g,a,q,L,k); % Return Re?
% % L.deltaL1
% % L.deltaL2
% end
% NRF.d=d; NRF.g=g; NRF.a=a; NRF.q=q; NRF.L=L;
% 
% % Finished loss loop so everything defined. The rest is gathering together
% % outputs
% % Hacking together arrays to plot
% % if ~isnan(g.Nb(:))
% %     N.phis(pp,ss) = d.phi;
% %     N.sigmas(pp,ss) = d.sigma;
% %     N.etas(pp,ss) = d.eta;  % Fan efficiency 
% % %     N.Rs(pp,ss) = d.Reaction; % Reaction
% %     N.rpms(pp,ss) = d.rpm1; % rpm
% %     N.N1(pp,ss) = g.Nb(1); % Blade count row 1
% %     N.N2(pp,ss) = g.Nb(2); % Blade count row 2
% %     N.Frs(pp,ss) = d.Fr; % Propulsive efficiency
% %     N.psis(pp,ss) = d.psi;
% %     N.FOMs(pp,ss) = d.Mf;
% %     Ltot = L.L1 + L.L2; 
% %     N.Ltot(pp,ss) = Ltot;
% %     N.L_eta(pp,ss) = L.eta_s;
% %     N.profL(pp,ss) = (L.Loss1.prof + L.Loss2.prof)/Ltot;
% %     N.baseL(pp,ss) = (L.Loss1.base + L.Loss2.base)/Ltot;
% %     N.tipL(pp,ss) = (L.Loss1.tip + L.Loss2.tip)/Ltot;
% %     N.endwallL(pp,ss) = (L.Loss1.endwall + L.Loss2.endwall)/Ltot;
% %     N.DH1(pp,ss) = d.DH(1);
% %     N.DH2(pp,ss) = d.DH(2);
% % end
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

% Finished loss loop so everything defined. The rest is gathering together
% outputs
% Hacking together arrays to plot

% if ~isnan(g.Nb(:))
%     C.phis(pp,ss) = d.phi;
%     C.sigmas(pp,ss) = d.sigma;
%     C.etas(pp,ss) = d.eta;  % Fan efficiency 
% %     C.Rs(pp,ss) = d.Reaction; % Reaction
%     C.rpms(pp,ss) = d.rpm1; % rpm
%     C.N1(pp,ss) = g.Nb(1); % Blade count row 1
%     C.N2(pp,ss) = g.Nb(2); % Blade count row 2
%     C.Frs(pp,ss) = d.Fr; % Propulsive efficiency
%     C.psis(pp,ss) = d.psi;
%     C.FOMs(pp,ss) = d.Mf;
%     Ltot = L.L1 + L.L2; 
%     C.Ltot(pp,ss) = Ltot;
%     C.L_eta(pp,ss) = L.eta_s;
%     C.profL(pp,ss) = (L.Loss1.prof + L.Loss2.prof)/Ltot;
%     C.baseL(pp,ss) = (L.Loss1.base + L.Loss2.base)/Ltot;
%     C.tipL(pp,ss) = (L.Loss1.tip + L.Loss2.tip)/Ltot;
%     C.endwallL(pp,ss) = (L.Loss1.endwall + L.Loss2.endwall)/Ltot;
%     C.DH1(pp,ss) = d.DH(1);
%     C.DH2(pp,ss) = d.DH(2);
% end

%% Store outputs in structures
% design(pp,ss).phisig = [d.phi, d.sigma];
% design(pp,ss).NRF = NRF; design(pp,ss).NRF.phisig = [d.phi, d.sigma];
% design(pp,ss).CRF = CRF; design(pp,ss).CRF.phisig = [d.phi, d.sigma];
% 
% mtip(pp,ss) = d.Mtip*2;

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
% display(['phi = ' num2str(d.phi) '  sigma = ' num2str(d.sigma) '   Time = ' num2str(toc,3) ]);
%     end
% end

%% Plotting
% Limits used for plotting, ensures effective comparisons
lims.Nb = [0 25]; % Sets limit of colormap for blade numbers
lims.L = [0 0.6]; % Sets loss limits between 0 and 1
lims.Ltot = [min(min(N.Ltot(:)), min(C.Ltot(:))) max(max(N.Ltot(:)), max(C.Ltot(:)))];

% Need to improve de Haller plotting
PlotCharts(N,lims);
% PlotCharts(C,lims);
% PlotVels(design(24),1)
% PlotVels(blah)
% velocity triangles, profiles (hub tip mean), chord

% Plotting dH spanwise
% figure; box on; hold on;
% for i=1
% y=linspace(0,1,6);
% x=design(1,i+1).NRF.d.DHall;
% plot(x,y)
% end
% legend;
% plot(0.72*ones(6),y,'r','LineWidth',2)
% xlabel('dH'); ylabel('% of span')
% 
% 
% figure;
% contourf(C.phis,C.sigmas,C.DH1); colorbar;
