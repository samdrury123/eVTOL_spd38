%% ---- FUNCTION: PLOTTING ---- %%
function PlotCharts(D,lim)

% Function plots charts of phi vs sigma and phi vs psi for design
% properties of interest
% Inputs:
% - D: design to be plotted (Normal- or Contra-rotating)
% - D.etas: fan efficiency
% - D.Frs:  propulsive efficiency (overall = fan*propulsive)
% - D.FOMs: figure of merits
% - D.rpms
% - D.N1:   number of rotor blades
% - D.N2:   number of stator blades
% - D.psis: flow coefficient
% - D.Ltot:     total loss
% - D.L_eta:    efficiency from dh/(dh+Tds)
% - D.profL:    fractional profile loss
% - D.baseL:    fractional base loss
% - D.tipL:     fractional tip loss
% - D.endwallL: fractional endwall loss

%% 
phis=D.phis; sigmas=D.sigmas; psis=D.psis;

figure;set(gcf, 'color', 'w'); grid off; box on;
set(gcf,'Position',[20 50 1300 700]);
subplot(2,4,1)
contourf(phis,sigmas,D.etas,35, 'edgecolor','none'); colorbar; caxis(caxis); hold on;
% [C,h]=contour(phis,sigmas,D.DH1,6, 'r:', 'LineWidth',1.1, 'ShowText','on');
[C1,h1]=contour(phis,sigmas,D.DH1, 'r:', 'LineWidth',1.1, 'LabelSpacing',300);
clabel(C1,h1); h1.LevelList=0.68:0.04:max(D.DH1(:)); %h.LevelStep=0.02;
% contour(phis,sigmas,D.DH1, [0.72 0.72], 'r--', 'LineWidth',2);
title("Fan efficiency \eta_a with de Haller number (>0.72)"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,2)
contourf(phis,sigmas,D.Frs,25, 'edgecolor','none'); colorbar
title("Propulsive efficiency \eta_p"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,3)
contourf(phis,sigmas,D.etas.*D.Frs,25, 'edgecolor','none'); colorbar
title("Overall efficiency \eta_{ov}"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,4)
contourf(phis,sigmas,D.FOMs,25, 'edgecolor','none'); colorbar
title("Figure of Merit"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,5)
contourf(phis,sigmas,D.rpms,25, 'edgecolor','none'); colorbar
title("RPM"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,6)
contourf(phis,sigmas,D.N1,25, 'edgecolor','none'); colorbar; caxis(lim.Nb);
title("Rotor blades"); xlabel("\phi"); ylabel("\sigma");
% text(1.1,1.3, ['CRDF = ' num2str(CRDF_flag)],'color','b')
subplot(2,4,7)
contourf(phis,sigmas,D.N2,25, 'edgecolor','none'); colorbar; caxis(lim.Nb);
title("Stator blades"); xlabel("\phi"); ylabel("\sigma");
subplot(2,4,8)
contourf(phis,sigmas,D.psis,25, 'edgecolor','none'); colorbar
title("Psi \psi"); xlabel("\phi"); ylabel("\sigma");

figure;set(gcf, 'color', 'w'); grid off; box on;
set(gcf,'Position',[20 50 1300 700]);
subplot(2,3,1)
contourf(phis,sigmas,D.profL,25, 'edgecolor','none'); colorbar; caxis(lim.L);
title("Fractional profile loss"); xlabel("\phi"); ylabel("\sigma");
subplot(2,3,2)
contourf(phis,sigmas,D.baseL,25, 'edgecolor','none'); colorbar; caxis(lim.L);
title("Fractional base pressure loss"); xlabel("\phi"); ylabel("\sigma");
subplot(2,3,4)
contourf(phis,sigmas,D.tipL,25, 'edgecolor','none'); colorbar; caxis(lim.L);
title("Fractional tip loss"); xlabel("\phi"); ylabel("\sigma");
subplot(2,3,5)
contourf(phis,sigmas,D.endwallL,25, 'edgecolor','none'); colorbar; caxis(lim.L);
title("Fractional endwall loss"); xlabel("\phi"); ylabel("\sigma");
subplot(2,3,3)
contourf(phis,sigmas,D.Ltot,25, 'edgecolor','none'); colorbar; caxis(lim.Ltot);
if strcmp(D.name,'NRF')
    title("NRF Total Loss");
else
    title("CRF Total Loss");
end
xlabel("\phi"); ylabel("\sigma");
subplot(2,3,6)
contourf(phis,sigmas,D.L_eta,25, 'edgecolor','none'); colorbar
title("Entropy efficiency"); xlabel("\phi"); ylabel("\sigma");

%Sigmas
figure;set(gcf, 'color', 'w'); grid off; box on;
set(gcf,'Position',[200 50 1200 600]);
subplot(1,2,1)
contourf(phis,sigmas,D.DH1,25, 'edgecolor','none'); colorbar; caxis(caxis);
hold on;
contour(phis,sigmas,D.DH1, [0.72 0.72], 'r--', 'LineWidth',2);
[C2,h2] = contour(phis,sigmas,D.L_eta, 'r:', 'LineWidth',1.5);
clabel(C2,h2); h2.LabelSpacing = 500; 
title("Blade 1 de Haller number ( > 0.72 for minimal loss) with contours of efficiency"); xlabel("\phi"); ylabel("\sigma"); hold off;
subplot(1,2,2)
contourf(phis,sigmas,D.DH2,25, 'edgecolor','none'); colorbar; caxis(caxis);
hold on;
contour(phis,sigmas,D.DH2, [0.72 0.72], 'r--', 'LineWidth',2);
[C3,h3] = contour(phis,sigmas,D.L_eta, 'r:', 'LineWidth',1.5);
clabel(C3,h3); h3.LabelSpacing = 500;
title("Blade 2 de Haller number ( > 0.72 ) with contours of efficiency"); xlabel("\phi"); ylabel("\sigma"); hold off;

% Psis
% figure;set(gcf, 'color', 'w'); grid off; box on;
% set(gcf,'Position',[200 50 1200 600]);
% subplot(1,2,1)
% contourf(phis,psis,D.DH1,25, 'edgecolor','none'); colorbar; caxis(caxis)
% hold on;
% contour(phis,psis,D.DH1, [0.72 0.72], 'w--', 'LineWidth',2);
% contour(phis,psis,D.L_eta,25, 'r:', 'LineWidth',1.5);
% title("Blade 1 de Haller number ( > 0.72 for minimal loss) with contours of efficiency"); xlabel("\phi"); ylabel("\psi"); hold off;
% subplot(1,2,2)
% contourf(phis,psis,D.DH2,25, 'edgecolor','none'); colorbar; caxis(caxis)
% hold on;
% contour(phis,psis,D.DH2, [0.72 0.72], 'w--', 'LineWidth',2);
% contour(phis,psis,D.L_eta,25, 'r:', 'LineWidth',1.5);
% title("Blade 2 de Haller number ( > 0.72 ) with contours of efficiency"); xlabel("\phi"); ylabel("\psi"); hold off;



figure;set(gcf, 'color', 'w'); grid off; box on;
set(gcf,'Position',[20 50 1000 600]);
subplot(2,3,1)
contourf(phis,psis,D.etas,25, 'edgecolor','none'); colorbar
title("Fan efficiency \eta_a"); xlabel("\phi"); ylabel("\psi");
subplot(2,3,2)
contourf(phis,psis,D.Frs,25, 'edgecolor','none'); colorbar
title("Propulsive efficiency \eta_p"); xlabel("\phi"); ylabel("\psi");
if max(D.Frs(:)) == 0
    subplot(2,3,3)
    contourf(phis,psis,D.FOMs,25, 'edgecolor','none'); colorbar
    title("Figure of Merit"); xlabel("\phi"); ylabel("\psi");
else
    subplot(2,3,3)
    contourf(phis,psis,D.etas.*D.Frs,25, 'edgecolor','none'); colorbar
    title("Overall efficiency \eta_{ov}"); xlabel("\phi"); ylabel("\psi");
end
subplot(2,3,4)
contourf(phis,psis,D.rpms,25, 'edgecolor','none'); colorbar
title("RPM"); xlabel("\phi"); ylabel("\psi");
% text(1.1,0.8, ['CRDF = ' num2str(CRDF_flag)],'color','b')
subplot(2,3,5)
contourf(phis,psis,D.N1,25, 'edgecolor','none'); colorbar; caxis(lim.Nb);
title("Rotor blades"); xlabel("\phi"); ylabel("\psi");
subplot(2,3,6)
contourf(phis,psis,D.N2,25, 'edgecolor','none'); colorbar; caxis(lim.Nb);
title("Stator blades"); xlabel("\phi"); ylabel("\psi");
