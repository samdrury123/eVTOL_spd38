function [x] = PlotFlow(rc, rh, radius, sections, phi, psi, psi_ts, delta, R, S, carter, V, ang, omega, reaction, contourflag, meshv)
global crdfflag
%% Plot PSI/PHI Distirubutions

x = (radius-rh)./(rc-rh);

%% 3D Plot blades
figure(2); hold on; grid on; box on; axis equal; title('Blade Geometry');
mesh(squeeze(R.XYZ(:,1,:)),squeeze(R.XYZ(:,3,:)),squeeze(R.XYZ(:,2,:)));
mesh(squeeze(S.XYZ(:,1,:)),squeeze(S.XYZ(:,3,:)),squeeze(S.XYZ(:,2,:)));
plot3(S.centerline(:,1,:),S.centerline(:,3,:),S.centerline(:,2,:));

%% Velocity Triangles
figure(1); subplot(1,2,1); title('Blade Sections');
subplot(1,2,2); hold on; axis equal; 
VT_rad = sections(linspace(1,end,3));
dx = 4*V.sec.x(end);
dy = omega.*VT_rad(end);
XSHIFT = [dx dx dx; dx dx dx];
YSHIFT = [-dy 0 dy; -dy 0 dy];
title('Velocity Triangles'); xlim([-.2*dx, 1.5*dx]);
text(dx/2,dy,'CASING');text(dx/2,0,'MID');text(dx/2,-dy,'HUB');

% ROTOR INLET
X1.X(1:2,:)     = [zeros(1,3); V.sec.x(linspace(1,end,3))];
X1.Y(1:2,:)     = [zeros(1,3); zeros(1,3)] + YSHIFT;
plot(X1.X, X1.Y, 'k');

REL1.X(1:2,:)   = [zeros(1,3); V.sec.x(linspace(1,end,3))];
REL1.Y(1:2,:)   = [zeros(1,3); omega.*VT_rad] + YSHIFT;
plot(REL1.X, REL1.Y, 'k');

U1.X(1:2,:)   = [V.sec.x(linspace(1,end,3)); V.sec.x(linspace(1,end,3))];
U1.Y(1:2,:)   = [zeros(1,3); omega.*VT_rad] + YSHIFT;
plot(U1.X, U1.Y, 'k');

% STATOR INLET
X2.X(1:2,:)     = [zeros(1,3); V.sec.x(linspace(1,end,3))] + XSHIFT;
X2.Y(1:2,:)     = [zeros(1,3); zeros(1,3)] + YSHIFT;
plot(X2.X, X2.Y, 'b');

V2.X(1:2,:)   = [zeros(1,3); V.sec.x(linspace(1,end,3))] + XSHIFT;
V2.Y(1:2,:)   = [zeros(1,3); -V.sec.theta2(linspace(1,end,3))] + YSHIFT;
plot(V2.X, V2.Y, 'b');

REL2.X(1:2,:)   = [zeros(1,3); V.sec.x(linspace(1,end,3))] + XSHIFT;
REL2.Y(1:2,:)   = [zeros(1,3); omega.*VT_rad-V.sec.theta2(linspace(1,end,3))] + YSHIFT;
plot(REL2.X, REL2.Y, 'b');

U2.X(1:2,:)   = [V.sec.x(linspace(1,end,3)); V.sec.x(linspace(1,end,3))] + XSHIFT;
U2.Y(1:2,:)   = [-V.sec.theta2(linspace(1,end,3)); omega.*VT_rad-V.sec.theta2(linspace(1,end,3))] + YSHIFT;
plot(U2.X, U2.Y, 'b');

%% PSI/PHI Spanwise Distribution
figure(3); set(gcf,'Position', [0 320 1500 1200]);

subplot(2,6,1);
hold on; title('\phi (Flow Coefficient)'); ylabel('% of Span'); xlabel('\phi');
plot(phi.span,x, 'k');

subplot(2,6,2);
hold on; title('\psi (Total-Total Stage Loading)'); ylabel('% of Span'); xlabel('\psi');
plot(psi.span,x, 'k');

subplot(2,6,3);
hold on; title('\psi_{ts} (Total-Static Stage Loading)'); ylabel('% of Span'); xlabel('\psi_{ts}');
plot(psi_ts.span,x, 'k');

%% Deviation
subplot(2,6,7);
hold on; title('\delta (Deviation)'); ylabel('% of Span'); xlabel('\delta (DEVIATION)');
plot(delta.span.R,x,'b');plot(delta.span.S,x,'r');legend('Rotor','Stator');

%% Pitch-Chord
subplot(2,6,8);
hold on; title('s/l (Pitch-Chord)'); ylabel('% of Span'); xlabel('s/l');
plot(R.span.pitchchord,x, 'b');plot(S.span.pitchchord,x, 'r');legend('Rotor','Stator');

%% Carter coefficient 'm'
subplot(2,6,9);
hold on; title('m (Carter variable)'); ylabel('% of Span'); xlabel('m');
plot(carter.span.R,x,'b'); 

if crdfflag == 1
    plot(carter.span.S, x,'r'); legend('Rotor','Stator');
else
    plot([carter.span.S carter.span.S],[0 1],'r'); legend('Rotor','Stator');
end


%% Velocities
subplot(2,6,4);
hold on; title('V_x'); ylabel('% of Span'); xlabel('V / ms^{-1}');
plot(V.span.x,x, 'k');

subplot(2,6,5);
hold on; title('V_{\theta 2}'); ylabel('% of Span'); xlabel('V / ms^{-1}');
plot(V.span.theta2,x, 'k');
plot(radius.*omega,x, 'b');
legend('V_{\theta 2}', 'U');

subplot(2,6,6);
hold on; title('Reaction'); ylabel('% of Span'); xlabel('Reaction');
plot(reaction.span,x, 'k');

%% Chord
subplot(2,6,10);
hold on; title('Rotor Chord'); ylabel('% of Span'); xlabel('Chord / mm');
plot(R.span.chord * 1000,x, 'b'); plot(S.span.chord * 1000,x, 'r');

%% Flow angles
subplot(2,6,11);
hold on; title('Rotor exit angles'); ylabel('% of Span'); xlabel('Degrees');
plot(ang.span.b2,x, 'k'); plot(-ang.span.chi2,x, 'b'); legend('\beta_2','\chi_2');

%% Omega/Mf Contour Plots

switch contourflag
    case 1
        figure(7); hold on;
        contourf(meshv.phi, sqrt(meshv.phi.^2./(2*meshv.psi)), meshv.power, 'LineColor','none'); colorbar; %cb = colorbar; %cb.TickLabels = {'6,000','8,000','10,000','12,000','14,000','16,000','18,000'};
        xlabel('\phi'); ylabel('\psi'); set(gcf, 'Position', [0 0 300 300]);
%         plot(linspace(0.0, 1, 100), 0.5*linspace(0.0, 1, 100).^2,'k', 'LineWidth', 2); 
%         plot(linspace(0.0, 1, 100), 0.5*linspace(0.0, 1, 100).^2/1.44,'-.k', 'LineWidth', 2);
%         legend('RPM', 'Diffuser = 1', 'Diffuser length limit', 'Location', 'northwest');
        
        figure(8); hold on;
        contourf(meshv.phi, meshv.psi, (meshv.propmass.chassis./meshv.propmass.total), 'LineColor','none'); cb = colorbar; %cb.Ticks = [0.5 1 1.5 2 2.5 3 3.5];%cb.TickLabels = {'6,000','8,000','10,000','12,000','14,000','16,000','18,000'};
        xlabel('\phi'); ylabel('\psi'); set(gcf, 'Position', [0 0 300 300]);
        plot(linspace(0.0, 1, 100), 0.5*linspace(0.0, 1, 100).^2,'k', 'LineWidth', 2); 
        plot(linspace(0.0, 1, 100), 0.5*linspace(0.0, 1, 100).^2/1.44,'-.k', 'LineWidth', 2);
        plot(linspace(0.0, 1, 100), 2*linspace(0.0, 1, 100).^2,'-.r', 'LineWidth', 2); ylim([0 0.5]);
        legend('M_F', 'Diffuser = 1', 'Diffuser length limit', 'Propeller max', 'Location', 'northwest');
    case 0
end
%% Report Figures
[phi_free, psi_free, ~, ~] = Distributions('free', phi, psi, radius,sqrt(0.5*(rc^2 + rh^2)), 2, radius);
[phi_cang, psi_cang, ~, ~] = Distributions('constangle', phi, psi, radius,sqrt(0.5*(rc^2 + rh^2)), 1.089, radius);
[phi_mixd, psi_mixd, ~, ~] = Distributions('custom', phi, psi, radius,sqrt(0.5*(rc^2 + rh^2)), 1.7, radius);

figure(4);

subplot(2,4,1);grid on; box on;
hold on; ylabel('% of Span', 'FontSize', 12); xlabel('\phi', 'FontSize', 12);
plot(phi_free,x, 'k');
plot(phi_cang,x, 'b');
plot(phi_mixd,x, 'r');
% legend('Free','Const. Angle','Mixed');

subplot(2,4,2);grid on; box on;
hold on; xlabel('\psi', 'FontSize', 12);yticklabels([' ',' ']);
plot(psi_free,x, 'k');
plot(psi_cang,x, 'b');
plot(psi_mixd,x, 'r');
% legend('Free','Const. Angle','Mixed');

subplot(2,4,3);grid on; box on;
hold on; xlabel('$V_x/(V_x)_m$', 'Interpreter', 'Latex', 'FontSize', 12);yticklabels([' ',' ']);
plot(phi_free.*omega.*radius/V.m.x,x, 'k');
plot(phi_cang.*omega.*radius/V.m.x,x, 'b');
plot(phi_mixd.*omega.*radius/V.m.x,x, 'r');
% legend('Free','Const. Angle','Mixed');

subplot(2,4,4);grid on; box on;
hold on; xlabel('$\alpha_2$', 'Interpreter', 'Latex', 'FontSize', 12);yticklabels([' ',' ']);
plot(atand(psi_free.*omega.*radius./(phi_free.*omega.*radius)),x, 'k');
plot(atand(psi_cang.*omega.*radius./(phi_cang.*omega.*radius)),x, 'b');
plot(atand(psi_mixd.*omega.*radius./(phi_mixd.*omega.*radius)),x, 'r');
legend('Free (A = -1, B = -2)','Const. Angle (A = B = -1.089)','Mixed (A = -1.029, B = -1.7)', 'Position', [1 0.5 0.15 0.0869]);

set(gcf, 'Position', [0 0 600 500]);
% subplot(1,4,4);
% hold on; title('\alpha_2 = tan^{-1}(V_{\theta 2} / V_x)'); xlabel('\alpha_2 / Degrees', 'FontSize', 12);
% plot(sqrt((psi_free.*omega.*radius).^2 + (phi_free.*omega.*radius).^2),x, 'k');
% plot(sqrt((psi_cang.*omega.*radius).^2 + (phi_cang.*omega.*radius).^2),x, 'b');
% plot(sqrt((psi_mixd.*omega.*radius).^2 + (phi_mixd.*omega.*radius).^2),x, 'r');
% legend('Free','Const. Angle','Mixed');

% figure(5);
% 
% subplot(1,2,1);
% hold on; title('Chord'); ylabel('% of Span'); xlabel('mm');
% plot(1000*R.span.chord,x, 'k');
% legend('Free','Mixed');
% 
% subplot(1,2,2);
% hold on; title('Pitch-Chord'); ylabel('% of Span'); xlabel('mm');
% plot(R.span.pitchchord,x, 'k');
% legend('Free','Mixed');

% save('ReportFigs','R','x');
%% CLEAR
% clear x;
end