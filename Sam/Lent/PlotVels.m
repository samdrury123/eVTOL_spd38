%% ---- FUNCTION: PLOTTING ---- %%
function PlotVels(D,n)

% Function plots, for a given design:
% - velocity profiles at hub, mean line and tip
% - velocity triangles at hub, mean line and tip
% - spanwise chord distribution
% - spanwise distribution factor

% Inputs:
% - D: design to be plotted, D=design(i) (Normal- or Contra-rotating)
% - n: 1 = NRF, 2 = CRF

%% Unpack variables
if n == 1
    a = D.NRF.a;
    g = D.NRF.g;
    d = D.NRF.d;
elseif n == 2
    a = D.CRF.a;
    g = D.CRF.g;
    d = D.CRF.d;
else 
    disp('Incorrect indice entered - 1=NRF, 2=CRF');
    return
end

phi = D.phisig(1);
sigma = D.phisig(2);
psi = d.psi;
DF = a.DFspan;
c = g.c;
Vss = a.Vss;
Vps = a.Vps;
alpha = a.alpha;
v = a.v;
r = g.r;
span = (r - r(1)) / (r(end) - r(1));
imid = a.imid;
x = linspace(0,1,size(Vss,2));

%% Plotting

figure;set(gcf, 'color', 'w'); grid off; box on;
set(gcf,'Position',[20 150 1500 600]);

% Velocity profile
Nr = size(r,1);
rind = [1,imid,Nr];
% Vss(rind,:,1)'
% DF
cols = {'-o';'-x';'-sq'};
% legend('B1', 'B2','Location','SouthOutside')

% Concatonate arrays
V1 = [Vss(rind,:,1) Vps(rind,:,1)]';
V2 = [Vss(rind,:,2) Vps(rind,:,2)]';
xs = [x x]';
lim = [0 round(max(max(V1(:)), max(V2(:)))+50,-2)];
for ii = 1:3 % radial position
    subplot(2,6,1:2); hold on; box on
    plot(xs,V1(:,ii),cols{ii});
    xlabel('x/{c_x}')
    ylabel('Velocity (m/s)')
    ylim(lim)
    xlim([-0.1 1.1])
    
    %     set(gca,'xtick',[])
    subplot(2,6,7:8); hold on; box on
    plot(xs,V2(:,ii),cols{ii});
    xlabel('x/{c_x}')
    ylabel('Velocity (m/s)')
    ylim(lim)
    xlim([-0.1 1.1])
    %     set(gca,'ytick',[])
end

subplot(2,6,1:2); hold on; grid on
set(gca,'xtick',[0 1], 'ygrid', 'off', 'xMinorTick', 'on', 'xMinorGrid','on')
legend('Hub', 'Mean Line', 'Tip')
title("Velocity profiles for Blade 1")

subplot(2,6,7:8); hold on; grid on
set(gca,'xtick',[0 1], 'ygrid', 'off', 'xMinorTick', 'on', 'xMinorGrid','on')
legend('Hub', 'Mean Line', 'Tip')
title("Velocity profiles for Blade 2")

% Diffusion factor
subplot(2,6,[5,11]); hold on; grid off; box on;
plot(DF,span);
plot([0.45 0.45],[0 1],'--r')
xlabel('Diffusion Factor')
legend('B1', 'B2','Location','SouthOutside')
xlim([0 0.6])

% Chord
subplot(2,6,[6,12]);  hold on; grid off; box on;
plot(c,span);
xlabel('Chord (m)')
legend('B1', 'B2','Location','SouthOutside')
xlim([0 1.2*max(c(:))])

% Velocity triangles
cols2 = {'-ko','-ro','-bo' ; '-kx','-rx','-bx' ; '-ksq','-rsq','-bsq'};
for rr = 1:3 % radial position
  for ss = 1:3 % axial position
     for ff = 1:3 % reference frame

        stx = ss*100;
        sty = 0;
        endx = stx+v(rind(rr),ss,ff)*cosd(alpha(rind(rr),ss,ff));
        endy = sty+v(rind(rr),ss,ff)*sind(alpha(rind(rr),ss,ff));

        subplot(2,6,[3:4,9:10]); hold on; box on
        plot([sty endy],[-stx -endx],cols2{rr,ff})
        xlabel('Velocity (m/s)')
        set(gca,'ytick',[])
        xlim([-100 100])
     end   
  end
end
title(['Phi: ' num2str(phi) '  Sigma: ' num2str(sigma) '  Psi: ' num2str(psi,3)])


%% Testing area

% Blade flow/metal angles
% figure;
% % Blade 1 relative flow and metal angles
%    subplot(4,6,[15 21]); hold on;grid off; box on;
%    plot(alpha(:,1,2),span)
%    plot(alpha(:,2,2),span)
%    plot(g.chi(:,2,1),span)
%    xlim([-70 20])
%    xlabel('Rel Blade 1 Flow Angle (deg)')
%    legend('Up B1', 'Dwn B1','DwnB1 Metal','Location','SouthOutside')
%       
%    % Blade 2 relative flow and metal angles
%    subplot(4,6,[16 22]); hold on; grid off; box on;
%    plot(alpha(:,2,3),span)
%    plot(alpha(:,3,3),span)
%    plot(g.chi(:,2,2),span)
%    xlim([-20 70])
%    xlabel('Rel Blade 2 Flow Angle (deg)')
%    legend('Up B2', 'Dwn B2','DwnB2 Metal','Location','SouthOutside')

