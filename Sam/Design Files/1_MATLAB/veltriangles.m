%% ---- FUNCTION: VELOCITY TRIANGLES ---- %%
function [a,g,d] = veltriangles(d,g,k,g2,dev,p,fignum)
% Function calculates velocities and angles for conventional or
% counter-rotating fan. 
% Currently works with alpha1=alpha3=0, single values of phi, psitt
% Inputs:
% - d.phi: flow coefficient 
% - d.psitt: stage loading coefficient
% - d.Umid: Meanline blade speed
% - d.rpm: shaft speed in rpm
% - k: rpm split between blade A and blade B, k=1 is conventional
% - g: geometry of fan
% - p: if p==1 plot velocity triangles
% Outputs:
% - a: Structure with aerodynamic velocities and angles, plus relative reaction
%      3d arrays of velocities and angles have dimensions: 
%      - radius
%      - position (upB1, dwnB1, dwnB2)
%      - reference frame (absolute, B1, B2)
% - g: geometry updated with spanwise chord, pitch, blade metal angles, etc
%      chi are arrays with blade metal angles and have dimensions:
%      - radius
%      - 1: LE, 2: TE
%      - 1: Blade 1 2: Blade 2
% - d: updated with blade rpms

% Set up radial positions
Nr = 31;
r = linspace(g.rhub,g.rcas,Nr)';
g.r=r;

% Split midspan blade velocities between blade blade rows
Uf = r ./ g.rmid * d.Umid;
Ua = k*Uf;
Ub = -(1-k)*Uf;
d.rpm1 = d.rpm*k; d.rpm2 = -d.rpm*(1-k);

% Set up spanwise stage loading coefficient, enthalpy and vtheta
psitt = d.psitt .* (g.rmid ./ r).^ (g2.n+1);
imid = 20 ; % THIS IS A HACK TO GET A POSITION CLOSE TO MIDLINE RATHER THAN MIDSPAN (Nr + 1)/2;
g.n = g2.n;

% Calculate Vx from midline relative blade speed
vx(1:Nr,1) = ones(size(r)) .* d.phi .* d.Umid;

% Calculate velocities and angles upstream of blade 1, n.b. alpha1 = 0
vth(1:Nr,1,1) = 0;
vth(:,1,2) = vth(:,1,1) - Ua;
vth(:,1,3) = vth(:,1,1) - Ub;

% Calculate velocities and angles downstream of blade 1, including axial
% velocity distribution from radial equilibrium equation
vth(:,2,1) = Ua .* psitt / k;
vth(:,2,2) = vth(:,2,1) - Ua;
vth(:,2,3) = vth(:,2,1) - Ub;

% Use Euler work equation to calculate change in enthalpy
dh0 = Ua .* vth(:,2,1);

% Numerically solve Simple Radial Equilibrium equation
% Boundary Condition
vx(imid,2)  = d.Umid * d.phi;

% Radial increment
dr = r(2)-r(1);

% Solve
for ii = imid:(Nr-1)  
   vx(ii+1,2)     = vx(ii,2) + (dh0(ii+1) - dh0(ii))/vx(ii,2) - (vth(ii,2,1)/vx(ii,2)) * (vth(ii+1,2,1) - vth(ii,2,1)) - ((vth(ii,2,1)).*(vth(ii+1,2,1)) / vx(ii,2)) * (dr/r(ii));
end
 
for ii = imid:-1:2   
   vx(ii-1,2)     = -(-vx(ii,2) + (dh0(ii) - dh0(ii-1))/vx(ii,2) - (vth(ii,2,1)/vx(ii,2)) * (vth(ii,2,1) - vth(ii-1,2,1)) - ((vth(ii,2,1)).*(vth(ii-1,2,1)) / vx(ii,2)) * (dr/r(ii)));
end

% Calculate velocities and angles downstream of blade 2, n.b. alpha3 = 0
vth(:,3,1) = 0;
vth(:,3,2) = vth(:,3,1) - Ua;
vth(:,3,3) = vth(:,3,1) - Ub;

vx(:,3,:)=vx(:,1,:);

% Calculate angles and velocities at three stations
for ss = 1:3
   for tt = 1:3
      vx(:,ss,tt) = vx(:,ss,1);
      alpha(:,ss,tt) = atand(vth(:,ss,tt) ./ vx(:,ss,tt));
      v(:,ss,tt) = vx(:,ss,tt) ./ cosd(alpha(:,ss,tt));
   end
end

% Calculate relative Reaction
R = -(vth(:,1,2)+vth(:,2,2)) ./ (2*Uf);

for bb=1:2

    
   % --- Work Out Blade Count --- %
   vup = v(:,bb,bb+1);
   vdn = v(:,bb+1,bb+1);
   vupth = vth(:,bb,bb+1);
   vdnth = vth(:,bb+1,bb+1);
   scDF(1:Nr,bb) = (g2.DF - (1 - vdn./vup)) * 2 .* vup  ./ abs((vupth - vdnth));
   g.DF = g2.DF;
   
   if g2.blades(bb) ~= 0
      % User input for blade count
      g.N(bb) = g2.blades(bb);   
   else
   
   % Calculate blade count based on midspan pitch-to-chord
   chordmid = (g.rcas-g.rhub)/g2.AR;
   pitchmid(bb) = scDF(round(Nr/2),bb) * chordmid;
   N(bb) = 2*pi*g.rmid ./ pitchmid(bb);
   g.N(bb) = round(N(bb));
   g.AR = g2.AR;
   end 
   % Calculate pitch with actual blade count
   g.s(1:Nr,bb) = r .* 2*pi ./ g.N(bb);
   
   % --- Calculate chords ---
   % Multiplier to get angles in right direction while doing deviation calcs
   bl = round(bb-1.5);
   
   % Get angles into "cascade" frame, i.e. all positive, for deviation calcs
   aup = bl * alpha(:,bb,bb+1);
   adwn =  bl * alpha(:,bb+1,bb+1);
   
   if dev==0 
       % Calculate chord and deviation to give constant diffusion factor across span
       % Calculate chord based on Lieblein
       g.c(1:Nr,bb) = g.s(:,bb) ./ scDF(:,bb);

       % Calculate blade metal angles with deviation using Howell and Carter
       % Blade LE metal angle, zero incidence 
       chi(:,1,bb) = aup;
       g.chi(:,1,bb) = bl .* chi(:,1,bb);
       
       % Howell and Carter
       ac = 0.5; % Maximum camber position to chord ratio
       m = 0.23*(2*ac).^2 + adwn / 500;
       chi(:,2,bb) = (adwn - m .* aup .* (g.s(:,bb) ./ g.c(:,bb)).^ 0.5) ./ (1 - m.*(g.s(:,bb) ./ g.c(:,bb)).^0.5);
       g.chi(:,2,bb) = bl .* chi(:,2,bb);
       devi(:,bb) = chi(:,2,bb) - adwn; 
       
   else
       % Calculate chord to given constant deviation across span based on midline value
       % Blade LE metal angle, zero incidence
       chi(:,1,bb) = aup;
       g.chi(:,1,bb) = bl .* chi(:,1,bb);
       
       % Calculate midspan chord based on Lieblein
       cmid = g.s(imid,bb) ./ scDF(imid,bb) ;
       
       % Calculate midspan deviation
       ac=0.5;
       m = 0.23*(2*ac).^2 + adwn(imid) / 500;
       chimid = (adwn(imid) - m .* aup(imid) .* (g.s(imid,bb) ./ cmid).^ 0.5) ./ (1 - m.*(g.s(imid,bb) ./ cmid).^0.5);
       devi(bb) =  (adwn(imid) - chimid);
       
       % Calculate chord based on constant deviation
       m = 0.23*(2*ac).^2 + adwn / 500;
       g.c(1:Nr,bb) = g.s(:,bb) .* (devi(bb) ./ (m .* (aup + devi(bb) - adwn) )).^-2;
       chi(:,2,bb) = adwn - devi(bb);
       g.chi(:,2,bb) = bl .* chi(:,2,bb);
       
   end
   
   if isfield(g2,'cman')==1
      % Adjust chord with manual input
      g.c(:,bb) = g.c(:,bb) .* linspace(g2.cman(bb,1),g2.cman(bb,2), size(g.c,1))';
   end
   % Recalculate deviation and blade metal angles
   % Howell and Carter
   ac = 0.5; % Maximum camber position to chord ratio
   m = 0.23*(2*ac).^2 + adwn / 500;
   chi(:,2,bb) = (adwn - m .* aup .* (g.s(:,bb) ./ g.c(:,bb)).^ 0.5) ./ (1 - m.*(g.s(:,bb) ./ g.c(:,bb)).^0.5);
   g.chi(:,2,bb) = bl .* chi(:,2,bb);
      
   % Recalculate Diffusion Factor based on chords calculated 
   DF(1:Nr,bb) = (1 - vdn./vup) + abs(vupth - vdnth) ./ (2.*vup) .* (g.s(:,bb) ./ g.c(:,bb));
   
end

if p==1
   % Plot Velocity Triangles 
   figure; set(gcf,'Position',[20 40 1200 600]);
   set(gcf, 'color', 'w'); 
   
   cols = {'-ko','-ro','-bo' ; '-kx','-rx','-bx' ; '-ksq','-rsq','-bsq'};
   rind = [1,round(Nr/2),Nr];
   for rr = 1:3
      for ss = 1:3
         for ff = 1:3
            
            stx = ss*100;
            sty = 0;
            endx = stx+v(rind(rr),ss,ff)*cosd(alpha(rind(rr),ss,ff));
            endy = sty+v(rind(rr),ss,ff)*sind(alpha(rind(rr),ss,ff));
            
            subplot(4,6,[1:2,7:8]); hold on; box on
            plot([sty endy],[-stx -endx],cols{rr,ff})
            xlabel('Velocity (m/s)')
            set(gca,'ytick',[])
            xlim([-120 120])
         end   
      end
   end
        
   % Plot spanwise distributions
   span = (r - g.rhub)./(g.rcas-g.rhub);
   
   % Vx
   subplot(4,6,[14 20]); hold on;grid off; box on;
   plot(vx(:,1,1),span)
   plot(vx(:,2,1),span)
   plot(vx(:,3,1),span)
   xlim([25 35]);
   ylabel('Span'); xlabel('Vx (m/s)');
   legend('Up B1', 'Dwn B1','DwnB2','Location','SouthOutside')
  
   % Reaction
   subplot(4,6,[13 19]); hold on;grid off; box on;
   plot(R,span)
   xlabel('Relative Reaction');
   xlim([0 1])
   
   % Blade 1 relative flow and metal angles
   subplot(4,6,[15 21]); hold on;grid off; box on;
   plot(alpha(:,1,2),span)
   plot(alpha(:,2,2),span)
   plot(g.chi(:,2,1),span)
   xlim([-70 20])
   xlabel('Rel Blade 1 Flow Angle (deg)')
   legend('Up B1', 'Dwn B1','DwnB1 Metal','Location','SouthOutside')
      
   % Blade 2 relative flow and metal angles
   subplot(4,6,[16 22]); hold on; grid off; box on;
   plot(alpha(:,2,3),span)
   plot(alpha(:,3,3),span)
   plot(g.chi(:,2,2),span)
   xlim([-20 70])
   xlabel('Rel Blade 2 Flow Angle (deg)')
   legend('Up B2', 'Dwn B2','DwnB2 Metal','Location','SouthOutside')

   % Diffusion factor
   subplot(4,6,[17 23]); hold on;grid off; box on;
   plot(DF,span);
   plot([0.45 0.45],[0 1],'--r')
   xlabel('Diffusion Factor')
   legend('B1', 'B2','Location','SouthOutside')
   xlim([0 0.6])
   
   % Chord
   subplot(4,6,[18 24]);  hold on;grid off; box on;
   plot(g.c,span);
   xlabel('Chord (m)')
   legend('B1', 'B2','Location','SouthOutside')
   xlim([0.005 0.025])
   % Display blade counts
   text(0.17,0.2,['B1: ' num2str(g.N(1))])
   text(0.17,0.1,['B2: ' num2str(g.N(2))])
   
end

%%SAM CREATING NEW GRAPH FOR TMR
if false
    
    % Plot spanwise distributions
    span = (r - g.rhub)./(g.rcas-g.rhub);

    % Vx
    figure(6); set(gcf,'Position',[20 40 200 450]);
    hold on;grid off; box on;
    plot(vx(:,1,1),span)
    plot(vx(:,2,1),span)
    plot(vx(:,3,1),span)
    xlim([25 35]);
    ylabel('Span'); xlabel('Vx (m/s)');
    legend('Up B1', 'Dwn B1','DwnB2','Location','SouthOutside')

    % Blade 1 relative flow and metal angles
    figure(7); set(gcf,'Position',[20 40 200 450]);
    hold on;grid off; box on;
    plot(alpha(:,1,2),span)
    plot(alpha(:,2,2),span)
    plot(g.chi(:,2,1),span)
    xlim([-70 20])
    xlabel('Rel Blade 1 Flow Angle (deg)')
    legend('Up B1', 'Dwn B1','DwnB1 Metal','Location','SouthOutside')

    % Blade 2 relative flow and metal angles
    figure(8); set(gcf,'Position',[20 40 200 450]);
    hold on; grid off; box on;
    plot(alpha(:,2,3),span)
    plot(alpha(:,3,3),span)
    plot(g.chi(:,2,2),span)
    xlim([-20 70])
    xlabel('Rel Blade 2 Flow Angle (deg)')
    legend('Up B2', 'Dwn B2','DwnB2 Metal','Location','SouthOutside')

    % Diffusion factor
    figure(9); set(gcf,'Position',[20 40 200 450]);
    hold on;grid off; box on;
    plot(DF,span);
    plot([0.45 0.45],[0 1],'--r')
    xlabel('Diffusion Factor')
    legend('B1', 'B2','Location','SouthOutside')
    xlim([0 0.6])

    % Chord
    figure(10); set(gcf,'Position',[20 40 200 450]);
    hold on;grid off; box on;
    plot(g.c,span);
    xlabel('Chord (m)')
    legend('B1', 'B2','Location','SouthOutside')
    xlim([0.005 0.025])
    % Display blade counts
    text(0.17,0.2,['B1: ' num2str(g.N(1))])
    text(0.17,0.1,['B2: ' num2str(g.N(2))])
end

% Gather together aerodynamic outputs into structure
a.v = v; a.vth = vth; a.vx = vx;
a.alpha = alpha; a.dev = dev;
a.R = R;

end