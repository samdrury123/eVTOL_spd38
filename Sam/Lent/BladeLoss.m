%% ---- FUNCTION: LOSS EVALUATION ---- %%
function [g,L] = BladeLoss(d,g,a,q,L,k)

%%REDOO intro
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

%% Should unpack some stuff - gets messy

% Set up angles and speeds for each row
% ang1 = [a.a1r a.a2r(2)];
% ang2 = [a.a2r(1) a.a3r];
% vel1 = [a.v1r a.v2r(2)];
% vel2 = [a.v2r(1) a.v3r];
% vel1t = [a.vt1r a.vt2r(2)];
% vel2t = [a.vt2r(1) a.vt3r];
% vel1x = [a.vx1 a.vx2];
imid=4; %%CHANGE THISSS
ang1 = [a.alpha(imid,1,2) a.alpha(imid,2,3)];
ang2 = [a.alpha(imid,2,2) a.alpha(imid,3,3)];
vel1 = [a.v(imid,1,2) a.v(imid,2,3)];
vel2 = [a.v(imid,2,2) a.v(imid,3,3)];
vel1t = [a.vth(imid,1,2) a.vth(imid,2,3)];
vel2t = [a.vth(imid,2,2) a.vth(imid,3,3)];
vel1x = [a.vx(imid,1,1) a.vx(imid,2,1)]; 
%vel2x = [a.vx2 a.vx3];
Temp1 = [q.T1 q.T2];
rho = [q.ro1 q.ro2];

Nb =        [0,0];
BLT_SS =    [0,0];
BLT_PS =    [0,0];
prof1 =     [0,0];
prof2 =     [0,0];
base =      [0,0];
tip =       [0,0];
endwall =   [0,0];

% Loop through each blade row
for rr=1:2
  % Blade Geometry
  C = (g.rc - g.rh)/g.AR; % chord %%THIS WILL NEED CHANGED FOR RADIAL EQ
  stag = (ang1(rr) + ang2(rr))/2; % stagger angle
  Cx = C * cosd(stag); % axial chord

  % Use Lieblein to get pitch-to-chord
  delvt = abs(vel1t(rr) - vel2t(rr));
  
  sc_rat = (d.DF - (1 - vel2(rr)./vel1(rr))) * 2 .* vel1(rr)  ./ delvt;
  s = C*sc_rat;
  Nb(rr) = round(2*pi*g.rm / s);
  if Nb(rr)<1 || Nb(rr)>100
      Nb(rr) = NaN;
  end
  s = 2 * pi * g.rm / Nb(rr);
  solid = C/s;

  % Shape factors for BL equations
  H = 2;            % DT / MT - **How was this selected?**
  He = 4*H/(3*H-1); % ET / MT

  %Dissipation Coefficient as function of Reynolds
  Re = vel1(rr)*C*rho(rr) / 18e-6;
  Cd = d.Cdnom*(Re/500000)^-0.2; % Correction due to effect of low Re (Dickens p154)

  % First guess in boundary layer displacement thickness (BLDT)
  BLDT_SS = 0.0001; BLDT_SS_g = 0;
  BLDT_PS = 0.0001; BLDT_PS_g = 0;

  % Iterate velocity profiles and boundary layers. This is explained in Dickens thesis
  while abs((BLDT_SS - BLDT_SS_g)/BLDT_SS) > 1e-10 && abs((BLDT_PS - BLDT_PS_g)/BLDT_PS) > 1e-10
    BLDT_SS_g = BLDT_SS;
    BLDT_PS_g = BLDT_PS;
   
    vTE = vel2(rr) * (1 - solid * ( (g.tTE + (BLDT_SS + BLDT_PS)/C) / cosd(ang2(rr))) )^-1; % Dickens 9.2
    vLE = vel1(rr) * (0.98 + 0.5 * solid * g.tmax * vel1x(rr) / delvt); % Dickens 9.3
   
    ATD = 0; % Roof-top parameter set to 0 - Dickens p147/8
    delv = 3 * vel1x(rr) * delvt / (solid * cosd(stag) * ((2*vLE*(1+2*ATD))+vTE*(1-ATD))); % Dickens 9.1 
    vLE_SS = vLE + delv;
    vLE_PS = vLE - delv;
   
    profilefun_SS1 =  @(z) (z*0 + vLE_SS).^3;
    profilefun_SS2 = @(z) (((vLE_SS - (ATD * vTE)) / (1-ATD)) - (z/C*(vLE_SS - vTE)/(1-ATD))).^3; % Dickens D.5
    BLET_SS = (2/vTE^3) * Cd * integral(profilefun_SS1,0,ATD*C) + (2/vTE^3) * Cd * integral(profilefun_SS2,ATD*C,C); % Dickens 9.10
    % Check Dickens 9.12 to improve this, if using the density ratio - not
    % an issue for low M, but 8% variation possibly
    BLDT_SS = BLET_SS/(3*He - 4); % Dickens 9.8
    BLMT_SS = BLET_SS/He;
   
    profilefun_PS1 =  @(z) (z*0 + vLE_PS).^3;
    profilefun_PS2 = @(z) (((vLE_PS - (ATD * vTE))/(1-ATD)) - (z/C*(vLE_PS - vTE)/(1-ATD))).^3; % Dickens D.6
    BLET_PS= (2/vTE^3) * Cd * integral(profilefun_PS1,0,ATD*C) + (2/vTE^3) * Cd * integral(profilefun_PS2,ATD*C,C); % 9.10
    BLDT_PS = BLET_PS/(3*He - 4);
    BLMT_PS = BLET_PS/He;
    
    BLT_SS(rr)=BLDT_SS;
    BLT_PS(rr)=BLDT_PS;
  end
  % Calculate blade profile and mixing losses - in Dickens (9.3) and also
  % Denton 1993 A3.7
  prof1(rr) = 2 * solid * ((BLMT_PS + BLMT_SS)/C) / cosd(ang2(rr)) * (0.5 * vel2(rr)^2 / Temp1(rr));  % mixing losses
  prof2(rr) = (solid*((BLDT_PS + BLDT_SS)/C + g.tTE) /cosd(ang2(rr)))^2 * (0.5 * vel2(rr)^2 / Temp1(rr)); % profile loss
  % Note - static temps used, but T0 ~ T with low M (~2% error for M=.3)
  % Note - can merge prof1 and prof2...
  
  % Calculate base pressure loss - Dickens (9.11) / Denton
  Cpb = -0.2; % Dickens p156 - based on turbines!
  base(rr) = - Cpb * solid * g.tTE / cosd(ang2(rr)) * (0.5 * vel2(rr)^2 / Temp1(rr)); % Base pressure loss

  % Shroud Loss - From Sungho Yoon Paper based on shroud losses for turbines in Denton
  hb = g.rc - g.rh; % span
  mm = (g.gap * d.Cc / hb) * ((tand(ang1(rr)))^2 - (tand(ang2(rr)))^2)^0.5; % Yoon 5 %(g * Cc / hn) * (abs((1/cosd(ang2))^2 - (tand(ang1))^2 ))^0.5 % # mm = (g_s * Cc / H_b) * np.sqrt( (1/np.cos(a2))**2 - (np.tan(a2))**2 ) 
  tip(rr) =  2 * mm * (1 - (tand(ang1(rr)) * sind(ang2(rr)) * cosd(ang2(rr))) ) * (0.5 * vel2(rr)^2 / Temp1(rr)); %
  % #tip =  mm * V2**2 * abs( 1 - (np.tan(a1) * np.sin(a2) * np.cos(a2)) ) / T #valid for compressible Yoon 6

  % Endwall Loss - comes from cd*V^3 argument but need to think about some more - important mechanisms I think
  endwallfun = @(x) (((vLE_SS + (vTE - vLE_SS)*x/Cx)).^3 + ((vLE_PS + (vTE - vLE_PS)*x/Cx)).^3)/2;
  % Hub  
  endwall1 = Cd * ((g.rm - hb/2)/g.rm)^3 / (vel1(rr)*Temp1(rr)*hb*cosd(ang1(rr))) * integral(endwallfun,0,Cx);
  % Casing
  endwall2 = Cd * ((g.rm + hb/2)/g.rm)^3 / (vel1(rr)*Temp1(rr)*hb*cosd(ang1(rr))) * integral(endwallfun,0,Cx);
  % Total
  endwall(rr) = (endwall1 + endwall2);
end

% Sum up losses
prof = prof1 + prof2;
Lrcalc = prof(1) + base(1) + tip(1) + endwall(1); %1st row
Lscalc = prof(2) + base(2) + 2*(1-k)*tip(2) + endwall(2); %2nd row, note no gap for stator - change this to toggle for CRDF

% Calculating efficiency from entropy
L.eta_s = q.cp*(q.T03s - q.T01) / ( q.cp*(q.T03 - q.T01) + q.T01*(Lrcalc+Lscalc));

%% Change in losses for while loop, Repack structures
L.deltaLr = (Lrcalc - L.Lr)/L.Lr;
L.deltaLs = (Lscalc - L.Ls)/L.Ls;
L.Lr = Lrcalc;
L.Ls = Lscalc;
L.rLoss.prof=prof(1); L.rLoss.base=base(1); L.rLoss.tip=tip(1); L.rLoss.endwall=endwall(1);
L.sLoss.prof=prof(2); L.sLoss.base=base(2); L.sLoss.tip=2*(1-k)*tip(2); L.sLoss.endwall=endwall(2);
L.rLoss.BLT_SS=BLT_SS(1); L.rLoss.BLT_PS=BLT_PS(1); L.sLoss.BLT_SS=BLT_SS(2); L.sLoss.BLT_PS=BLT_SS(2); 

g.c=C; g.s=s; g.Nb=Nb; %Return Re??

