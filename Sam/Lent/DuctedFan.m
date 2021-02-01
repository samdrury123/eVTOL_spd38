clear; close;

%% SPD Adding CRDF and tidying
bluebear = 0;
CRDF_flag = 0;

if bluebear == 1
    Th = 1500; % Thrust
    u0 = 2*25.7; % Flight speed
    D = 1.4; % Casing diameter
    alt = 0; % Altitude in m
    htrat = 0.3; % Hub-to-tip ratio
    AR = 4; % Blade aspect ratio
    Cdnom = 0.002; % Nominal dissipation coefficient for BL loss
    tTE = 0.003; % Trailing edge thickness, this should be a % of chord
    tmax = 0.02; % Max blade thickness, this should be a % of chord
    g = 0.001; % Shroud clearance in m
    Cc = 0.6; % Contraction coefficient for shroud loss calculation
else % Whittle eVTOL
    Th = 10; % Thrust
    u0 = 0; % Flight speed
    D = 0.06*2; % Casing diameter
    alt = 0; % Altitude in m
    htrat = 1/3; % Hub-to-tip ratio
    AR = 1.8; % Blade aspect ratio
    Cdnom = 0.0002; % Nominal dissipation coefficient for BL loss (Dickens p154)
    tTE = 0.0005; % Trailing edge thickness, in m
    tmax = 0.003; % Max blade thickness, in m
    g = 0.0005; % Shroud clearance in m
    Cc = 0.6; % Contraction coefficient for shroud loss calculation - Yoon p5 quotes Cd = 0.7
end

sigmalist = [0.8 0.9 1 1.1 1.2]; 
% sigmalist = [0.7 0.8 0.9 1 1.1 1.2]; 
% sigmalist = [0.9 1 1.1]; 
philist = [0.5 0.6 0.7 0.8 0.9];
% philist = [0.7 0.8];
DF=0.35; % Diffusion factor for Lieblein

[Tatm,patm,roatm,~] = Altitude(alt); % At sea level

% Gas Constants
gam = 1.4;
R = 287;
cp = gam*R / (gam-1);

Lr_count = zeros(size(philist,2),size(sigmalist,2));
Ls_count = zeros(size(philist,2),size(sigmalist,2));
deltaT_count = zeros(size(philist,2),size(sigmalist,2));
BL_count = zeros(size(philist,2),size(sigmalist,2),2);
deltavx_count = zeros(size(philist,2),size(sigmalist,2));


for pp = 1:size(philist,2)
    for ss = 1:size(sigmalist,2)
    tic
phi = philist(pp);
sigma = sigmalist(ss);

% Calculate geometry
rc = D/2; 
rh = rc*htrat;
rm = ((rc^2+rh^2)/2)^0.5;
A1 = pi*(rc^2 - rh^2);
A4 = A1 * sigma;

% Flight Mach Number and stagnation quantities
M0 = u0 / (gam*R*Tatm)^0.5;
p00 = patm ./ (1 + (gam-1)/2.*M0.^2).^-(gam/(gam-1));
T00 = Tatm ./ (1 + (gam-1)/2.*M0.^2).^-1;

% Set up Mach Number vector for lookups
M = 0.001:0.001:1;

% Non-dim impulse function for lookups
Fnd = (gam-1).^0.5 ./ gam .* (1+gam.*M.^2) ./ M .* (1 + (gam-1)/2.*M.^2).^-0.5;
% Non-dim mass flow functions for lookups
Mndp = gam / (gam-1)^0.5 .* M .* (1 + (gam-1)/2.*M.^2).^0.5;
Mndp0 = gam / (gam-1)^0.5 .* M .* (1 + (gam-1)/2.*M.^2).^(-0.5*(gam+1)/(gam-1));

% Assume exit pressure is atmospheric
p4 = patm;

% Guess T03 and initiate losses for iteration loops
T03 = 293; 
T04 = T03;
Lr = 0;     % Entropy change in rotor
Ls = 0;     % Entropy change in stator
deltaLr = 1;
deltaLs = 1;
Lr_c=0; Ls_c=0;
% Iteation loop for losses, improve while criterion (i.e. % change)
while abs(deltaLr) > 0.001/100 && abs(deltaLs) > 0.01/100
    Lr_c = Lr_c+1;
    Ls_c = Ls_c+1;

%% Iteration loop to find T03 and T04, again improve so not absolute value
deltaT=1;
deltaT_c = 0;
while abs(deltaT) > 0.001/100 %0.005 Now changing to percentagedeltaT_count
    deltaT_c = deltaT_c+1;
    % Non-dimensional momentum equation with impulse function (see 3A3 EP2,Q1)
    Thcalc = A4.* (Fnd .* Mndp .* p4 - patm) - Mndp .* (cp*T04).^-0.5 .* p4 .* A4 .* u0;
    M4 = interp1(Thcalc,M,Th);
    Fnd4 = (gam-1).^0.5 ./ gam .* (1+gam.*M4.^2) ./ M4 .* (1 + (gam-1)/2.*M4.^2).^-0.5;
    
    % Exit conditions 
    p04 = p4 ./ (1 + (gam-1)/2.*M4.^2).^-(gam/(gam-1));
    T4 = T04 ./ (1 + (gam-1)/2.*M4.^2);
    u4 = M4.*(gam*R*T4)^0.5;
    
    % Inlet and Exit duct isentropic
    p03 = p04;
    p01 = p00;
    T01 = T00;

    % Calculate inlet Mach Number, M1
    Mndp4 = gam / (gam-1)^0.5 .* M4 .* (1 + (gam-1)/2.*M4.^2).^0.5;
    Mndp01 = Mndp4 .* p4./p01 .* (T01./T04).^0.5 .* sigma;
    M1 = interp1(Mndp0,M,Mndp01);
    T1 = T01 .* (1 + 0.5*(gam-1).*M1.^2).^-1;
    u1 = M1 .* (gam*R*T1)^0.5;
    
    % Blade speed from flow coefficient
    Um = u1/phi;
    U1 = Um*(1-CRDF_flag/2);
    U2 = -Um*CRDF_flag/2;
    rpm = U1 / (2*pi*rm) * 60;
    Mtip = U1 * rc/rm /(gam*R*T1)^0.5; % U1 = -U2, Mtip equal

    % Stage loading
    psi = cp*(T04-T01) / Um^2;

    % Use loss to get T03
    T03calc = T01* exp((Lr + Ls + R*log(p03/p01))/cp);
    
    % Isentropic case
    T03s = T01 .* (p03/p01).^((gam-1)/gam);
    
    % Calculate change for while loop
    deltaT = (T03 - T03calc)/T03; %% Add relative term
    T03 = T03calc;
    T04 = T03;   
end
deltaT_count(pp,ss) = deltaT_c;
%% Meanline velocity triangles
% Upstream Rotor
v1 = u1;
vx1 = v1;
vt1 = 0;
vt1r = -U1;
v1r = (vx1^2 + vt1r^2)^0.5;
c1 = (gam*R*T1)^0.5;
M1r = v1r / c1;
a1 = 0;
a1r = atand(vt1r / vx1);
p1 = p01*(1+0.5*(gam-1)*M1^2)^(-gam/(gam-1));
p01r = p1*(1+0.5*(gam-1)*M1r^2)^(gam/(gam-1));
ro1 = p1/(R*T1); 
T01r = T1+0.5*v1r^2/cp;
Mndp01r = gam / (gam-1)^0.5 .* M1r .* (1 + (gam-1)/2.*M1r.^2).^(-0.5*(gam+1)/(gam-1));

% Downstream Rotor - Slightly tricky because A is constant not vx
vt2 = U1 * psi;% / (1-CRDF_flag/2);
vt2r(1) = vt2 - U1;
vt2r(2) = vt2 - U2;
p02r(1)=p01r*exp(-Lr/R);
T02r=T01r;
if deltaLr == 1
    vx2 = vx1; %this is a first guess
end
deltavx = 1;
% Iterate to find vx2 with constant area
% This loop seems to be quite unstable (hence the 0.01 smoothing) and could
% be improved with bisection method, or Newton-Raphson, etc. This should speed
% up programme and make it more stable. 

% SPD - updated to use the bisection method, starting with boundaries 2%
% either side of vx1 - stable, negligible speed change from adjusting
% convergence toleratnce for deltavx. Could set vx2_high = vx1, compressors
% increase static density so for const A, vx2<vx1 - negligible speed
% difference, only 1 iteration

deltavx_c = 0;
bisect_fac = 2/100;
vx2_low = (1-bisect_fac)*vx1; vx2_high = (1+bisect_fac)*vx1;

while abs(deltavx) > 0.0001/100 %0.00001
    vx2 = 0.5*(vx2_low + vx2_high);
    deltavx_c = deltavx_c+1;
   a2 = atand(vt2 / vx2);
   a2r(1) = atand(vt2r(1) / vx2);
   v2 = (vx2^2 + vt2^2)^0.5;
   v2r(1) = (vx2^2 + vt2r(1)^2)^0.5;
   Mndp02r(1) = Mndp01r * p01r/p02r(1) * cosd(a1r)/cosd(a2r(1));
   M2r(1) = interp1(Mndp0,M,Mndp02r(1));
   T2 = T02r - 0.5*v2r(1)^2 / cp;
   vx2calc = (M2r(1)^2*gam*R*T2 - vt2r(1)^2)^0.5;
   if (vx2calc - vx2) > 0
       vx2_low = vx2;
   else
       vx2_high = vx2;
   end
   deltavx = (vx2-vx2calc)/vx2;
%    deltavx2 = (vx2-vx2calc);
%    vx2 = vx2-0.1*deltavx2;
%    vx2 = vx2calc;
end
deltavx_count(pp,ss) = deltavx_c;

% deltavx2
% vx2
T02 = T2 + 0.5*v2^2/cp;
c2 = (gam*R*T2)^0.5;
M2 = v2/c2;
p2 = p02r(1)*(1+(gam-1)*0.5*M2r(1)^2)^(-gam/(gam-1));
p02 = p2*(1+(gam-1)*0.5*M2^2)^(gam/(gam-1));
ro2 = p2/(R*T2);

a2r(2) = atand(vt2r(2) / vx2);
v2r(2) = (vx2^2 + vt2r(2)^2)^0.5;
M2r(2) = v2r(2)/c2;
p02r(2) = p2*(1+(gam-1)*0.5*M2r(2)^2)^(gam/(gam-1));
T03r=T02r;
Mndp02r(2) = gam / (gam-1)^0.5 .* M2r(2) .* (1 + (gam-1)/2.*M2r(2).^2).^(-0.5*(gam+1)/(gam-1));


if CRDF_flag == 0 % Downstream Stator
    v2r(2)=v2; 
    p03=p02*exp(-Ls/R);
    a3 = 0; a3r=0;     % Change this to allow non-zero?
    vt3 = 0; vt3r=0;
    Mndp03 = Mndp4 * (p4 / p03) * sigma;
    M3 = interp1(Mndp0,M,Mndp03);
    T3 = T03 .* (1 + (gam-1)/2.*M3.^2).^-1;
    v3 = M3 * (gam*R*T3)^0.5;
    vx3 = v3; v3r = v3;
    p3 = p03*(1+0.5*(gam-1)*M3^2)^(-gam/(gam-1));
    ro3 = p3/(R*T3);
else % Downstream Rotor
    p03r=p02r(2)*exp(-Ls/R);
    a3 = 0;     % Change this to allow non-zero?
    vt3r = -U2;
    vt3 = 0;
    
    deltavx_2 = 1;
    vx3_low = (1-bisect_fac)*vx2; vx3_high = (1+bisect_fac)*vx2;
    while abs(deltavx_2) > 0.0001/100 %0.00001
        vx3 = 0.5*(vx3_low + vx3_high);
       v3 = vx3;
       v3r = (vx3^2 + vt3r^2)^0.5;
       a3r = atand(vt3r / vx3);
       Mndp03r = Mndp02r(2) * p02r(2)/p03r * cosd(a2r(2))/cosd(a3r);
       M3r = interp1(Mndp0,M,Mndp03r);
       T3 = T03r - 0.5*v3r^2 / cp;
       vx3calc = (M3r^2*gam*R*T3 - vt3r^2)^0.5;
       if (vx3calc - vx3) > 0
           vx3_low = vx3;
       else
           vx3_high = vx3;
       end
       deltavx_2 = (vx3-vx3calc)/vx3;
    end
    M3 = M3r*v3/v3r;
    p3 = p03r*(1+0.5*(gam-1)*M3r^2)^(-gam/(gam-1));
    p03 = p3*(1+(gam-1)*0.5*M3^2)^(gam/(gam-1));
    ro3 = p3/(R*T3);
end

%% Evaluate analytical velocity profiles and loss for each blade row, plus calculate shroud clearance and endwall losses

% Set up angles and speeds for each row
ang1 = [a1r a2r(2)];
ang2 = [a2r(1) a3r];
vel1 = [v1r v2r(2)];
vel2 = [v2r(1) v3r];
vel1t = [vt1r vt2r(2)];
vel2t = [vt2r(1) vt3r];
vel1x = [vx1 vx2]; 
vel2x = [vx2 vx3];
Temp1 = [T1 T2];
rho = [ro1 ro2];

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
  C = (rc - rh)/AR; % chord
  stag = (ang1(rr) + ang2(rr))/2; % stagger angle
  Cx = C * cosd(stag); % axial chord

  % Use Lieblein to get pitch-to-chord
  delvt = abs(vel1t(rr) - vel2t(rr));
  
  sc_rat = (DF - (1 - vel2(rr)./vel1(rr))) * 2 .* vel1(rr)  ./ delvt;

%   solid = 1/sc_rat;
  s = C*sc_rat;
%   Nb(rr) = 2*pi*rm / s;
  Nb(rr) = round(2*pi*rm / s);
%   s = 2 * pi * rm / (round(2*pi*rm / s));
  s = 2 * pi * rm / Nb(rr);
  solid = C/s;

  % Shape factors for BL equations
  H = 2;            % DT / MT - **How was this selected?**
  He = 4*H/(3*H-1); % ET / MT

  %Dissipation Coefficient as function of Reynolds
  Re = vel1(rr)*C*rho(rr) / 18e-6;
  Cd = Cdnom*(Re/500000)^-0.2; % Correction due to effect of low Re (Dickens p154)

  % First guess in boundary layer displacement thickness (BLDT)
  BLDT_SS = 0.0001; BLDT_SS_g = 0;
  BLDT_PS = 0.0001; BLDT_PS_g = 0;

  % Iterate velocity profiles and boundary layers. This is explained in Dickens thesis
  BL_c = 0;
  while abs((BLDT_SS - BLDT_SS_g)/BLDT_SS) > 1e-10 && abs((BLDT_PS - BLDT_PS_g)/BLDT_PS) > 1e-10
    BL_c = BL_c + 1;
    BLDT_SS_g = BLDT_SS;
    BLDT_PS_g = BLDT_PS;
   
    vTE = vel2(rr) * (1 - solid * ( (tTE + BLDT_SS + BLDT_PS) / (cosd(ang2(rr)) * C) ))^-1; % Dickens 9.2
    vLE = vel1(rr) * (0.98 + 0.5 * solid * tmax/C * vel1x(rr) / delvt); % Dickens 9.3
   
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
    BL_count(pp,ss,rr) = BL_c;
  % Calculate blade profile and mixing losses - in Dickens (9.3) and also
  % Denton 1993 A3.7
  prof1(rr) = 2 * solid * ((BLMT_PS + BLMT_SS)/C) / cosd(ang2(rr)) * (0.5 * vel2(rr)^2 / Temp1(rr));  % mixing losses
  prof2(rr) = (solid*((BLDT_PS + BLDT_SS)/C + tTE/C)/cosd(ang2(rr)))^2 * (0.5 * vel2(rr)^2 / Temp1(rr)); % profile loss
  % Note - static temps used, but T0 ~ T with low M (~2% error for M=.3)
  % Note - can merge prof1 and prof2...
  
  % Calculate base pressure loss - Dickens (9.11) / Denton
  Cpb = -0.2; % Dickens p156 - based on turbines!
  base(rr) = - Cpb * solid * (tTE/C) / cosd(ang2(rr)) * (0.5 * vel2(rr)^2 / Temp1(rr)); % Base pressure loss

  % Shroud Loss - From Sungho Yoon Paper based on shroud losses for turbines in Denton
  hb = rc - rh;
  mm = (g * Cc / hb) * ((tand(ang1(rr)))^2 - (tand(ang2(rr)))^2)^0.5; % Yoon 5 %(g * Cc / hn) * (abs((1/cosd(ang2))^2 - (tand(ang1))^2 ))^0.5 % # mm = (g_s * Cc / H_b) * np.sqrt( (1/np.cos(a2))**2 - (np.tan(a2))**2 ) 
  tip(rr) =  2 * mm * (1 - (tand(ang1(rr)) * sind(ang2(rr)) * cosd(ang2(rr))) ) * (0.5 * vel2(rr)^2 / Temp1(rr)); % Yoon 6
   %       #tip =  mm * V2**2 * abs( 1 - (np.tan(a1) * np.sin(a2) * np.cos(a2)) ) / T #valid for compressible
 
  % Endwall Loss - comes from cd*V^3 argument but need to think about some more - important mechanisms I think
  endwallfun = @(x) (((vLE_SS + (vTE - vLE_SS)*x/Cx)).^3 + ((vLE_PS + (vTE - vLE_PS)*x/Cx)).^3)/2;
  % Hub  
  endwall1 = Cd * ((rm - hb/2)/rm)^3 / (vel1(rr)*Temp1(rr)*hb*cosd(ang1(rr))) * integral(endwallfun,0,Cx);
  % Casing
  endwall2 = Cd * ((rm + hb/2)/rm)^3 / (vel1(rr)*Temp1(rr)*hb*cosd(ang1(rr))) * integral(endwallfun,0,Cx);
  % Total
  endwall(rr) = (endwall1 + endwall2);
end

% Sum up losses
prof = prof1 + prof2;
Lrcalc = prof(1) + base(1) + tip(1) + endwall(1); %1st row
Lscalc = prof(2) + base(2) + tip(2) + endwall(2); %2nd row

% Change in losses for while loop
deltaLr = (Lrcalc - Lr)/Lr;
deltaLs = (Lscalc - Ls)/Ls;
Lr = Lrcalc;
% vx2 = vx2-0.01*deltavx;
Ls = Lscalc;

end
Lr_count(pp,ss) = Lr_c; Ls_count(pp,ss) = Ls_c;

display(['phi = ' num2str(phi) '  sigma = ' num2str(sigma) '   Time = ' num2str(toc,3) '   Tip Loss = ' num2str(endwall./[Lr Ls])])


% Finished loss loop so everything defined. The rest is gathering together
% outputs

% Power
mdot = Mndp01 * A1 * p01 / (cp*T01)^0.5;
P = mdot * cp*(T03-T01);

% Reaction
Reaction = (T2-T1)/(T3-T1);

% Tip Mach Number
Mtip2 = (rpm/60*2*pi*rc) / (gam*R*T1)^0.5;

% Propulsive efficiency and hover figure of merit
Fr = 2*u0 / (u0+u4);
Mf = Th / P * sqrt(Th / (2*roatm*A1));

% Fan Efficiency
eta = (T03s - T01) / (T03 - T01);

% Inlet streamtube capture 
% Currently assumes isentropic inlet duct
Mndp00 = gam / (gam-1)^0.5 .* M0 .* (1 + (gam-1)/2.*M0.^2).^(-0.5*(gam+1)/(gam-1));
A0 = A1 * Mndp01 / Mndp00;
CAR = A0/A1;
rnrat = 1.2;
Ain = pi*(rc*rnrat)^2;
MFR = A0/Ain;

% Advance ratio and loading coefficient
J = u0 / (rpm/60 * 2*rc);
CT = Th / (roatm * (rpm/60)^2 * (2*rc)^4);


% Hacking together arrays to plot
phis(pp,ss) = phi;
sigmas(pp,ss) = sigma;
etas(pp,ss) = eta;  % Fan efficiency 
Rs(pp,ss) = Reaction; % Reaction
rpms(pp,ss) = rpm; % rpm
Nr(pp,ss) = Nb(1); % Blade count row 1
Ns(pp,ss) = Nb(2); % Blade count row 2
Frs(pp,ss) = Fr; % Propulsive efficiency
psis(pp,ss) = psi;
FOMs(pp,ss) = Mf;

    
    end
end

counts = [Lr_count,Ls_count;deltaT_count,deltavx_count;BL_count(:,:,1), BL_count(:,:,2)]
disp("Lr_count,Ls_count")
disp("deltaT_count,deltavx_count")
disp("BL_count(:,:,1), BL_count(:,:,2)")

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
if u0 == 0
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





