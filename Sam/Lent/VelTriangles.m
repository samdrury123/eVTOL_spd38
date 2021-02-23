%% ---- FUNCTION: VELOCITY TRIANGLES ---- %%
function [d,g,a,q,L] = VelTriangles(d,g,q,L,k,dev)

% Function calculates velocities and angles for conventional or
% counter-rotating fan, as well as the losses in each blade row
% Currently works with alpha1=alpha3=0, single design point
% Inputs:
% - d.phi: flow coefficient 
% - d.psi: stage loading coefficient
% - d.Um: Meanline blade speed
% - k: rpm split between blade A and blade B, k=1 is conventional
% - g: geometry of fan
% - q: aerodynamic quantities (M,T,p)
% - L: L1,L2 entropy gain over blade 1,2
% Outputs:
% - d: updated with blade rpms
% - g: geometry updated with spanwise chord, pitch, blade metal angles, etc
%      chi are arrays with blade metal angles and have dimensions:
%      - radius
%      - 1: LE, 2: TE
%      - 1: Blade 1 2: Blade 2
% - a: Structure with aerodynamic velocities and angles, plus relative reaction
%      3d arrays of velocities and angles have dimensions: 
%      - radius
%      - position (upB1, dwnB1, dwnB2)
%      - reference frame (absolute, B1, B2)
% - q: updated with M,T,p,rho
% - L: updated L1,L2 for iteration, breakdown of loss

%% Unpack
u1=d.u1; Um=d.Um; psi=d.psi; sigma=d.sigma;
rm=g.rm; rh=g.rh; rc=g.rc;
gam=q.gam; R=q.R; cp=q.cp;
p01=q.p01; p4=q.p4; T03=q.T03; T1=q.T1; M1=q.M1; Mndp4=q.Mndp4;

% Set up Mach Number vector for lookups
Mach = 0.001:0.001:1;
Mndp0 = gam / (gam-1)^0.5 .* Mach .* (1 + (gam-1)/2.*Mach.^2).^(-0.5*(gam+1)/(gam-1));

% Set up radial positions, trying to make rm close to a point in r
Nr = 6; % Number of radial slices
r = linspace(rh,rc,Nr)';
[~,imid] = min(abs((r-rm)/(rc-rh)));
r(imid) = rm;
g.r=r;
% thresh=0.1; rmin=1;
% while rmin>thresh
%     Nr=Nr+1;
%     r = linspace(rh,rc,Nr)';
%     [rmin,i] = min(abs((r-rm)/(rc-rh)));
% end

% Split blade speeds between rows (CRDF)
Um = r ./ rm * Um;
U1 = Um*k;      rpm1 = d.rpm*k;
U2 = -Um*(1-k); rpm2 = -d.rpm*(1-k);

% Set up spanwise stage loading coefficient, enthalpy and vtheta
psi = psi .* (rm ./ r).^ (g.n+1);

% Set Vx constant across span at inlet
vx(1:Nr,1) = ones(size(r)) * u1;

% Calculate velocities and angles upstream of blade 1
vth(1:Nr,1,1) = 0; % alpha1 = 0
vth(:,1,2) = vth(:,1,1) - U1;
vth(:,1,3) = vth(:,1,1) - U2;

% Calculate velocities and angles downstream of blade 1, including axial
% velocity distribution from radial equilibrium equation
vth(:,2,1) = U1 .* psi / k; %(Um*psi)
vth(:,2,2) = vth(:,2,1) - U1;
vth(:,2,3) = vth(:,2,1) - U2;

% Calculate velocities and angles downstream of blade 2
vth(:,3,1) = 0; % alpha3 = 0
vth(:,3,2) = vth(:,3,1) - U1;
vth(:,3,3) = vth(:,3,1) - U2;

% Calculate angles and velocities upstream
for tt = 1:3 % reference frames
  vx(:,1,tt) = vx(:,1,1);
  alpha(:,1,tt) = atand(vth(:,1,tt) ./ vx(:,1,tt));
  v(:,1,tt) = vx(:,1,tt) ./ cosd(alpha(:,1,tt));
end
% disp('starting vx2');
%% Solve for vx2 - have radial positions for all (incl M) but p,p0,T,T0,rho
% Assume radial velocities are negligible
% Upstream Rotor
% v1 = u1;
% vx1 = v1;
% vt1 = 0;
% vt1r = -U1;
% v1r = (vx1^2 + vt1r^2)^0.5;
c1 = (gam*R*T1)^0.5;
M(1:Nr,1,1) = M1;
M(1:Nr,1,2) = v(:,1,2)/c1;
M(1:Nr,1,3) = v(:,1,3)/c1;
% M1r = v1r / c1;
% a1 = 0;
% a1r = atand(vt1r / vx1);
p1 = p01*(1+0.5*(gam-1)*M(imid,1,1)^2)^(-gam/(gam-1));
p01r = p1*(1+0.5*(gam-1)*M(imid,1,2)^2)^(gam/(gam-1));
ro1 = p1/(R*T1); 
T01r = T1+0.5*v(imid,1,2)^2/cp;
Mndp01r = gam / (gam-1)^0.5 .* M(imid,1,2) .* (1 + (gam-1)/2.*M(imid,1,2).^2).^(-0.5*(gam+1)/(gam-1));

% Downstream Rotor - Slightly tricky because A is constant not vx
% vt2 = Um * psi; 
% vt2r(1) = vt2 - U1;
% vt2r(2) = vt2 - U2;
p02r(1)=p01r*exp(-L.L1/R);
T02r=T01r;
if L.deltaL1 == 1
    vx(:,2,:) = vx(:,1,:); %this is a first guess
end

% Iterate to find vx2 with constant area
% This loop seems to be quite unstable (hence the 0.01 smoothing) and could
% be improved with bisection method, or Newton-Raphson, etc. This should speed
% up programme and make it more stable. 

% SPD - updated to use the bisection method, starting with boundaries 2%
% either side of vx1 - stable, negligible speed change from adjusting
% convergence toleratnce for deltavx. Could set vx2_high = vx1, compressors
% increase static density so for const A, vx2<vx1 - negligible speed
% difference, only 1 iteration

bisect_fac = 2/100;
vx2_low = (1-bisect_fac)*vx(:,1,1); vx2_high = (1+bisect_fac)*vx(:,1,1);
deltavx2 = 1;

while abs(deltavx2) > 0.01/100 %0.00001
%     disp('vx2 loop')
    vx(:,2,1) = 0.5*(vx2_low + vx2_high);
    alpha(:,2,1) = atand(vth(:,2,1) ./ vx(:,2,1));
    alpha(:,2,2) = atand(vth(:,2,2) ./ vx(:,2,1));
    v(:,2,1) = (vx(:,2,1).^2 + vth(:,2,1).^2) .^0.5;
    v(:,2,2) = (vx(:,2,1).^2 + vth(:,2,2).^2) .^0.5;
    Mndp02r(1) = Mndp01r * p01r/p02r(1) * cosd(alpha(imid,1,2))/cosd(alpha(imid,2,2));
    M(:,2,2) = interp1(Mndp0,Mach,Mndp02r(1)) .* v(:,2,2)./v(imid,2,2);
    T2 = T02r - 0.5*v(imid,2,2)^2 / cp;
    vx2calc = (M(:,2,2).^2*gam*R*T2 - vth(:,2,2).^2).^0.5;
    if (vx2calc(imid) - vx(imid,2,1)) > 0
       vx2_low = vx(:,2,1);
    else
       vx2_high = vx(:,2,1);
    end
    deltavx2 = (vx2calc(imid) - vx(imid,2,1))/vx(imid,2,1);
%     if mean(vx2calc - vx(:,2,1)) > 0
%        vx2_low = vx(:,2,1);
%     else
%        vx2_high = vx(:,2,1);
%     end
%     deltavx2 = mean(vx(:,2,1)-vx2calc)/mean(vx(:,2,1));
end


%% Radial Equilibrium between blade rows (Position 2)

% Use Euler work equation to calculate change in enthalpy
dh0 = U1 .* vth(:,2,1);

% Initialise temperature, pressure, density arrays
T(1:Nr,1) = T1;
T0(1:Nr,1,1) = q.T01;
T0(:,3,1) = T03; T0(:,1,2) = T01r; T0(:,2,3) = T02r;
p(1:Nr,1) = p1;
p0(1:Nr,1,1) = p01; p0(:,1,2) = p01r; p0(:,2,2) = p02r(1); 
ro(1:Nr,1) = ro1;

% Add in second row quantities
p(imid,2) = p02r(1)*(1+(gam-1)*0.5*M(imid,2,2)^2)^(-gam/(gam-1));
T(imid,2) = T2;
ro(imid,2) = p(imid,2)/(R*T(imid,2));

% Numerically solve Simple Radial Equilibrium equation

for ii = imid:(Nr-1)
    dr = r(ii+1) - r(ii); % Radial increment
    vx(ii+1,2,1)     = vx(ii,2,1) + (dh0(ii+1) - dh0(ii))/vx(ii,2,1) - (vth(ii,2,1)/vx(ii,2,1)) * (vth(ii+1,2,1) - vth(ii,2,1)) - ((vth(ii,2,1)).*(vth(ii+1,2,1)) / vx(ii,2,1)) * (dr/r(ii));
    p(ii+1,2)        = p(ii,2) + ro(ii,2) * (vth(ii,2,1)*vth(ii+1,2,1)) * (dr/r(ii));
    T(ii+1,2)        = T(ii,2) + (p(ii+1,2) - p(ii,2)) / (ro(ii,2)*cp);
    ro(ii+1,2)       = p(ii+1,2) / (R*T(ii+1,2));
end
for ii = imid:-1:2
    dr = r(ii) - r(ii-1);
    vx(ii-1,2,1)     = -(-vx(ii,2,1) + (dh0(ii) - dh0(ii-1))/vx(ii,2,1) - (vth(ii,2,1)/vx(ii,2,1)) * (vth(ii,2,1) - vth(ii-1,2,1)) - ((vth(ii,2,1)).*(vth(ii-1,2,1)) / vx(ii,2,1)) * (dr/r(ii)));
    p(ii-1,2)        = -(-p(ii,2) + ro(ii,2) * (vth(ii,2,1)*vth(ii-1,2,1)) * (dr/r(ii)));
    T(ii-1,2)        = -(-T(ii,2) + (p(ii,2) - p(ii-1,2)) / (ro(ii,2)*cp));
    ro(ii-1,2)       = -(-p(ii-1,2) / (R*T(ii-1,2)));
end
vx(:,2,2) = vx(:,2,1); vx(:,2,3) = vx(:,2,1);

% Reset alpha and v according to vx spanwise distribution
alpha(:,2,:) = atand(vth(:,2,:) ./ vx(:,2,:));
v(:,2,:) = (vx(:,2,:).^2 + vth(:,2,:).^2) .^0.5;


%% Solve for vx3
% alpha(:,2,3) = atand(vth(:,2,3)./vx(:,2,1));
% v(:,2,3) = (vx(:,2,1).^2 + vth(:,2,3).^2).^0.5;
% v(:,2,1) = (vx(:,2,1).^2 + vth(:,2,1).^2) .^0.5;
% a2r(2) = atand(vt2r(2) / vx2);
% v2r(2) = (vx2^2 + vt2r(2)^2)^0.5;
% M2r(2) = v2r(2)/c2;
% p02r(2) = p2*(1+(gam-1)*0.5*M(imid,2,3)^2)^(gam/(gam-1));
Mndp02r(2) = gam / (gam-1)^0.5 .* M(imid,2,3) .* (1 + (gam-1)/2.*M(imid,2,3).^2).^(-0.5*(gam+1)/(gam-1));

% Correct stagnation quantities
T0(:,2,:) = T(:,2) + 0.5*v(:,2,:).^2./cp;
c2 = (gam*R*T(:,2)).^0.5;
M(:,2,1) = v(:,2,1)./c2;
M(:,2,3) = v(:,2,3)./c2;
p0(:,2,1:3) = p(:,2).*(1+(gam-1).*0.5.*M(:,2,:).^2).^(gam/(gam-1));

% Zero exit swirl, hence there is no static or stagnation variation across the exit plane
alpha(:,3,1) = 0; 

if rpm2 == 0 % Downstream Stator
    p0(:,3,1)=p0(imid,2,1)*exp(-L.L2/R); p0(:,3,3)=p0(:,3,1); % Fudge this to be equal across exit plane, loss evaluated at mean line
    T0(:,3,3) = T0(:,3,1);
    alpha(:,3,3) = 0;
    Mndp03 = Mndp4 * (p4 / p0(imid,3,1)) * sigma;
    M(:,3,1) = interp1(Mndp0,Mach,Mndp03);
    M(:,3,3) = M(:,3,1);
    T(:,3) = T0(:,3,1) .* (1 + (gam-1)/2.*M(:,3,1).^2).^-1;
    c3 = (gam*R*T(:,3)).^0.5;
    v(:,3,1) = M(:,3,1) .* c3;
    v(:,3,3) = v(:,3,1);
    vx(:,3,1) = v(:,3,1); vx(:,3,2) = vx(:,3,1); vx(:,3,3) = vx(:,3,1);
    p(:,3) = p0(:,3,1).*(1+0.5.*(gam-1).*M(:,3,1).^2).^(-gam/(gam-1));
    ro(:,3) = p(:,3)./(R*T(:,3));
else % Downstream Rotor
    p0(:,3,3) = p0(imid,2,3)*exp(-L.L2/R); % Fudge this to be equal across exit plane, loss evaluated at mean line
%     p03r=p02r(2)*exp(-L.L2/R);
%     T03r=T02r;
    T0(:,3,3) = T0(:,2,2); % This should be corrected for k
    
    deltavx3 = 1;
    vx3_low = (1-bisect_fac)*vx(:,2,1); vx3_high = (1+bisect_fac)*vx(:,2,1);
    while abs(deltavx3) > 0.0001/100 %0.00001
%         disp('vx3 loop')
        vx(:,3,1) = 0.5*(vx3_low + vx3_high);
        v(:,3,1) = vx(:,3,1);
        v(:,3,3) = (vx(:,3,1).^2 + vth(:,3,3).^2).^0.5;
        alpha(:,3,3) = atand(vth(:,3,3)./vx(:,3,1));
%         v3r = (vx3^2 + vt3r^2)^0.5;
%         a3r = atand(vt3r / vx3);
        Mndp03r = Mndp02r(2) * p0(imid,2,3)/p0(imid,3,3) * cosd(alpha(imid,2,3))/cosd(alpha(imid,3,3));
        M(:,3,3) = interp1(Mndp0,Mach,Mndp03r)* v(:,3,3)./v(imid,3,3);
        T(:,3) = T0(imid,3,3) .* (1 + (gam-1)/2.*M(imid,3,3).^2).^-1;%T0(imid,3,3) - 0.5.*v(imid,3,3).^2 ./ cp; % Forcing T3 constant across span
        vx3calc = (M(:,3,3).^2.*gam.*R.*T(:,3) - vth(:,3,3).^2).^0.5;
        if (vx3calc(imid) - vx(imid,3,1)) > 0
            vx3_low = vx(:,3,1);
        else
            vx3_high = vx(:,3,1);
        end
        deltavx3 = (vx3calc(imid) - vx(imid,3,1))/vx(imid,3,1);
    end
    vx(:,3,1) = v(:,3,1); vx(:,3,2) = vx(:,3,1); vx(:,3,3) = vx(:,3,1);
    c3 = (gam*R*T(:,3)).^0.5;
    p(:,3) = p0(imid,3,3)*(1+0.5*(gam-1)*M(imid,3,3)^2)^(-gam/(gam-1));
    p0(:,3,1) = p(imid,3)*(1+(gam-1)*0.5*M(imid,3,1)^2)^(gam/(gam-1));
    ro(:,3) = p(:,3)./(R*T(:,3));
end

% Backfilling arrays
alpha(:,3,2) = atand(vth(:,3,2) ./ vx(:,3,2));
v(:,3,2) = vx(:,3,2) ./ cosd(alpha(:,3,2));
M(:,3,2) = v(:,3,2)./c3;
M(:,3,1) = v(:,3,1)./c3;
T0(:,3,2) = T(:,3)+0.5*v(:,3,2).^2./cp;
T0(:,1,3) = T(:,1)+0.5*v(:,1,3).^2./cp;
p0(:,3,2) = p(:,3).*(1+(gam-1)*0.5*M(:,3,2).^2).^(gam/(gam-1));
p0(:,1,3) = p(:,1).*(1+(gam-1)*0.5*M(:,1,3).^2).^(gam/(gam-1));

% de Haller number = V2/V1 in relevant frame
d.DHall(:,1) = v(:,2,2)./v(:,1,2);
d.DHall(:,2) = v(:,3,3)./v(:,2,3);
d.DH(1) = v(imid,2,2)./v(imid,1,2);
d.DH(2) = v(imid,3,3)./v(imid,2,3);

%% Calculating losses, chords, chis for each blade row 11/02

% Initialise loss components to be stored
prof1 =     [0,0];
prof2 =     [0,0];
baseT =      [0,0];
tip =       [0,0];
endwall =   [0,0];

% Loop through each blade row
for rr=1:2
    % Gather inputs
    vup = v(:,rr,rr+1);
    vdwn = v(:,rr+1,rr+1);
    vupth = vth(:,rr,rr+1);
    vdwnth = vth(:,rr+1,rr+1);
    vxup = vx(:,rr,rr+1);
    Tup = T0(:,rr,rr+1);
    Tdwn = T0(:,rr+1,rr+1);
    % Multiplier to get angles in right direction while doing deviation calcs
    bl = round(rr-1.5);
    % Get angles into "cascade" frame, i.e. all positive, for deviation calcs
    aup = bl * alpha(:,rr,rr+1);
    adwn =  bl * alpha(:,rr+1,rr+1);

    % Blade Geometry
    C = (rc - rh)/g.AR; % initial guess, chord at mean line
    stag(:,rr) = (aup + adwn)/2;
    %   stag = (ang1(rr) + ang2(rr))/2; % stagger angle
    %   Cx = C * cosd(stag); % axial chord
    Cx(:,rr) = C * cosd(stag(:,rr)); % axial chord

    % Use Lieblein to get pitch-to-chord
    %   delvt = abs(vel1t(rr) - vel2t(rr));
    delvt(:,rr) = abs(vupth - vdwnth);

    %   sc_rat = (d.DF - (1 - vel2(rr)./vel1(rr))) * 2 .* vel1(rr)  ./ delvt;
    sc_rat(:,rr) = (d.DF - (1 - vdwn./vup)) * 2 .* vup  ./ delvt(:,rr);
    s(:,rr) = C*sc_rat(:,rr);
    g.Nb(rr) = round(2*pi*rm / s(imid,rr));
    if g.Nb(rr)<1
        g.Nb(rr) = NaN;
        disp('Negative s/c')
        L.deltaL1 = 1e-10; % satisfy criterion to exit loop
        L.deltaL2 = 1e-10;
        a = 0;
        return
    elseif g.Nb(rr)>60
        g.Nb(rr) = NaN;
        disp(['Number of blades in blade ' num2str(rr) ' exceeds 60'])
        L.deltaL1 = 1e-10; % satisfy criterion to exit loop
        L.deltaL2 = 1e-10;
        a = 0;
        return
    end
    %   s = 2 * pi * rm / Nb(rr);
    %   solid = C/s;
    g.s(:,rr) = 2 * pi * r / g.Nb(rr);
    solid(:,rr) = C./g.s(:,rr);

    %% Chi, Chord, DF
%     dev = 1; % d is whether to design with constant deviation across span 
    %        - d = 0: chord is set to given constant diffusion factor across span
    %        - d = 1:  chord is set to give constant deviation across span based on midspan value
    if dev==0 
       % Calculate chord and deviation to give constant diffusion factor across span
       % Calculate chord based on Lieblein
       g.c(1:Nr,rr) = g.s(:,rr) ./ sc_rat(:,rr);

       % Calculate blade metal angles with deviation using Howell and Carter
       % Blade LE metal angle, zero incidence 
       chi(:,1,rr) = aup;
       g.chi(:,1,rr) = bl .* chi(:,1,rr);

       % Howell and Carter
       ac = 0.5; % Maximum camber position to chord ratio
       m = 0.23*(2*ac).^2 + adwn / 500;
       chi(:,2,rr) = (adwn - m .* aup .* (g.s(:,rr) ./ g.c(:,rr)).^ 0.5) ./ (1 - m.*(g.s(:,rr) ./ g.c(:,rr)).^0.5);
       g.chi(:,2,rr) = bl .* chi(:,2,rr);
       devi(:,rr) = chi(:,2,rr) - adwn; 

    else
       % Calculate chord to given constant deviation across span based on midline value
       % Blade LE metal angle, zero incidence
       chi(:,1,rr) = aup;
       g.chi(:,1,rr) = bl .* chi(:,1,rr);

       % Calculate meanline chord based on Lieblein
       cmid = g.s(imid,rr) ./ sc_rat(imid,rr) ;

       % Calculate midspan deviation
       ac=0.5;
       m = 0.23*(2*ac).^2 + adwn(imid) / 500;
       chimid = (adwn(imid) - m .* aup(imid) .* (g.s(imid,rr) ./ cmid).^ 0.5) ./ (1 - m.*(g.s(imid,rr) ./ cmid).^0.5);
       devi(rr) =  (adwn(imid) - chimid);

       % Calculate chord based on constant deviation
       m = 0.23*(2*ac).^2 + adwn / 500;
       g.c(1:Nr,rr) = g.s(:,rr) .* (devi(rr) ./ (m .* (aup + devi(rr) - adwn) )).^-2;
       chi(:,2,rr) = adwn - devi(rr);
       g.chi(:,2,rr) = bl .* chi(:,2,rr);

    end

%     if isfield(g2,'cman')==1
%       % Adjust chord with manual input
%       g.c(:,rr) = g.c(:,rr) .* linspace(g2.cman(rr,1),g2.cman(rr,2), size(g.c,1))';
%     end

    % Recalculate deviation and blade metal angles
    % Howell and Carter
    ac = 0.5; % Maximum camber position to chord ratio
    m = 0.23*(2*ac).^2 + adwn / 500;
    chi(:,2,rr) = (adwn - m .* aup .* (g.s(:,rr) ./ g.c(:,rr)).^ 0.5) ./ (1 - m.*(g.s(:,rr) ./ g.c(:,rr)).^0.5);
    g.chi(:,2,rr) = bl .* chi(:,2,rr);

    % Recalculate Diffusion Factor based on chords calculated 
    a.DFspan(1:Nr,rr) = (1 - vdwn./vup) + abs(vupth - vdwnth) ./ (2.*vup) .* (g.s(:,rr) ./ g.c(:,rr));
    if a.DFspan > 0.6
        disp('Warning: DF > 0.6')
    end

    %% Calculate Losses
    
    % Shape factors for BL equations
    H = 2.7;            % DT / MT - **How was this selected?**
    He = 4*H/(3*H-1); % ET / MT

    %Dissipation Coefficient as function of Reynolds - TURNED OFF 11/02
%     Re = vel1(rr)*C*rho(rr) / 18e-6;
    Cd = d.Cdnom;%*(Re/500000)^-0.2; % Correction due to effect of low Re (Dickens p154)

    %% Profile loss
    % First guess in boundary layer displacement thickness (BLDT)
    %   BLDT_SS = 0.0001; BLDT_SS_g = 0;
    %   BLDT_PS = 0.0001; BLDT_PS_g = 0;


    % Iterate velocity profiles and boundary layers. This is explained in Dickens thesis
    for ii=1:Nr
    BLDT_SS = 0.0001; BLDT_SS_g = 0;
    BLDT_PS = 0.0001; BLDT_PS_g = 0;
    while abs((BLDT_SS - BLDT_SS_g)/BLDT_SS) > 1e-10 && abs((BLDT_PS - BLDT_PS_g)/BLDT_PS) > 1e-10
    %   while abs((BLDT_SS - BLDT_SS_g)/BLDT_SS) > 1e-10 && abs((BLDT_PS - BLDT_PS_g)/BLDT_PS) > 1e-10
    %     BLDT_SS_g = BLDT_SS;
    %     BLDT_PS_g = BLDT_PS;
    %    
    %     vTE = vel2(rr) * (1 - solid * ( (g.tTE + (BLDT_SS + BLDT_PS)/C) / cosd(ang2(rr))) ).^-1; % Dickens 9.2
    %     vLE = vel1(rr) * (0.98 + 0.5 * solid * g.tmax * vel1x(rr) / delvt); % Dickens 9.3
    %    
    %     ATD = 0; % Roof-top parameter set to 0 - Dickens p147/8
    %     delv = 3 * vel1x(rr) * delvt / (solid * cosd(stag) * ((2*vLE*(1+2*ATD))+vTE*(1-ATD))); % Dickens 9.1 
    %     vLE_SS = vLE + delv;
    %     vLE_PS = vLE - delv;
    %    
    %     profilefun_SS1 =  @(z) (z*0 + vLE_SS).^3;
    %     profilefun_SS2 = @(z) (((vLE_SS - (ATD * vTE)) / (1-ATD)) - (z/C*(vLE_SS - vTE)/(1-ATD))).^3; % Dickens D.5
    %     BLET_SS = (2/vTE^3) * Cd * integral(profilefun_SS1,0,ATD*C) + (2/vTE^3) * Cd * integral(profilefun_SS2,ATD*C,C); % Dickens 9.10
    %     % Check Dickens 9.12 to improve this, if using the density ratio - not
    %     % an issue for low M, but 8% variation possibly
    %     BLDT_SS = BLET_SS/(3*He - 4); % Dickens 9.8
    %     BLMT_SS = BLET_SS/He;
    %    
    %     profilefun_PS1 =  @(z) (z*0 + vLE_PS).^3;
    %     profilefun_PS2 = @(z) (((vLE_PS - (ATD * vTE))/(1-ATD)) - (z/C*(vLE_PS - vTE)/(1-ATD))).^3; % Dickens D.6
    %     BLET_PS= (2/vTE^3) * Cd * integral(profilefun_PS1,0,ATD*C) + (2/vTE^3) * Cd * integral(profilefun_PS2,ATD*C,C); % 9.10
    %     BLDT_PS = BLET_PS/(3*He - 4);
    %     BLMT_PS = BLET_PS/He;
    %     
    %     BLT_SS(rr)=BLDT_SS;
    %     BLT_PS(rr)=BLDT_PS;


    BLDT_SS_g = BLDT_SS;
    BLDT_PS_g = BLDT_PS;
    C = g.c(ii,rr);

    vTE = vdwn(ii) * (1 - solid(ii,rr) * ( (g.tTE + BLDT_SS + BLDT_PS)/C / cosd(adwn(ii))) ).^-1; % Dickens 9.2
    vLE = vup(ii) * (0.98 + 0.5 * solid(ii,rr) * g.tmax * vxup(ii) / delvt(ii,rr)); % Dickens 9.3
    
    ATD = 0; % Roof-top parameter set to 0 - Dickens p147/8
    delv = 3 * vxup(ii) * delvt(ii,rr) / (solid(ii,rr) * cosd(stag(ii,rr)) * ((2*vLE*(1+2*ATD))+vTE*(1-ATD))); % Dickens 9.1 
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
    
    end
    
    % Storing velocity profiles for plotting later
    x = linspace(0,C,10);
    a.Vss(ii,1:size(x,2),rr) = profilefun_SS2(x).^(1/3);
    a.Vps(ii,1:size(x,2),rr) = profilefun_PS2(x).^(1/3);
    
    BLT_SS(ii,rr)=BLDT_SS;
    BLT_PS(ii,rr)=BLDT_PS;
    p1(ii,rr) = 2 * solid(ii,rr) * ((BLMT_PS + BLMT_SS)/C) / cosd(adwn(ii)) * (0.5 * vdwn(ii)^2 / Tup(ii));  % mixing losses - this separation is explained just below Denton A3.7
    p2(ii,rr) = (solid(ii,rr)*(BLDT_PS + BLDT_SS + g.tTE)/C /cosd(adwn(ii)))^2 * (0.5 * vdwn(ii)^2 / Tup(ii)); % profile loss
    vLE_SSe(ii) = vLE_SS;
    vLE_PSe(ii) = vLE_PS;
    vTEe(ii) = vTE;
    
%     Reth_ps = vup(ii)*BLMT_PS*ro(ii,rr) / 18e-6
%     Reth_ss = vup(ii)*BLMT_SS*ro(ii,rr) / 18e-6
    end
    % Calculate blade profile and mixing losses - in Dickens (9.3) and also
    % Denton 1993 A3.7
    prof1(rr) = ((p1(1,rr) + p1(Nr,rr))/2 + sum(p1(2:(Nr-1),rr)) )/5; % finite volume method
    prof2(rr) = ((p2(1,rr) + p2(Nr,rr))/2 + sum(p2(2:(Nr-1),rr)) )/5; % finite volume method
    %   prof1(rr) = 2 * solid * ((BLMT_PS + BLMT_SS)/C) / cosd(ang2(rr)) * (0.5 * vel2(rr)^2 / Temp1(rr));  % mixing losses
    %   prof2(rr) = (solid*((BLDT_PS + BLDT_SS)/C + g.tTE) /cosd(ang2(rr)))^2 * (0.5 * vel2(rr)^2 / Temp1(rr)); % profile loss
    % Note - static temps used, but T0 ~ T with low M (~2% error for M=.3)
    % Note - can merge prof1 and prof2...
    %% Base loss - calcuated at radial points, averaged across span (such that casing and hub nodes have half the weighting)
    % Calculate base pressure loss - Dickens (9.11) / Denton Fig 10, Cpb
    % scales with Mach number and tTE/throat
    %   tw = g.tTE*solid/cosd(ang2(rr)); % see Fig 10
    %   Cpb = 0.2/0.15 *tw
    Cpb = -0.1; % Dickens p156 - based on turbines!
    %   base(rr) = - Cpb * solid * g.tTE / cosd(ang2(rr)) * (0.5 * vel2(rr)^2 / Temp1(rr)); % Base pressure loss
    base(:,rr) = - Cpb .* solid(:,rr) .* g.tTE./g.c(:,rr) ./ cosd(adwn) .* (0.5 .* vdwn.^2 ./ Tup(:)); % Base pressure loss
    baseT(rr) = ((base(1,rr) + base(Nr,rr))/2 + sum(base(2:(Nr-1),rr)) )/5; % finite volume method
    %% Tip/Shroud loss - evaluated using casing values
    % Shroud Loss - From Sungho Yoon Paper based on shroud losses for turbines in Denton
    %   hb = rc - rh; % span
    %   mm = (g.gap * d.Cc / hb) * ((tand(ang1(rr)))^2 - (tand(ang2(rr)))^2)^0.5; % Yoon 5 %(g * Cc / hn) * (abs((1/cosd(ang2))^2 - (tand(ang1))^2 ))^0.5 % # mm = (g_s * Cc / H_b) * np.sqrt( (1/np.cos(a2))**2 - (np.tan(a2))**2 ) 
    %   tip(rr) =  2 * mm * (1 - (tand(ang1(rr)) * sind(ang2(rr)) * cosd(ang2(rr))) ) * (0.5 * vel2(rr)^2 / Temp1(rr)); %
    hb = rc - rh; % span
    mm = (g.gap * d.Cc / hb) * ((tand(aup(Nr)))^2 - (tand(adwn(Nr))^2))^0.5; % Yoon 5 %(g * Cc / hn) * (abs((1/cosd(ang2))^2 - (tand(ang1))^2 ))^0.5 % # mm = (g_s * Cc / H_b) * np.sqrt( (1/np.cos(a2))**2 - (np.tan(a2))**2 ) 
    tip(rr) =  2 * mm * (1 - (tand(aup(Nr)) * sind(adwn(Nr)) * cosd(adwn(Nr))) ) * (0.5 * vdwn(Nr)^2 / Tup(Nr)); %
    if rpm2 == 0
        tip(2) = 0; % Stators have no gap, so no shroud loss
    end
    % #tip =  mm * V2**2 * abs( 1 - (np.tan(a1) * np.sin(a2) * np.cos(a2)) ) / T #valid for compressible Yoon 6
    %% Endwall loss - evaluated at hub and casing values !! vel depends on whether the rotor has a rotating casing shroud
    %   % Endwall Loss - comes from cd*V^3 argument but need to think about some more - important mechanisms I think
    %   endwallfun = @(x) (((vLE_SS + (vTE - vLE_SS)*x/Cx)).^3 + ((vLE_PS + (vTE - vLE_PS)*x/Cx)).^3)/2;
    %   % Hub  
    %   endwall1 = Cd * ((rm - hb/2)/rm)^3 / (vel1(rr)*Temp1(rr)*hb*cosd(ang1(rr))) * integral(endwallfun,0,Cx);
    %   % Casing
    %   endwall2 = Cd * ((rm + hb/2)/rm)^3 / (vel1(rr)*Temp1(rr)*hb*cosd(ang1(rr))) * integral(endwallfun,0,Cx);
    % Endwall Loss - comes from cd*V^3 argument but need to think about some more - important mechanisms I think
    endwallfun = @(x,ii) (((vLE_SSe(ii) + (vTEe(ii) - vLE_SSe(ii)).*x./Cx(ii,rr))).^3 + ((vLE_PSe(ii) + (vTEe(ii) - vLE_PSe(ii)).*x./Cx(ii,rr))).^3)./2;
    % Hub  
    endwall_h = Cd * ((rm - hb/2)/rm)^3 / (vup(1)*Tup(1)*hb*cosd(aup(1))) * integral(@(x) endwallfun(x,1),0,Cx(1,rr));
    % Casing
    endwall_c = Cd * ((rm + hb/2)/rm)^3 / (vup(Nr)*Tup(Nr)*hb*cosd(aup(Nr))) * integral(@(x) endwallfun(x,Nr),0,Cx(Nr,rr));
    % Total
    endwall(rr) = (endwall_h + endwall_c);
end

% Sum up losses
prof = prof1 + prof2;
L1calc = prof(1) + baseT(1) + tip(1) + endwall(1); %1st row
L2calc = prof(2) + baseT(2) + tip(2) + endwall(2); %2nd row, note no gap for stator - change this to toggle for CRDF

% Calculating efficiency from entropy
L.eta_s = q.cp*(q.T03s - q.T01) / ( q.cp*(q.T03 - q.T01) + q.T01*(L1calc+L2calc));
if L.eta_s < 0.4
    a=0;
    return
end

%% Change in losses for while loop, Repack structures
% dL1 = (L1calc - L.L1)/L.L1;
% dL2 = (L2calc - L.L2)/L.L2;
% if abs(dL1) > 0.8*abs(L.L1)
%     L.deltaL1 = 0.5*dL1;
%     L.L1 = L.L1 + 0.5*(L1calc-L.L1);
% %     disp('smoothing 1')
% else
%     L.deltaL1 = dL1;
%     L.L1 = L1calc;
% %     disp('not smoothing 1')
% end
% if abs(dL2) > 0.8*abs(L.L2)
%     L.deltaL2 = 0.5*dL2;
%     L.L2 = L.L2 + 0.5*(L2calc-L.L2);
% %     disp('smoothing 2')
% else
%     L.deltaL2 = dL2;
%     L.L2 = L2calc;
% %     disp('not smoothing 1')
% end
L.deltaL1 = (L1calc - L.L1)/L.L1;
L.deltaL2 = (L2calc - L.L2)/L.L2;
L.L1 = L1calc;
L.L2 = L2calc;
L.Loss1.prof=prof(1); L.Loss1.base=baseT(1); L.Loss1.tip=tip(1); L.Loss1.endwall=endwall(1);
L.Loss2.prof=prof(2); L.Loss2.base=baseT(2); L.Loss2.tip=tip(2); L.Loss2.endwall=endwall(2);
L.Loss1.BLT_SS=BLT_SS(:,1); L.Loss1.BLT_PS=BLT_PS(:,1); L.Loss2.BLT_SS=BLT_SS(:,2); L.Loss2.BLT_PS=BLT_PS(:,2); 


%% Calculating operating point quantities - MOVE THIS TO CVanalysis_comp

% Power
d.mdot = q.Mndp01 * g.A1 * p01 / (cp*q.T01)^0.5;
d.P = d.mdot * cp*(q.T03-q.T01);

% Reaction
d.Reaction = (T(:,2)-T(:,1))./(T(:,3)-T(:,1));

% Tip Mach Number
d.Mtip = (rpm1/60*2*pi*rc) / (gam*R*T(Nr,1))^0.5;

% Propulsive efficiency and hover figure of merit
d.Fr = 2*d.u0 / (d.u0+d.u4);
d.Mf = d.Th / d.P * sqrt(d.Th / (2*q.atm.ro*g.A1));

% Fan Efficiency
d.eta = (q.T03s - q.T01) / (q.T03 - q.T01);

%% Repacking into structs
d.U1=U1; d.U2=U2; d.rpm1=rpm1; d.rpm2=rpm2;

a.alpha=alpha; a.v=v; a.vx=vx; a.vth=vth;
a.imid = imid;

q.M = M;
q.T0=T0; q.T=T;
q.p0=p0; q.p=p;
q.ro=ro;
