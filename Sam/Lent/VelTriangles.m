%% ---- FUNCTION: VELOCITY TRIANGLES ---- %%
function [d,g,a,q] = VelTriangles(d,g,q,L,k)

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

%% Unpack
u1=d.u1; Um=d.Um; psi=d.psi; sigma=d.sigma;
gam=q.gam; R=q.R; cp=q.cp;
p01=q.p01; p4=q.p4; T03=q.T03; T1=q.T1; M1=q.M1; Mndp4=q.Mndp4;

% Set up Mach Number vector for lookups
Mach = 0.001:0.001:1;
Mndp0 = gam / (gam-1)^0.5 .* Mach .* (1 + (gam-1)/2.*Mach.^2).^(-0.5*(gam+1)/(gam-1));

% Set up radial positions, trying to make rm close to a point in r
Nr = 6; % Number of radial slices
r = linspace(g.rh,g.rc,Nr)';
[~,imid] = min(abs((r-g.rm)/(g.rc-g.rh)));
r(imid) = g.rm;
g.r=r;
% thresh=0.1; rmin=1;
% while rmin>thresh
%     Nr=Nr+1;
%     r = linspace(g.rh,g.rc,Nr)';
%     [rmin,i] = min(abs((r-g.rm)/(g.rc-g.rh)));
% end

% Split blade speeds between rows (CRDF)
U1 = Um*k;      rpm1 = U1/(2*pi*g.rm) * 60;
U2 = -Um*(1-k); rpm2 = U2/(2*pi*g.rm) * 60;

% Set up spanwise stage loading coefficient, enthalpy and vtheta
psi = psi .* (g.rm ./ r).^ (g.n+1);

% Set Vx constant across span at inlet
vx(1:Nr,1) = ones(size(r)) * u1;

% Calculate velocities and angles upstream of blade 1
vth(1:Nr,1,1) = 0; % alpha1 = 0
vth(:,1,2) = vth(:,1,1) - U1;
vth(:,1,3) = vth(:,1,1) - U2;

% Calculate velocities and angles downstream of blade 1, including axial
% velocity distribution from radial equilibrium equation
vth(:,2,1) = U1 * psi / k; %(Um*psi)
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
p02r(1)=p01r*exp(-L.Lr/R);
T02r=T01r;
if L.deltaLr == 1
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

while abs(deltavx2) > 0.0001/100 %0.00001
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
T02 = T2 + 0.5*v(imid,2,1)^2/cp;
c2 = (gam*R*T2)^0.5;
M(:,2,1) = v(:,2,1)./c2;
p2 = p02r(1)*(1+(gam-1)*0.5*M(imid,2,2)^2)^(-gam/(gam-1));
p02 = p2*(1+(gam-1)*0.5*M(imid,2,1)^2)^(gam/(gam-1));
ro2 = p2/(R*T2);

%% Radial Equilibrium between blade rows (Position 2)

% Use Euler work equation to calculate change in enthalpy
dh0 = U1 .* vth(:,2,1);

% Numerically solve Simple Radial Equilibrium equation

% Solve
for ii = imid:(Nr-1)
    dr = r(ii+1) - r(ii); % Radial increment
    vx(ii+1,2,1)     = vx(ii,2,1) + (dh0(ii+1) - dh0(ii))/vx(ii,2,1) - (vth(ii,2,1)/vx(ii,2,1)) * (vth(ii+1,2,1) - vth(ii,2,1)) - ((vth(ii,2,1)).*(vth(ii+1,2,1)) / vx(ii,2,1)) * (dr/r(ii));
end
for ii = imid:-1:2
    dr = r(ii) - r(ii-1);
    vx(ii-1,2,1)     = -(-vx(ii,2,1) + (dh0(ii) - dh0(ii-1))/vx(ii,2,1) - (vth(ii,2,1)/vx(ii,2,1)) * (vth(ii,2,1) - vth(ii-1,2,1)) - ((vth(ii,2,1)).*(vth(ii-1,2,1)) / vx(ii,2,1)) * (dr/r(ii)));
end
vx(:,2,2) = vx(:,2,1); vx(:,2,3) = vx(:,2,1);

%% Solve for vx3
alpha(:,2,3) = atand(vth(:,2,3)./vx(:,2,1));
v(:,2,3) = (vx(:,2,1).^2 + vth(:,2,3).^2).^0.5;
M(:,2,3) = v(:,2,3)/c2;
% a2r(2) = atand(vt2r(2) / vx2);
% v2r(2) = (vx2^2 + vt2r(2)^2)^0.5;
% M2r(2) = v2r(2)/c2;
p02r(2) = p2*(1+(gam-1)*0.5*M(imid,2,3)^2)^(gam/(gam-1));
Mndp02r(2) = gam / (gam-1)^0.5 .* M(imid,2,3) .* (1 + (gam-1)/2.*M(imid,2,3).^2).^(-0.5*(gam+1)/(gam-1));

alpha(:,3,1) = 0; % Zero exit swirl

if k == 1 % Downstream Stator
    p03=p02*exp(-L.Ls/R); p03r=p03;
    T03r=T03;
    alpha(:,3,3) = 0;
    Mndp03 = Mndp4 * (p4 / p03) * sigma;
    M(:,3,1) = interp1(Mndp0,Mach,Mndp03);
    M(:,3,3) = M(:,3,1);
    T3 = T03 .* (1 + (gam-1)/2.*M(imid,3,1).^2).^-1;
    v(:,3,1) = M(:,3,1) * (gam*R*T3)^0.5;
    v(:,3,3) = v(:,3,1);
    vx(:,3,1) = v(:,3,1); vx(:,3,2) = vx(:,3,1); vx(:,3,3) = vx(:,3,1);
    p3 = p03*(1+0.5*(gam-1)*M(imid,3,1)^2)^(-gam/(gam-1));
    ro3 = p3/(R*T3);
else % Downstream Rotor
    p03r=p02r(2)*exp(-L.Ls/R);
    T03r=T02r;
    
    deltavx3 = 1;
    vx3_low = (1-bisect_fac)*vx(:,1,1); vx3_high = (1+bisect_fac)*vx(:,1,1);
    while abs(deltavx3) > 0.0001/100 %0.00001
        vx(:,3,1) = 0.5*(vx3_low + vx3_high);
        v(:,3,1) = vx(:,3,1);
        v(:,3,3) = (vx(:,3,1).^2 + vth(:,3,3).^2).^0.5;
        alpha(:,3,3) = atand(vth(:,3,3)./vx(:,3,1));
%         v3r = (vx3^2 + vt3r^2)^0.5;
%         a3r = atand(vt3r / vx3);
        Mndp03r = Mndp02r(2) * p02r(2)/p03r * cosd(alpha(imid,2,3))/cosd(alpha(imid,3,3));
        M(:,3,3) = interp1(Mndp0,Mach,Mndp03r)* v(:,3,3)./v(imid,3,3);
        T3 = T03r - 0.5*v(imid,3,3)^2 / cp;
        vx3calc = (M(:,3,3).^2*gam*R*T3 - vth(:,3,3).^2).^0.5;
        if (vx3calc(imid) - vx(imid,3,1)) > 0
            vx3_low = vx(:,3,1);
        else
            vx3_high = vx(:,3,1);
        end
        deltavx3 = (vx3calc(imid) - vx(imid,3,1))/vx(imid,3,1);
    end
    vx(:,3,1) = v(:,3,1); vx(:,3,2) = vx(:,3,1); vx(:,3,3) = vx(:,3,1);
    c3 = (gam*R*T3)^0.5;
    M(:,3,1) = v(:,3,1)/c3;
    p3 = p03r*(1+0.5*(gam-1)*M(imid,3,3)^2)^(-gam/(gam-1));
    p03 = p3*(1+(gam-1)*0.5*M(imid,3,1)^2)^(gam/(gam-1));
    ro3 = p3/(R*T3);
end

% Backfilling arrays
alpha(:,3,2) = atand(vth(:,3,2) ./ vx(:,3,2));
v(:,3,2) = vx(:,3,2) ./ cosd(alpha(:,3,2));

%% Calculating operating point quantities

% Power
d.mdot = q.Mndp01 * g.A1 * p01 / (cp*q.T01)^0.5;
d.P = d.mdot * cp*(q.T03-q.T01);

% Reaction
d.Reaction = (T2-q.T1)/(T3-q.T1);

% Tip Mach Number
d.Mtip = (rpm1/60*2*pi*g.rc) / (gam*R*q.T1)^0.5;

% Propulsive efficiency and hover figure of merit
d.Fr = 2*d.u0 / (d.u0+d.u4);
d.Mf = d.Th / d.P * sqrt(d.Th / (2*q.atm.ro*g.A1));

% Fan Efficiency
d.eta = (q.T03s - q.T01) / (q.T03 - q.T01);

%% Repacking into structs
d.U1=U1; d.U2=U2; d.rpm1=rpm1; d.rpm2=rpm2;

a.alpha=alpha; a.v=v; a.vx=vx; a.vth=vth;

q.M = M;
q.T01r=T01r; q.T02r=T02r; q.T03r=T03r; q.T02=T02; q.T2=T2; q.T3=T3;
q.p01r=p01r; q.p02r=p02r; q.p03r=p03r; q.p02=p02; q.p03=p03; q.p1=p1; q.p2=p2; q.p3=p3;
q.ro1=ro1; q.ro2=ro2; q.ro3=ro3;
