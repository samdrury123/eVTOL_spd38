%% ---- FUNCTION: VELOCITY TRIANGLES ---- %%
function [a,q] = VelTriangles(d,q,L,CRDF)

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
u1=d.u1; U1=d.U1; U2=d.U2; Um=d.Um; psi=d.psi; sigma=d.sigma;
gam=q.gam; R=q.R; cp=q.cp;
p01=q.p01; p4=q.p4; T03=q.T03; T1=q.T1; M1=q.M1; Mndp4=q.Mndp4;

% Set up Mach Number vector for lookups
M = 0.001:0.001:1;
Mndp0 = gam / (gam-1)^0.5 .* M .* (1 + (gam-1)/2.*M.^2).^(-0.5*(gam+1)/(gam-1));

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
vt2 = Um * psi; %%%%CHECK THISSSS
vt2r(1) = vt2 - U1;
vt2r(2) = vt2 - U2;
p02r(1)=p01r*exp(-L.Lr/R);
T02r=T01r;
if L.deltaLr == 1
    vx2 = vx1; %this is a first guess
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
vx2_low = (1-bisect_fac)*vx1; vx2_high = (1+bisect_fac)*vx1;
deltavx = 1;

while abs(deltavx) > 0.0001/100 %0.00001
    vx2 = 0.5*(vx2_low + vx2_high);
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
end

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


if CRDF == 0 % Downstream Stator
    v2r(2)=v2; 
    p03=p02*exp(-L.Ls/R); p03r=p03;
    a3 = 0; a3r=0;     % Change this to allow non-zero?
    vt3 = 0; vt3r=0;
    Mndp03 = Mndp4 * (p4 / p03) * sigma;
    M3 = interp1(Mndp0,M,Mndp03); M3r=M3;
    T3 = T03 .* (1 + (gam-1)/2.*M3.^2).^-1;
    v3 = M3 * (gam*R*T3)^0.5;
    vx3 = v3; v3r = v3;
    p3 = p03*(1+0.5*(gam-1)*M3^2)^(-gam/(gam-1));
    ro3 = p3/(R*T3);
else % Downstream Rotor
    p03r=p02r(2)*exp(-L.Ls/R);
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
    c3 = (gam*R*T3)^0.5;
    M3 = v3/c3;
    p3 = p03r*(1+0.5*(gam-1)*M3r^2)^(-gam/(gam-1));
    p03 = p3*(1+(gam-1)*0.5*M3^2)^(gam/(gam-1));
    ro3 = p3/(R*T3);
end

%% Repacking into structs
a.v1=v1; a.v2=v2; a.v3=v3; a.vx1=vx1; a.vx2=vx2; a.vx3=vx3; a.v1r=v1r; a.v2r=v2r; a.v3r=v3r;
a.vt1=vt1; a.vt2=vt2; a.vt3=vt3; a.vt1r=vt1r; a.vt2r=vt2r; a.vt3r=vt3r;
a.a1=a1; a.a2=a2; a.a3=a3; a.a1r=a1r; a.a2r=a2r; a.a3r=a3r;

q.M1r=M1r; q.M2=M2; q.M3=M3; q.M2r=M2r; q.M3r=M3r;
q.T01r=T01r; q.T02r=T02r; q.T03r=T03r; q.T02=T02; q.T2=T2; q.T3=T3;
q.p01r=p01r; q.p02r=p02r; q.p03r=p03r; q.p02=p02; q.p03=p03; q.p1=p1; q.p2=p2; q.p3=p3;
q.ro1=ro1; q.ro2=ro2; q.ro3=ro3;
