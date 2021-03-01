function [] = RadialEquilibriumEquation()

%% Inputs 
n           = 999999;
r_h         = 11.6e-3;
r_c         = 35e-3;
r_m         = (r_h + r_c)/2;
s           = r_c - r_h;
r           = linspace(r_h,r_c,n)';
dr          = r(2) - r(1);
omega       = 18000 * (2 * pi / 60);
phi_m         = 0.6;
psi_m       = 0.35;

U           = (r * omega);

n_mid       = (n + 1)/2;

% Set arrays
dh_0        = ones(n,1);
V_theta2    = ones(n,1);
V_x          = ones(n,1);

% Boundary Condition
V_x(n_mid)   = r_m * omega * phi_m;

%% Stage Loading Distribution
% 1 - Free Vortex
% 2 - Forced Vortex
% 3 - Self-generated distribution
flag = 4;

if flag == 1
% Free Vortex
psi         = psi_m * (r_m./r).^2;
elseif flag == 2
% Forced Vortex
psi         = psi_m * ones(n,1);
elseif flag == 3
% Self-generated
a           = ceil(n*0.3);
psi_min     = 0.25;
psi_max     = psi_m;
for i = 1:a
    psi(i)         = psi_min + (i-1)/(a-1) * (psi_m-psi_min);
end
for i = a:n
    psi(i)         = psi_m + (i - a - 1)/(n - a - 1) * (psi_max - psi_m);
end
elseif flag == 4
psi         = psi_m * (r_m ./ r).^1.8;
end

%% Analytical Solution
% Free Vortex
if flag == 1
    K           = psi_m * r_m^2;
    psi_a       = K./r.^2;
    dh_0_a      = psi_a .* U.^2;
    V_theta2_a  = (dh_0_a ./ U);
    V_x_a       = phi_m * r_m * omega * ones(n,1);
    phi_a       = V_x_a ./ U;
elseif flag == 2
    U_m         = r_m * omega;
    K           = psi_m .* r.^2 ;
    psi_a       = K./ r.^2;
    dh_0_a      = psi_a .* U.^2;
    V_theta2_a  = (dh_0_a ./ U);
    k1          = psi_m * omega;
    k2          = 2 * r_m^2 * k1 * (omega - k1) - (phi_m * U_m)^2;
    V_x_a       = sqrt(2 * r.^2 * k1 * (omega - k1) - k2);
    phi_a       = V_x_a ./ U;
end
%% Solver

for i = 1:n
    dh_0(i)         = (r(i) * omega)^2 * psi(i);
    V_theta2(i)     = (r(i) * omega) * psi(i);
end

for i = n_mid:(n-1)  
    V_x(i+1)     = V_x(i) + (dh_0(i+1) - dh_0(i))/V_x(i) - (V_theta2(i)/V_x(i)) * (V_theta2(i+1) - V_theta2(i)) - ((V_theta2(i))^2 / V_x(i)) * (dr/r(i));
end

for i = n_mid:-1:2   
    V_x(i-1)     = -(-V_x(i) + (dh_0(i) - dh_0(i-1))/V_x(i) - (V_theta2(i)/V_x(i)) * (V_theta2(i) - V_theta2(i-1)) - ((V_theta2(i))^2 / V_x(i)) * (dr/r(i)));
end

%% Flow Coefficient

phi         = V_x ./ U;

%% Plot
figure(1)
subplot(2,2,1)
hold on
grid on
plot(V_x,(r-r_h)/s)
if flag == 1 | flag == 2
plot(V_x_a,(r-r_h)/s)
end
xlabel('V_x')
ylabel('(r-r_h)/s')
legend('Numerical Solution','Analytical Solution')
hold off

figure(1)
subplot(2,2,2)
hold on
grid on
plot(phi,(r-r_h)/s)
if flag == 1 or flag == 2
plot(phi_a,(r-r_h)/s)
end
xlabel('\phi')
ylabel('(r-r_h)/s')
legend('Numerical Solution','Analytical Solution')
hold off

figure(1)
subplot(2,2,3)
hold on
grid on
plot(dh_0,(r-r_h)/s)
if flag == 1 or flag == 2
plot(dh_0_a,(r-r_h)/s)
end
xlabel('\Deltah_0')
ylabel('(r-r_h)/s')
legend('Numerical Solution','Analytical Solution')
hold off

figure(1)
subplot(2,2,4)
hold on
grid on
plot(psi,(r-r_h)/s)
if flag == 1 or flag == 2
plot(psi_a,(r-r_h)/s)
end
xlabel('\psi')
ylabel('(r-r_h)/s')
legend('Numerical Solution','Analytical Solution')
hold off
