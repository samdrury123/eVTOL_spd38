%% Input Parameters
global flagplot

%% Radial Equilibrium
n               = 99999;
N_mid           = ceil(N/2);
[V_x,phi,psi]   = RadialEquilibriumFunction(r_h,r_c,omega,phi_mid,psi_mid,n,p);

%% Thrust and Power
% Mean-line blade speed
U_m             = r_m*omega;
% Change in stagnation enthalpy from stage loading
r_temp          = linspace(r_h,r_c,n)';
dr              = r_temp(2)-r_temp(1);
U               = r_temp.*omega;
dh_0            = psi.*U.^2;
% Reset values
mdot            = 0;
dH_0            = 0;
% Integrate for total mass flow rate and stagnation enthalpy change
for i = 1:n
    mdot        = mdot+2*pi*rho*V_x(i)*r_temp(i)*dr;
    dH_0        = dH_0+dh_0(i)*rho*V_x(i)*2*pi*r_temp(i)*dr;
end
% Reset values
phi_m           = 0;
psi_m           = 0;

% Integrate for mass averaged phi and psi
for i = 1:n
    psi_m       = psi_m+(1/mdot)*dh_0(i)/U(i)^2*rho*V_x(i)*2*pi*r_temp(i)*dr;
    phi_m       = phi_m+(1/mdot)*rho*V_x(i)^2/U(i)*2*pi*r_temp(i)*dr;        
end
% Inlet conditions
V0              = mdot/(rho*A0); % uniform inlet flow
P00             = Pinf+0.5*rho*Vinf^2;
% Mass averaged specific stagnation enthalpy change
dh_0_m          = dH_0/mdot; % mass averaged stagnation enthalpy change
% Mass averaged stagnation pressure change
dP_0_m          = rho*dh_0_m*eta.flow;
% Stagnation pressures
P01             = eta.inlet*P00;
P02             = P01+dP_0_m;
% Calculate nozzle ratio
if flag == 1
    % Remaining stagnation pressures
    P03             = eta.stator*P02;
    P04             = eta.nozzle*P03;
    % Jet presure and velocity
    P4              = Pinf;
    V4              = abs((2*(P04-P4)/rho))^0.5;
    % Absolute inlet/outlet velocities
    V1              = mdot/(rho*A1);    % assume no swirl at inlet to rotor
    V3              = mdot/(rho*A3);    % assume no swirl at exit to stator
    % Static pressures
    P3              = P03-0.5*rho*V3^2;
    P1              = P01-0.5*rho*V1^2;
    P0              = P00-0.5*rho*V0^2;
    % Nozzle Ratio
    NR              = V4/V3;
    % Exit Area (from nozzle ratio)
    A4              = A3/NR;
end
% Fixed Design (NR specified)
if flag == 2
    % Absolute velocities
    P4              = Pinf;
    V1              = mdot/(rho*A1); % assume no swirl at inlet to rotor
    V3              = mdot/(rho*A3); % assume no swirl at exit to stator
    V4              = V3/NR; % Set by continuity
    % Stagnation pressures
    P04             = P4+0.5*rho*V4^2;
    P03             = P04/eta.nozzle;
    P02             = P03/eta.stator;
    dP_0_m          = P02-P01;
    % Calculate isetropic propulsor efficiency
    eta.flow        = dP_0_m/(rho*dh_0_m);
    % Static pressures
    P3              = P03-0.5*rho*V3^2;
    P1              = P01-0.5*rho*V1^2;
    P0              = P00-0.5*rho*V0^2;
    % Exit Area
    A4              = A3*NR;
end

%% Thrust and Power
% Calculate thrust
Thrust          = mdot*(V4);
% Calculate shaft power and motor power
Pow_s           = dH_0;                     % Shaft power
Pow_in          = dH_0/eta.motor;           % Electrical power in

%% Flow Angles
r_range_b       = linspace(r_h,r_c,N);
[angles]        = FlowAngles(phi,psi,n,N);

%% Velocity Triangles
% Extract flow angles
alpha1          = angles(:,1);
beta1           = angles(:,2); % relative
alpha2          = angles(:,3);
beta2           = angles(:,4); % relative
Vx              = V_x(ceil(linspace(1,n,N)));
% Rotor inlet
V1_abs          = Vx ./ cos(alpha1);
V1_theta        = Vx .* tan(alpha1);
V1_rel          = Vx ./ cos(beta1);
V1_theta_rel    = Vx .* tan(beta1);

% Rotor outlet
V2_abs          = Vx./cos(alpha2);
V2_theta        = Vx.*tan(alpha2);
V2_rel          = Vx./cos(beta2);
V2_theta_rel    = Vx.*tan(beta2);
V2              = V2_abs(N_mid);
P2              = P02-0.5*rho*V2^2;

%% Reaction
R               = 1-0.5*(phi(ceil(n/2))^2/psi(ceil(n/2)))*(tan(alpha2(N_mid))^2)

%% Blade Count
[r.blades,s.blades,sc_r,sc_s,c_r,c_s]  = BladeNumber(AR_rotor,AR_stator,DF_r,DF_s,r_h,r_c,angles,N);

%% Plot Velocity Triangles
C_y             = 1.5*Vx(1);
j_range         = ceil(linspace(1,N,3));

for i = 1:3
    j                   = j_range(i);
    C_x                 = 2*(i-1)*U_m;
    Rotor_X(i,:)        = [0,V1_theta(j),0,V1_theta_rel(j),0] + C_x;
    Rotor_Y(i,:)        = [0,-Vx(j),-Vx(j),-Vx(j),0];
    Stator_X(i,:)       = [0,V2_theta(j),0,V2_theta_rel(j),0,0] + C_x;
    Stator_Y(i,:)       = [0,-Vx(j),-Vx(j),-Vx(j),0,-Vx(j)] - 1.5*C_y;
end
% Produce plots if plotting flag is set
if flagplot == 1
    figure(2)
    hold on
    axis equal
    axis ([-50 250 -100 20 ])
    for i = 1:3
        plot(Rotor_X(i,:),Rotor_Y(i,:))
    end
    for i = 1:3
        plot(Stator_X(i,:),Stator_Y(i,:))
    end
    legend('Hub','Mid Span','Casing')
    hold off
end
% Sections labelled + 1 from above
Section         = [1,2,3,4,5];
% Pre-plot calculations
P               = [P0,P1,P2,P3,P4];
P_0             = [P00,P01,P02,P03,P04];
V               = [V0,V1,V2,V3,V4];
qinf            = 0.5*rho*U_m^2;
CP              = (P-Pinf)/qinf;
CP0             = (P_0-Pinf)/qinf;
% Produce plots if plotting flag is set
if flagplot == 1
    % Plot static pressure coefficient across section
    figure(1)
    subplot(1,2,1)
    hold on
    grid on
    plot(Section,CP)
    xlabel('Section')
    ylabel('C_P')
    hold off
    % Plot stagnation pressure coefficient across section
    figure(1)
    subplot(1,2,2)
    hold on
    grid on
    plot(Section,CP0)
    xlabel('Section')
    ylabel('C_{P0}')
    hold off
end