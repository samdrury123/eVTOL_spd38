% Calculates the rotor and stator blade inlet and outlet metal angles

%% Geometry
r.chord         = c_r;     
span            = r_c - r_h;
r.tip_gap       = tip_gap_percent/100 * r_m * (1-HTR)/(1+tip_gap_percent/100*(1-HTR)); % ??? Incidence - due to hub spinning

%% Velocity triangles
% Convert all angles to degrees
alpha1          = rad2deg(alpha1);
beta1           = rad2deg(beta1);
beta2           = rad2deg(beta2);
% Calculate relavent velocity triangle components
tau             = linspace(r_h,r_c,N)';
Vx              = Vx(N_mid);
V2_rel          = Vx./cosd(beta2);
V2_theta_rel    = V2_rel.*sind(beta2);
V2_theta        = tau*omega+V2_theta_rel;
alpha2          = atand(V2_theta/Vx);

%% Rotor inlet/outlet angles
r.chi1          = abs(beta1)-r.I1;


m               = 0.23*(2*ac_r).^2 + abs(beta2)/500;        % m for Carter's rule

r.chi2          = (abs(beta2) - m .* r.chi1 .* sc_r .^ 0.5) ./ (1-m.*sc_r.^0.5);
% r.chi2          = abs(beta2);

%% Stator inlet/outlet angles
s.chord         = c_s;
s.radius        = r.radius;              
s.chi2          = alpha2-s.I2;        % Stator inlet based on absolute flow angle and incidence
s.chi3          = linspace(0,0,N)';     % Stator exit flow is axial

%% Deviation and camber
delta = abs(beta2)-r.chi2;                  % Calculate deviation
theta = r.chi1-r.chi2;                  % Calculate camber

%% Save design parameters
save(file_name{1},'r','s');
