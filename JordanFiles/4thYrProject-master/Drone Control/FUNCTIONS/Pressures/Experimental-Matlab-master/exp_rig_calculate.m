function [e,V] = exp_rig_calculate(e,r_mid)
% EXP_RIG_CALCULATE  Calculate operating operating point of a rig with read data
%
%   e = EXP_RIG_CALCULATE(e,r_mid)
%
%   e - input data structure with rig instrument readings
%   r_mid - optional midheight radius for flow coefficient


%% Process readings into useful data

% Recall air properties
exp_air

% Use rig temperature if atmospheric not logged
if isfield(e,'Ta') == 0
    e.Ta = e.T(e.N.To);
end

% Extract shaft speed and atmospheric readings
w = e.w; Pa = e.Pa; Ta = e.Ta; Ha = e.Ha;

% Calculate gas properties of humid air using AMJD / KRE approach
[R,cp,ga] = humidityCorrection(Pa,Ta,Ha);

% Calculate reference flow conditions at the rig inlet
P = mean(e.P(:,e.N.P),2) - e.P(:,e.N.Pa) + Pa;
Po = mean(e.P(:,e.N.Po),2) - e.P(:,e.N.Pa) + Pa;
To = mean(e.T(:,e.N.To),2);

% Inlet velocity from compressible flow relations, ignore swirl or bellmouth contraction
M = ( ((Po ./ P).^((ga-1)./ga) - 1) .* (2./(ga - 1)) ).^0.5;
V = M .* (ga .* R .* To).^0.5;

% Calculate static temperature from compressible flow relation
T = To ./ (1 + M.^2 .* (ga -1) / 2);

% Density from ideal gas law
e.ro = P ./ (R .* T);

% Calculate midheight blade speed
e.U_mid = w * r_mid;


%% Calculate operating point
    
% Calculate flow coefficient
e.phi = V ./ e.U_mid;
    
% Calculate non-dimensional flow and speed functions
e.mrtp = e.ro .* V .* (cp .* To).^0.5 ./ Po;
e.nrt = w ./ (cp .* To).^0.5;

% Take only real parts to improve robustness for zero rig speed
e.phi = real(e.phi); e.mrtp = real(e.mrtp); e.nrt = real(e.nrt);

% Record gas properties
e.R = R; e.cp = cp; e.ga = ga;


end