function [rpm, thr, pwr, fom, re] = exp2plot(data, proporfan, chord)
% Return the thrust and power associated with experimental data of RPM
% CURRENT and PAYLOAD (kg)
% Load relevant variables

load('efficiencies.mat');
load('vars', 'W', 'kv', 'r_cCRIT', 'r_h');
load('APC_FIT.mat', 'APC_FIT');

rhoAIR      = 1.225;                        % kg/m^3
Ax = pi*(r_cCRIT^2 - r_h^2);
%% Get metrics

rpm = data(:,1);

switch proporfan
    case 'prop'
        % Total efficiency
        totaleta = eta.control * eta.motor;
        area = pi*((10*0.0254/2).^2 - (10e-3)^2);
        emptyweight = W.body + W.prop;
        
    case 'fan'
        % Total efficiency
        totaleta = eta.total;
        area = pi*(r_cCRIT.^2 - r_h.^2);
        emptyweight = W.total;
end

% totaleta = 1; 

% Thrust
thr = emptyweight + data(:,3).*9.81/W.numberofprop;

% Power (Electrical)
pwr = 15 .* data(:,2) ./ W.numberofprop;

fom = (thr./pwr) .* sqrt(thr ./ (2 * rhoAIR * area));

vx = sqrt(thr ./ (2*rhoAIR*area));

U = rpm.*(pi/30).*(r_cCRIT + r_h)./2;
v = sqrt(vx.^2 + U.^2);
re = rhoAIR .* v .* chord ./ (18.5e-6);
end