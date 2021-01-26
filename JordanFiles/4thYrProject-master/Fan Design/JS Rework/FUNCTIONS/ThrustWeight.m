function [omega] = ThrustWeight(phi, psi, rc, rh, mass)
%% Thrust = Weight
% - Equate thrust and weight to determine omega required

rhoAIR = 1.225;
sigma = sqrt(phi.^2 ./ (2*psi));

omega = sqrt( 2*sigma.*(mass.total/4).*9.81 ./ (rhoAIR*pi*(phi.^2).*(rc^4 - rh^4)));