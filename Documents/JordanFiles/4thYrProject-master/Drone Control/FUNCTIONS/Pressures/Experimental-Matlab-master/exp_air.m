% Generate air property data

% Perfect gas constants
air.cp = 1003.308; air.ga = 1.40047;
air.R = 287.058;
air.cv = air.cp / air.ga; air.R = air.cp - air.cv;

% Viscosity constants for Sutherland's Formula
air.mu_ref = 18.27e-6; air.T_ref = 291.15; air.C = 120;

% Water vapour constants for ideal gas mixture equations
vap.R = 461.495;