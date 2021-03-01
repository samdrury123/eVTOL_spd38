function P = calc_pressurerise(P)
%% P = calc_pressurerise(P)
% Calculate the pressure rise across various components
% P: struct containing a DSA struc

% Static Pressure rise across...
% Rotor
P.delP_rotor = P.DSA.P12 - P.DSA.P1;
% Stator
P.delP_stator = P.DSA.P2 - P.DSA.P12;
% Stage
P.delP_stage = P.DSA.P2 - P.DSA.P1;

% Diffuser
P.delP_diffuser = P.DSA.PDE - P.DSA.PDI;
P.Cpr = (P.DSA.PDE - P.DSA.P2) ./ (0.5*1.225*P.Vx^2);
P.diff_in_loss = (P.DSA.PDI - P.DSA.P2) ./ (0.5*1.225*P.Vx^2);
end