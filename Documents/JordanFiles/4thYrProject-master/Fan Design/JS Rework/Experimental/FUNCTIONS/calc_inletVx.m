function P = calc_inletVx(P)
%% P = calc_inletVx(P)
% Calculate inlet axial velocity (Vx) at casing
% P: struct containing DSA struct
rho = 1.225;
P.Vx = sqrt(2*(-P.DSA.P1)/rho);
end