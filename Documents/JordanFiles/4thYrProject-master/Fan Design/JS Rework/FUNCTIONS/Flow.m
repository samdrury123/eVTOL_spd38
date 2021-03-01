function [angle, V, reaction] = Flow(r, omega, phi, psi)
global crdfflag
% Return velocities
V.x         = (r .* omega) .* phi;
V.theta2    = (r .* omega) .* psi;

V.rel1      = sqrt((r.*omega) .^ 2 + V.x .^ 2);

% Return flow angles
if crdfflag == 1
    V.rel1      = sqrt((r.*(omega/2)) .^ 2 + V.x .^ 2);
    angle.a1    = atand(1./(2*phi));
    angle.a3    = atand(1./(2*phi));
    
    angle.b1    = atand(-1 ./ (2*phi));
    angle.a2    = atand((2*psi + 1) ./ (2*phi) );
    angle.b2    = atand(1./(2*phi) - tand(angle.a2));

    angle.b3    = angle.a3;
    angle.i1    = 0;
    angle.i3    = 0;
else 
    V.rel1      = sqrt((r.*omega) .^ 2 + V.x .^ 2);
    angle.a1    = 0;
    angle.a3    = 0;
    
    angle.b1    = atand(-omega .* r ./ V.x);
    angle.a2    = atand(V.theta2 ./ V.x);
    angle.b2    = atand((V.theta2 - r .* omega) ./ V.x);

    angle.b3    = 0;
    angle.i1    = 0;
    angle.i3    = 0;
end

% Other velocities
V.abs1      = V.x ./ cosd(angle.a1);
V.abs2      = V.x ./ cosd(angle.a2);
V.rel2      = V.x ./ cosd(angle.b2);
V.rel3      = V.x ./ cosd(angle.b3);
% Reaction
reaction = 1-0.5*(phi.^2./psi).*(tand(angle.a2).^2);
end