function P = calc_phi(P)
rho = 1.225;
P.omega = P.rpmmean*pi/30;

A = -(0.25^2 / 0.8^2)*(2-P.p) - 1;
B = -P.p;
P.phi_c = sqrt(-P.DSA.P1 ./ (0.5*rho*(P.omega .* P.rc).^2));
P.phi_m = sqrt(-P.DSA.P1 ./ (0.5*rho*(P.omega .* P.rc).^2)) .* (P.rc / sqrt( (P.rc^2 + 20e-3^2)/2 ))^(-A);
P.psi_m = (-P.DSA.P2 ./ (rho*(P.omega .* P.rc).^2) + 0.5.*P.phi_c.^2) .* (P.rc / sqrt( (P.rc^2 + 20e-3^2)/2 ))^(-B);
P.sig_exp = P.phi_m.^2 ./ (2*P.psi_m);
end