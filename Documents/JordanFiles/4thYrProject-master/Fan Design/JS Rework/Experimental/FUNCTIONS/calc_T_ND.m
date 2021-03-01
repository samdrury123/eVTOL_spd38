function P = calc_T_ND(P, area)
rho = 1.225;
U = (P.rpmmean*pi/30)* sqrt((P.rc + 20e-3^2)/2);
P.T_ND = P.T ./ (0.5*rho*U.^2 .* area);
end