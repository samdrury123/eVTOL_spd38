function [S] = calcs(S)
% Return calculated values for specified dataset
%

r_cCRIT = 57.69e-3; r_h = 20e-3; AR = 1.5; rho = 1.225; mu = 18.13e-6;
c = (r_cCRIT - r_h) / AR; Ax = pi*(r_cCRIT^2 - r_h^2);

if strcmp(S.type,'fan')
    S.U = (pi/30) .* S.rpmmean .* 0.5 .* (r_cCRIT + r_h); 
    S.Vx = sqrt(S.T .* S.sig ./ (rho * Ax)); 
    S.V = sqrt(S.U.^2 + S.Vx.^2); 
    S.Re = rho .* S.V .* c ./ mu;
    
    S.phi = sqrt(S.T .* S.sig ./ (rho.*Ax.*S.U.^2));

    S.psi = S.T  ./ (2 * rho .* Ax .* S.sig .* S.U.^2);

    S.sigpredicted = sqrt(S.phi.^2 ./ (2 .* S.psi));
else
    if strcmp(S.name, 'APC')
        S.U = (pi/30) .* S.rpmmean .* 5*25.4e-3; S.Re = rho .* S.U .* S.c ./ mu;
    end
end
end