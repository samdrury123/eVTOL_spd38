function [phidist, psidist, psi_tsdist, pc] = Distributions(p, phi, psi, r,rm,pp, radius)
%% Return coefficient distributions based on vortex design choice
% 'p' is a string = 'free', 'forced', 'constangle'

switch p
    case 'free'
        % Distributed values
        psidist        = psi.m .* (rm ./ r) .^ 2;
%         psidist(psidist > 1) = 1;
        phidist        = phi.m .* (rm ./ r);
        psi.p = 0;
        pc = 0;
    case 'custom' % To get Jonnies Design
        % Distributed values
        psidist        = psi.m .* (rm ./ r) .^ pp;
        
        psir           = psi.m .* (rm ./ radius) .^ pp; % psi(r)
        dpsidr         = (-pp * psi.m .* rm .^ pp) .* (radius) .^(-pp -1); % d(psi)/dr
        RHS            = -radius.*psir.*dpsidr - 2.*psir.^2; % RHS of radial equilibrium
        
        for i = 1:length(radius)
            del = 1;
            c1 = 0; c2 = 1; c3 = 2.5; 
            while del > 0.000001    
                % Equation to solve
                diff1 = (1-c1).* (phi.m .* rm.^c1 .* radius(i).^(-c1)).^2 - RHS(i);
                diff2 = (1-c2).* (phi.m .* rm.^c2 .* radius(i).^(-c2)).^2 - RHS(i);
                diff3 = (1-c3).* (phi.m .* rm.^c3 .* radius(i).^(-c3)).^2 - RHS(i);
                if diff1 * diff2 > 0
                    c1 = c2;
                    c2 = 0.5*(c2+c3);
                    del = min(abs(diff2), abs(diff3));
                    c = c2;
                elseif diff2*diff3 > 0
                    c3 = c2;
                    c2 = 0.5*(c2+c1);
                    del = min(abs(diff2), abs(diff1));
                    c = c2;
                end
            end
            pc(i) = c; 
        end

        phidistl       = phi.m .* (rm ./ radius) .^ pc; % Solve this with a bisection method. 
        
        if length(r) > 6
            phidist = phidistl;
        else 
            phidist = phidistl(ceil(linspace(1,end, 5)));
            pc = pc(ceil(linspace(1,end, 5)));
        end
        
        psi.p = pp;
        phi.p = pc;
    case 'forced'
        % Distributed values
        psidist        = psi.m .* r .^ 0;
        phidist        = phi.m .* (rm ./ r);
        psi.p = 0;
        phi.p = 1;
        pc = 0;
    case 'constangle'
        % Solution of form phi = phi.m.*r^(-B/A)
        A = 1 + (psi.m/phi.m).^2; B = 1 + 2*(psi.m/phi.m).^2;
        % Distributed
        phidist = phi.m .* (rm ./ r) .^ (B/A);
        psidist = psi.m .* (rm ./ r) .^ (B/A);
        psi.p = B/A
        phi.p = B/A;
        pc = 0;
end

psi_tsdist      = 2.*psidist - phidist.^2;
end