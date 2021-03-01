function [V, angle, R, S, carter, delta] = Deviation(TYPE, V, angle, phi, psi, R, S, rc, rh, rm)
%% Calculate Deviation (type = 'airfoil' or 'blade')
% Deviation estimated using Carter's rule

%% Inlet metal angle
switch TYPE
    case 'blade'        
        % Carter's Rule
        [angle.m, V.m, R.m, S.m, delta.m, carter.m] = Carter(angle.m, V.m, rc, rh, rm); % Midline
        [angle.sec, V.sec, R.sec, S.sec, delta.sec, carter.sec] = Carter(angle.sec, V.sec, rc, rh, rm); % Section
        [angle.span, V.span, R.span, S.span, delta.span, carter.span] = Carter(angle.span, V.span, rc, rh, rm); % Span
        
    case 'airfoil'
        %% Fan aerofoil assumption
        % psi_{optimum} = (phi/sqrt(2))*(l/s)*C_L
        % C_L = 2*pi*sin(theta/2)
        % C_L approx = pi*theta
        
        %% Angles, no delta apply Runge-Kutta condition
        % Midspan
        angle.m.chi1 = angle.m.b1+angle.i1;
        angle.m.chi2 = angle.m.b2;
        angle.m.chi3 = angle.m.a2+angle.i3;
        angle.m.chi4 = 0;
        R.m.sc = (phi.m ./ (sqrt(2).*psi.m)).*pi.*(angle.m.b2 - angle.m.b1);
        
        % Distributed values
        angle.span.chi1 = angle.span.b1+angle.i1;
        angle.span.chi2 = angle.span.b2;
        angle.span.chi3 = angle.span.a2+angle.i3;
        angle.span.chi4 = 0;
        R.span.sc = (phi.span ./ (sqrt(2).*psi.span)).*pi.*(angle.span.b2 - angle.span.b1);
        
        % Section values
        angle.sec.chi1 = angle.sec.b1+angle.i1;
        angle.sec.chi2 = angle.sec.b2;
        angle.sec.chi3 = angle.sec.a2+angle.i3;
        angle.sec.chi4 = 0;
        R.sec.sc = (phi.sec ./ (sqrt(2).*psi.sec)).*pi.*(angle.sec.b2 - angle.sec.b1);
        
%         % Calculate DF
%         DF.span = (1 - V.span.rel2 ./ V.span.rel1) + (V.span.theta2 ./ (2.*V.span.rel1)) .* R.span.sc;
%         DF.sec = (1 - V.sec.rel2 ./ V.sec.rel1) + (V.sec.theta2 ./ (2.*V.sec.rel1)) .* R.sec.sc;
        
        m.span = 0.23 + abs(angle.span.b2) ./ 500;
        angle.span.chi2 = (angle.span.b2 - angle.span.chi1 .* m.span .* R.span.sc .* (0.5)) ./ (1 - m.span.*R.span.sc.*(0.5));
        
        delta = angle.span.chi2 - angle.span.b2;

end
end