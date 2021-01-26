function [V_x,phi,psi] = RadialEquilibriumFunction(r_h,r_c,omega,phi_mid,psi_mid,n,p)

    global flagplot

    %% Geometry
    r_m         = (r_h+r_c)/2;            % Mean radius (m)
    span        = r_c-r_h;                % Blade span (m)
    r           = linspace(r_h,r_c,n)';
    
    %% Inital calculations
    dr          = r(2)-r(1);
%     omega       = rpm*(2*pi/60);
    U           = (r*omega);
    n_mid       = (n+1)/2;
    
    % Set arrays
    dh_0        = ones(n,1);
    V_theta2    = ones(n,1);
    V_x         = ones(n,1);
    
    % Boundary Condition
    V_x(n_mid)  = r_m*omega*phi_mid;

    %% Stage Loading Distribution
    % p = 2 - Free Vortex
    % p = 0 - Forced Vortex
    % Stage loading distribtion from radial equilibrium equation
    psi         = psi_mid*(r_m./r).^p;

    %% Solver
    for i = 1:n
        % Use stage loading coefficient to obtain stagnation enthalpy
        % change
        dh_0(i)         = (r(i)*omega)^2*psi(i);
        V_theta2(i)     = (r(i)*omega)*psi(i);
    end
    for i = n_mid:(n-1)
        % Forward Euler method
        V_x(i+1)        = V_x(i)+(dh_0(i+1)-dh_0(i))/V_x(i)-(V_theta2(i)/V_x(i))*(V_theta2(i+1)-V_theta2(i))-((V_theta2(i))^2/V_x(i))*(dr/r(i));
    end
    for i = n_mid:-1:2   
        % Backward Euler method
        V_x(i-1)        = -(-V_x(i)+(dh_0(i)-dh_0(i-1))/V_x(i)-(V_theta2(i)/V_x(i))*(V_theta2(i)-V_theta2(i-1))-((V_theta2(i))^2/V_x(i))*(dr/r(i)));
    end

    %% Flow Coefficient
    phi         = V_x./U;

    %% Plot
    if flagplot == 1
        % Plot variation in Phi across the span
        figure(4)
        subplot(1,2,1)
        grid on
        hold on
        plot((r-r_h)/span,phi)
        hold off
        ylabel('\phi')
        xlabel('(r-r_h)/s')
        % Plot variation in Psi across the span
        figure(4)
        subplot(1,2,2)
        grid on
        hold on
        plot((r-r_h)/span,psi)
        hold off
        ylabel('\psi')
        xlabel('(r-r_h)/s')
        hold off
    end
end