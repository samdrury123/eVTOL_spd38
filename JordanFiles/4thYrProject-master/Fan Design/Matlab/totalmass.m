% Required paylod for fan to be preferred method (lower absolute power use)

% Number of test points
n = 100;

% Densities
rhoAIR      = 1.225;                        % kg/m^3

% Payload
payload = 0;    % kilograms

massrange = linspace(1, 100, n);

for i = 1:n
    disp(i);
    OptimiseMass(massrange(i));
    
    load('vars');

    % Values
    Ax = pi*(r_cCRIT^2 - r_h^2);
    U = rpm*(pi/30)*(r_cCRIT + r_h)/2;
    sigma = sqrt((phi^2)/(2*psi));

    weighttotal = W.total + (payload/W.numberofprop)*9.81;     % Newtons

    %% Thrust and Power
    T = phi*sqrt(2*psi)*rhoAIR*Ax*U^2;

    % Power wrt thrust
    fanpower(i) = (T^(3/2))/sqrt(4*sigma*rhoAIR*Ax);
    
    [proprpm, propweight, propp] = proprpm();
    
    proppower(i) = propp;
    
    if fanpower(i) < proppower(i)
        disp('Limit of Fan Superiority');
        disp(massrange(i));
        break
    end
end

