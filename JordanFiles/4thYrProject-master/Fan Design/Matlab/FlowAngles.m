function [ angles ] = FlowAngles(phi,psi,n,N)
    %Calculates flow angles assuming axial inlet velocity
    index           = ceil(linspace(1,n,N));
    phi             = phi(index);
    psi             = psi(index);
    
    for i = 1:N
        a1(i)           = 0; % axial inlet
%         a1_rel(i)       = -atan(1/phi(i) - tan(a1(i)));
        a1_rel(i)       = atan(-1 / phi(i));
%         a2(i)           = atan(tan(a1(i)) + psi(i) / (phi(i)));
        a2(i)           = atan(psi(i) / phi(i));
%         a2_rel(i)       = -atan(1/phi(i) - tan(a2(i)));
        a2_rel(i)       = atan((psi(i) - 1)/phi(i));
        
        % Assemble flow angle matrix
        angles          = [a1(:), a1_rel(:), a2(:), a2_rel(:)];
    end
end
