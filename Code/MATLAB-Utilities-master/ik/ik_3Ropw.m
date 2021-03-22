function theta123 = ik_3Ropw( robot_parameter, C )
%theta = ik_3Ropw( robot_parameter, C )
%   Computes the inverse kinematics of serial manipulators with an
%   ortho-parallel basis.
%
%   robot_parameter ... [a1, a2, b1, c1, c2, c3, c4]
%   C ... [cx0, cy0, cz0]
%   theta is a matrix with all four possible solutions.
%   (c) Mathias Brandstötter | IACE | UMIT | 2014 | Version 1.4

a_1 = robot_parameter(1);
a_2 = robot_parameter(2);
b_1 = robot_parameter(3);
c_1 = robot_parameter(4);
c_2 = robot_parameter(5);
c_3 = robot_parameter(6);
% c4 is not used in this function

cx0 = C(1);
cy0 = C(2);
cz0 = C(3);

nx1 = sqrt(cx0^2+cy0^2-b_1^2)-a_1;
s12 = nx1^2+(cz0-c_1)^2; % s12 = s1^2
s22 = (nx1+2*a_1)^2+(cz0-c_1)^2; % s22 = s2^2
k = sqrt(a_2^2+c_3^2);

theta123 = ones(3,4)*1i; % alloction

if isreal(nx1)
    if isreal(nx1+a_1)
        theta123(1,1) = atan2(cy0,cx0)-atan2(b_1,nx1+a_1);
        theta123(1,2) = theta123(1,1);
        theta123(1,3) = atan2(cy0,cx0)+atan2(b_1,nx1+a_1)-pi;
        theta123(1,4) = theta123(1,3);
    end % if

        theta123(2,1) = -acos((s12+c_2^2-k^2)/(2*s12^0.5*c_2))+atan2(nx1,cz0-c_1);
        theta123(2,2) =  acos((s12+c_2^2-k^2)/(2*s12^0.5*c_2))+atan2(nx1,cz0-c_1);
    if isreal(nx1+2*a_1)
        theta123(2,3) = -acos((s22+c_2^2-k^2)/(2*s22^0.5*c_2))-atan2(nx1+2*a_1,cz0-c_1);
        theta123(2,4) =  acos((s22+c_2^2-k^2)/(2*s22^0.5*c_2))-atan2(nx1+2*a_1,cz0-c_1);
    end % if

    theta123(3,1) =  acos((s12-c_2^2-k^2)/(2*c_2*k))-atan2(a_2,c_3);
    theta123(3,2) = -acos((s12-c_2^2-k^2)/(2*c_2*k))-atan2(a_2,c_3);
    theta123(3,3) =  acos((s22-c_2^2-k^2)/(2*c_2*k))-atan2(a_2,c_3);
    theta123(3,4) = -acos((s22-c_2^2-k^2)/(2*c_2*k))-atan2(a_2,c_3);

end % if

% output cosmetics
% theta=theta+pi;
% theta=mod(real(theta),2*pi)+1i*imag(theta);
% theta=theta-pi;

end