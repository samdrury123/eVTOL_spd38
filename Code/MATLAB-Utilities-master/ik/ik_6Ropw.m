function [ theta ] = ik_6Ropw( robot_parameter, in1, in2, phi)
%theta = ik_6Ropw( robot_parameter, EE ) for matrix-representation or
%theta = ik_6Ropw( robot_parameter, u, d ) for direction representation
%theta = ik_6Ropw( robot_parameter, u, d, phi ) for direction-angle representation
%
%   Computes the inverse kinematics of serial manipulators with an
%   ortho-parallel basis and a sperical wrist.
%
%   theta ... a matrix with all four possible solutions
%   robot_parameter ... [a1, a2, b1, c1, c2, c3, c4]
%   EE ... homogeneous end-effector matrix (first row: [1,0,0,0])
%   u ... vector from the basis to the end-effector
%   d ... vector of the direction of the end-effector
%   phi ... angle about z-axis [rad]
%   (c) Arthur Angerer and Mathias Brandstötter | IACE | UMIT | 2014 | Version 1.2
%
%   Associated paper: M. Brandstötter, A. Angerer, and M. Hofbaur, "An
%   Analytical Solution of the Inverse Kinematics Problem of Industrial
%   Serial Manipulators with an Ortho-parallel Basis and a Spherical
%   Wrist," Austrian Robotics Workshop, 2014.

if nargin < 2
    error('Error: Not correct number of inputs arguments.')
end % if

switch nargin
    case 2 % The end-effector matrix is given
        if (size(in1,1))==4
            EE = in1;
        else
            error('Error: EE is not a homogeneous 4x4 matrix.')
        end % if
        
    case 3 % The end-effector position and direction is given
        if size(in1) ~= 3
            error('Error: Not correct number of elements in input variable "u".')
        end % if
        
        if size(in2) ~= 3
            error('Error: Not correct number of elements in input variable "d".')
        end % if
        
        d=in2;
        
        if d(1)==0 && d(2)==0 && d(3)==0
            error('Error: d is not a direction.')
        end % if
        
        u=in1;
        
        if d(1)==0 && d(2)==0
            
            EE = [1, 0, 0, 0; u(1), d(3), 0, 0; u(2), 0, 1, 0; u(3), 0, 0, d(3)]; % phi=0          
            
        else
        
            d=d/norm(d);
            dn=sqrt(d(1)^2+d(2)^2);

            EE = [1,0,0,0;
                  u(1), (d(1)*d(3))/dn, -d(2)/dn, d(1);
                  u(2), (d(2)*d(3))/dn,  d(1)/dn, d(2);
                  u(3), -dn, 0, d(3)];
          
        end % if
          
    case 4 % The end-effector position, direction, and twist is given
        if size(in1) ~= 3
            error('Error: Not correct number of elements in input variable "u".')
        end % if
        
        if size(in2) ~= 3
            error('Error: Not correct number of elements in input variable "d".')
        end % if

        d=in2;
        
        if d(1)==0 && d(2)==0 && d(3)==0
            error('Error: Input variable "d" is not a direction.')
        end % if
        
        u=in1;
        
        if d(1)==0 && d(2)==0
            
            EE = [1, 0, 0, 0; u(1), cos(phi), -sin(phi), 0; u(2), sin(phi), cos(phi), 0; u(3), 0, 0, 1];           
            
        else
            d=d/norm(d);
            dn=sqrt(d(1)^2+d(2)^2);
         
            EE = [1, 0, 0, 0;
             u(1),  d(1)*d(3)/dn*cos(phi)-d(2)/dn*sin(phi), -d(1)*d(3)/dn*sin(phi)-d(2)/dn*cos(phi), d(1);
             u(2),  d(1)/dn*sin(phi)+d(2)*d(3)/dn*cos(phi),  d(1)/dn*cos(phi)-d(2)*d(3)/dn*sin(phi), d(2);
             u(3), -dn*cos(phi), dn*sin(phi), d(3)];
        end % if
         
    otherwise
        error('Error: Not correct number of inputs arguments.')
end % switch

if abs(1-det(EE)) > 1e-4
    error('Error: EE is not a special orthogonal matrix. The determinant of EE has to be +1.')
end % if

if size(robot_parameter) ~= 7
    error('Error: Not correct number of elements in input variable "robot_parameter".')
end % if


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN COMPUTATION OF THE INVERSE KINEMATICS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculation of point C

c4 = robot_parameter(7);

C = [1 0 0 0; EE(2:4,1)-EE(2:4,4)*c4, EE(2:4,2:4)];

% Computing the first three angles

theta123 = ik_3Ropw( robot_parameter, C(2:4,1) );
ind = find([isreal(theta123(:,1)) isreal(theta123(:,2)) isreal(theta123(:,3)) isreal(theta123(:,4))]==1); % indizes where theta123 is real
% Computing the second three angles

c23=cos(theta123(2,ind)+theta123(3,ind));
s23=sin(theta123(2,ind)+theta123(3,ind));

r11=EE(2,2); r12=EE(2,3); r13=EE(2,4);
r21=EE(3,2); r22=EE(3,3); r23=EE(3,4);
r31=EE(4,2); r32=EE(4,3); r33=EE(4,4);

M1 = ones(3,4)*1i; % alloction 

if not(isempty(ind))
    M1(1,ind) = atan2(-r13*sin(theta123(1,ind)) + r23*cos(theta123(1,ind)), r13*cos(theta123(1,ind)).*c23 + r23*sin(theta123(1,ind)).*c23 - r33*s23);
    M1(2,ind) = atan2(sqrt(1-(r13*cos(theta123(1,ind)).*s23+r23*sin(theta123(1,ind)).*s23+r33*c23).^2), r13*cos(theta123(1,ind)).*s23+r23*sin(theta123(1,ind)).*s23+r33*c23);
    M1(3,ind) = atan2(r12*cos(theta123(1,ind)).*s23 + r22*sin(theta123(1,ind)).*s23 + r32*c23, -r11*cos(theta123(1,ind)).*s23 - r21*sin(theta123(1,ind)).*s23 - r31*c23);
end % if

M2(1,:) =  M1(1,:) + pi;
M2(2,:) = -M1(2,:);
M2(3,:) =  M1(3,:) - pi;

% Combine all possible solutions

theta = [theta123(:,1) theta123(:,1) theta123(:,2) theta123(:,2) theta123(:,3) theta123(:,3) theta123(:,4) theta123(:,4);
         M1(:,1) M2(:,1) M1(:,2) M2(:,2) M1(:,3) M2(:,3) M1(:,4) M2(:,4)];

% Output cosmetics
theta=theta+pi;
theta=mod(real(theta),2*pi)+1i*imag(theta);
theta=theta-pi;

end % function