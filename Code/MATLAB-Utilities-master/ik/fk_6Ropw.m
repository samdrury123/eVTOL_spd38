function EE = fk_6Ropw( robot_parameter, theta )
%theta = fk_6Ropw( robot_parameter, theta)
%   Computes the forward kinematics of serial manipulators with an
%   ortho-parallel basis.
%   theta is a vector or matrix with joint angles [theta1, ..., theta6]'.
%   (c) Mathias Brandstötter, Arthur Angerer | IACE | UMIT | 2014 | Version 1.2
%
%   Associated paper: M. Brandstötter, A. Angerer, and M. Hofbaur, "An
%   Analytical Solution of the Inverse Kinematics Problem of Industrial
%   Serial Manipulators with an Ortho-parallel Basis and a Spherical
%   Wrist," Austrian Robotics Workshop, 2014.

a1 = robot_parameter(1);
a2 = robot_parameter(2);
b = robot_parameter(3);
c1 = robot_parameter(4);
c2 = robot_parameter(5);
c3 = robot_parameter(6);

q1 = theta(1,:);
q2 = theta(2,:);
q3 = theta(3,:);
q4 = theta(4,:);
q5 = theta(5,:);
q6 = theta(6,:);

% forward kinematics (positioning part)

psi3 = atan2(a2,c3);
k = sqrt(a2^2+c3^2);

cx1 = c2*sin(q2)+k*sin(q2+q3+psi3)+a1;
cy1 = b;
cz1 = c2*cos(q2)+k*cos(q2+q3+psi3);

cx0 = cx1.*cos(q1) - cy1.*sin(q1);
cy0 = cx1.*sin(q1) + cy1.*cos(q1);
cz0 = cz1 + c1;

% forward kinematics (orientation part)

s1 = sin(q1);
s2 = sin(q2);
s3 = sin(q3);
s4 = sin(q4);
s5 = sin(q5);
s6 = sin(q6);

c1 = cos(q1);
c2 = cos(q2);
c3 = cos(q3);
c4 = cos(q4);
c5 = cos(q5);
c6 = cos(q6);

R_0c = [c1*c2*c3-c1*s2*s3  -s1  c1*c2*s3+c1*s2*c3;
        s1*c2*c3-s1*s2*s3   c1  s1*c2*s3+s1*s2*c3
        -s2*c3-c2*s3           0  -s2*s3+c2*c3];
    
R_ce = [c4*c5*c6-s4*s6  -c4*c5*s6-s4*c6  c4*s5;
         s4*c5*c6+c4*s6  -s4*c5*s6+c4*c6  s4*s5;
        -s5*c6            s5*s6           c5];
    
R_0e = R_0c*R_ce;

c4 = robot_parameter(7); % dont be confused with c4=cos(theta4)
u = [cx0 cy0 cz0]'+ c4*R_0e*[0 0 1]';

EE = [1 0 0 0; u R_0e];

end