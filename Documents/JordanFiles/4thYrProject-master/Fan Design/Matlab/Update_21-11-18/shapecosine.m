function [x,z] = shapecosine(tmax,xmax,A,t_te)
% Function that finds the shape class of the type
% S = A * cos(2pi*x) + Bx + C
% that has maximum thickness tmax at x=xmax
% _______________________05/09/2014 Andrea Maffioli________________________
xle = linspace(0,0.05,81);
xbody = linspace(0.05,1,96);
x = [xle xbody(2:end)];
% Thickness at Trailing Edge
zte=t_te/2;
% Class function C(x) and Maximum thickness function F(x)
Cx = sqrt(x).*(1-x);
F = (tmax/2-x*zte)./Cx;
% Values of S(x) and dS/dx at x=xmax
Sxmax = (tmax/2-xmax*zte)/(sqrt(xmax)*(1-xmax));
dSdx = -zte/(sqrt(xmax)*(1-xmax))-(tmax/2-xmax*zte)*(1-3*xmax)/(2*xmax^(3/2)*(1-xmax)^2);
% Cosine wavenumber=1 in xmax range = [0.27 0.39], k=1/2 outside this range
k=1;
if xmax>0.39 %|| xmax<0.26
    k=1/2;
    A=-A;
end

% Find coefficients of shape function
B = dSdx + 2*pi*A*sin(2*pi*k*xmax);
C = Sxmax - A*cos(2*pi*k*xmax) - B*xmax;
S = A*cos(2*pi*k*x) + B*x + C;
% Find thickness distribution and plot
z = S.*Cx+x*zte;
% figure
%figure, plot(x,z), hold on, plot(x,-z), axis equal
%figure, plot(x,S), grid on, hold on%, plot(x,F), axis([0 1 0 0.3])
end