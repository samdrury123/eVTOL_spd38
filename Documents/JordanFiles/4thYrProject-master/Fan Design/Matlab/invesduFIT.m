function [x] = invesduFIT(y)
%invesduFIT(\sigma) = [hi/L]
%   Fit to data of x = L/hi; y = Ae/A3 = sigma
%   Relationship between the diffusion factor of the conical exit diffusor
%   and the aspect ratio of the diffusor as defined by the separation limit
%   in ESDU 75026
a = 0.3313;
b = 0.5275;
c = 1.06;

x = ((y-c)/a)^(-1/b);

end