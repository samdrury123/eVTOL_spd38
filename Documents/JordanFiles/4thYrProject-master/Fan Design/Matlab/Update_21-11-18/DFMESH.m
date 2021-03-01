function [df] = DFMESH(phi, psi)
%DF(phi, psi)
%   Fit to data of x = L/hi; y = Ae/A3 = sigma
%   Relationship between the diffusion factor of the conical exit diffusor
%   and the aspect ratio of the diffusor as defined by the separation limit
%   in ESDU 75026
% a = 0.2863;
% b = 0.5627;
% c = 1.06;
global n

a = 0.3313;
b = 0.5275;
c = 1;

sigma = sqrt((phi.^2)./(2.*psi));
df = zeros(n,n);

for idx = 1:numel(sigma)
    if sigma(idx) > 1
%         disp('Break');
        df(idx) = ((sigma(idx) - c)./a).^(-1/b);
    else
        df(idx) = 9999999999;
    end
end

% sigma = 1.4;
% df = 17.2776 - 0.00270343*(14326.2*sqrt(2.05239e14*sigma^2-9.72083e14*sigma+1.15761e15)-2.05239e11*sigma+4.86041e11)^(1/3)+29878.6/(14326.2*sqrt(2.05239e14*sigma^2-9.72083e14*sigma+1.15761e15)-2.05239e11*sigma+4.86041e11)^(1/3);
% df = 9.23225*(sigma - 1.06)^1.777145903678692;
end