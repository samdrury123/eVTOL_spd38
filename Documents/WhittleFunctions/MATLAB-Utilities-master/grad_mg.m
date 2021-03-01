%Function to find the central difference gradient for 2 equal length
%vectors x,y with single sided at the extremes
%
%  dydx = grad_mg(x,y)
%
function [dydx, d2ydx2] = grad_mg(x,y)
 
m = length(x);

dydx = zeros(size(x));
if length(x) ~= length(y)
    error('vectors must be same length')
end
dydx(1) = (y(2) - y(1))/(x(2) - x(1));
for ii = 2:m-1
    d1 = x(ii+1) - x(ii);
    d2 = x(ii) - x(ii-1);
    % 1st derivative
    dydx(ii) = ( y(ii+1) - (y(ii-1)*((d1./d2).^2)) - (y(ii)*(1-((d1/d2)^2))) )./(d1.*(1+(d1/d2)));
            
    % 2nd derivative
    d2ydx2(ii) = ( y(ii+1) +  (y(ii-1)*(d1./d2).^3) - (y(ii)*(1+((d1/d2)^3)))...
        -  dydx(ii)*d1*(1-((d1/d2)^2))  )./(((d1.^2)/2).*(1+((d1/d2))));
    
    % 3rd derivative
    
    
end
%dydx(2:m-1) = (y(3:m)-y(1:m-2))./(x(3:m)-x(1:m-2));
dydx(m) = (y(m) - y(m-1))/(x(m) - x(m-1));

d2ydx2(1) = d2ydx2(2);
d2ydx2(m) = d2ydx2(m-1);


end