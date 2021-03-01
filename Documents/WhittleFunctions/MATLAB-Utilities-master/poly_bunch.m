function s = poly_bunch(s1,s2,ds1,ds2,ni,er)
% POLY_BUNCH Bunch a points distribution to a specified spacing and value

% Default to cubic fit with free intiial expansion ratio
if exist('er','var') == 0; er = 0; end;

% Point index as x-values
x1 = 1; x2 = ni;

% Determine order of polynomial fit
if er == 0
    
    % Fit a polynomial through four available points
    p = polyfit([1 2 ni-1 ni],[s1 s1+ds1 s2-ds2 s2],3);
    
else
    
    % Fit a higher order polynomial that changes more slowly
    p = polyfit([1 2 3 ni-2 ni-1 ni],[s1 s1+ds1 s1+(1+er)*ds1 s2-(1+er)*ds2 s2-ds2 s2],5);
    
end

% Evaluate polynomial at desired resolution
s = polyval(p,x1:x2);


end

