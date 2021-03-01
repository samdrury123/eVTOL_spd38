function S = hyperbolic_bunch( nPoints, Sp1, Sp2 )
%Generates a hyperbolically bunched points distribution

nPoints = double(nPoints);

A = sqrt(Sp1/Sp2);
b = fminbnd(@(b) abs(sinh(b)-b/((nPoints-1)*sqrt(Sp1*Sp2))),0,100);

S = zeros(1,nPoints); 
for i = 1:nPoints
    R = (i-1)/(nPoints-1) - 0.5;
    u = 1+tanh(b*R)/tanh(b/2);
    S(i) = u/(2*A + (1-A)*u);
end

end