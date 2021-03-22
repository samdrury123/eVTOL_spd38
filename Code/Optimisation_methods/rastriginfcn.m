function f = rastriginfcn(x)
    n = size(x,2); % gets dimension of x
    A = 10; % constant of function
    f = A*n + sum((x.^2 - A * cos(2*pi*x)));
end