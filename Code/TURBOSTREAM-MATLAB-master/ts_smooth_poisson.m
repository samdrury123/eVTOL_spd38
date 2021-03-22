function g = ts_smooth_poisson(g)
% TS_SMOOTH_POISSON  Smooth out any wiggles in a poisson solution of the wall distance
%
%   g = TS_SMOOTH_POISSON(g)
%
%   g - cell array of TURBOSTREAM solution

% Remove negative poisson results
for bid = ts_get_bids(g)
    g{bid+1}.phi(g{bid+1}.phi < 0) = 0;
end

% Smooth poisson solution keeping wall values the same
for bid = ts_get_bids(g)
    v = g{bid+1}.phi;
    g{bid+1}.phi = smooth3(v,'box',5);
    q = g{bid+1}.mwall == 0;
    g{bid+1}.phi(q) = v(q);
end
    
end