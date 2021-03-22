function g = ts_calc_egen(g,bids,P_dyn)
% TS_CALC_EGEN  Calculate entropy generation and work input rates for individual cells

% Default to all blocks
if exist('bids','var') == 0 || isempty(bids) == 1
    bids = ts_get_bids(g);
end

% Calculate secondary variables
g = ts_secondary(g);

% Loop over all blocks 
for bid = bids
    
    % Calculate volumes and areas
    [V,A_norm] = ts_calc_vol(g{bid+1}.x,g{bid+1}.r,g{bid+1}.rt);
    
    % Extract face centred average values
    g{bid+1}.rovt = g{bid+1}.rorvt ./ g{bid+1}.r; 
    varnames = {'ro' 'rovx' 'rovr' 'rovt' 's' 'ho' 'Po_rel'};
    for v = 1:length(varnames)
        x = g{bid+1}.(varnames{v});
        x = reshape(cat(4,x(1:end-1,1:end-1,1:end-1),x(2:end,1:end-1,1:end-1),...
            x(1:end-1,2:end,1:end-1), x(1:end-1,1:end-1,2:end),...
            x(2:end,2:end,1:end-1), x(2:end,1:end-1,2:end),...
            x(1:end-1,2:end,2:end), x(2:end,2:end,2:end)),[],8);
        x = cat(2, sum(x(:,[1 4 7 3]),2), sum(x(:,[2 5 8 6]),2), sum(x(:,[1 2 6 4]),2), ...
            sum(x(:,[3 7 8 5]),2), sum(x(:,[1 3 5 2]),2), sum(x(:,[6 8 7 4]),2)) / 4;
        f.(varnames{v}) = reshape(x,[size(V) 1 6]);
    end
    
    % Calculate mass fluxes
    mdot = sum(cat(4,f.rovx,f.rovr,f.rovt) .* A_norm,4);
    
    % Calculate mass in the cell
    ro = mean(f.ro,5);
    m = ro .* V;
    
    % Calculate entropy generation rate
    g{bid+1}.sdot = zeros(size(g{bid+1}.x));
    g{bid+1}.sdot(1:end-1,1:end-1,1:end-1) = sum(mdot .* f.s,5) ./ m;
    
    % Calculate pressure loss coefficient
    yp = abs(sum(mdot .* f.Po_rel,5) ./ ((0.5 * sum(abs(mdot),5)).^2 * P_dyn));
    
    % Distribute back to nodes
    g{bid+1}.yp = zeros(size(g{bid+1}.x));
    g{bid+1}.yp(1:end-1,1:end-1,1:end-1) = yp;
%     g{bid+1}.yp(2:end-1,2:end-1,2:end-1) = 0.125 * ...
%         (yp(1:end-1,1:end-1,1:end-1) + yp(2:end,1:end-1,1:end-1) + ...
%         yp(1:end-1,2:end,1:end-1,:) + yp(1:end-1,1:end-1,2:end,:) + ...
%         yp(2:end,2:end,1:end-1,:) + yp(2:end,1:end-1,2:end,:) + ...
%         yp(1:end-1,2:end,2:end,:) + yp(2:end,2:end,2:end,:));
    
    
end



end