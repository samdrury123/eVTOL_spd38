function g = ts_remove_probes(g)
% TS_REMOVE_PROBES  Remove all probe patches from a TURBOSTREAM mesh
%
%   g = TS_REMOVE_PROBES(g)
%
%   g - cell array of TURBOSTREAM solution
%
%   Note: Only works when probes are final patches on a block

% Loop over all blocks
for bid = ts_get_bids(g)
    for pid = flip(ts_get_pids(g,bid),2)
    
        % Check kind of patch and remove
        if g{bid+1}.patch{pid+1}.attribute.kind == 8
            g{bid+1}.patch(pid+1) = [];
            g{bid+1}.attribute.np = g{bid+1}.attribute.np - 1;
        end
        
    end
end

end