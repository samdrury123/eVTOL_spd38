function bid_os = ts_find_os(g,r)
% TS_FIND_OS  Find block numbers of all O-mesh blocks in a TURBOSTREAM mesh

% Find radii if unspecified
if exist('r','var') == 0
    rmax = -inf;
    rmin = inf;

    for bid = ts_get_bids(g)
        rmax = max(max(max(max(g{bid+1}.r))),rmax);
        rmin = min(min(min(min(g{bid+1}.r))),rmin);
    end
else
    rmax = r(2);
    rmin = r(1);
end

% Initialise o-mesh vector
bid_os = [];

% Loop over all blocks and find o-mesh that is peridiodic to itself and spans 50% of the radii given
for bid = ts_get_bids(g)
    for pid = ts_get_pids(bid,g)
        if g{bid+1}.patch{pid+1}.attribute.kind == 5
            if g{bid+1}.patch{pid+1}.attribute.nxbid == g{bid+1}.patch{pid+1}.attribute.bid
                if g{bid+1}.patch{pid+1}.attribute.ist == g{bid+1}.attribute.ni - 1 && ...
                        g{bid+1}.patch{pid+1}.attribute.ien == g{bid+1}.attribute.ni
                    if (g{bid+1}.r(1,end,1) - g{bid+1}.r(1,1,1)) > 0.5*(rmax-rmin)
                        bid_os = [bid_os bid];
                    end
                end
            end
        end
    end
end



end