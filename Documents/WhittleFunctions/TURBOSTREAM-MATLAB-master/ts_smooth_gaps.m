function g = ts_smooth_gaps(g,n_rows,xr_hub,xr_cas)
% TS_SMOOTH_GAPS  A function to smooth mesh coordinates in mid passage blocks

% Get o-mesh bids
[~,bid_os] = ts_get_row_bids(g,xr_hub,xr_cas);

% Loop over all rows and o-mesh patches
for n = n_rows
    for pid = ts_get_pids(g,bid_os(n))
        
        % Find mid passage blocks to smooth
        if g{bid_os(n)+1}.patch{pid+1}.attribute.idir == 0 && ...
                g{bid_os(n)+1}.patch{pid+1}.attribute.ken == g{bid_os(n)+1}.attribute.nk && ...
                g{bid_os(n)+1}.patch{pid+1}.attribute.kst == g{bid_os(n)+1}.attribute.nk-1
            bid = g{bid_os(n)+1}.patch{pid+1}.attribute.nxbid;
            
            % Loop over all internal mesh points
            for i = 2:g{bid+1}.attribute.ni-1
                for k = 2:g{bid+1}.attribute.nk-1
                    
                    % Apply smoothing to axial and tangential coordinates
                    g{bid+1}.x(i,:,k) = smooth(g{bid+1}.r(i,:,k),g{bid+1}.x(i,:,k));
                    g{bid+1}.rt(i,:,k) = smooth(g{bid+1}.r(i,:,k),g{bid+1}.rt(i,:,k));
                end
            end
        end        
    end
end


end