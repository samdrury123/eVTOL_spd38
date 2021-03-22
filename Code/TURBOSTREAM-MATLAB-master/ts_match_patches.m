function g = ts_match_patches(g)
% Alter mesh coordinates on periodic patches to ensure continuity

% Loop over all blocks and patches
for bid = ts_get_bids(g)
    for pid = ts_get_pids(g,bid)
        if g{bid+1}.patch{pid+1}.attribute.kind == 5
        
            % Pull out current and matching blocks and patches
            b0 = g{bid+1};
            p0 = b0.patch{pid+1}.attribute;
            b1 = g{p0.nxbid+1};
            p1 = b1.patch{p0.nxpid+1}.attribute;

            % Current coordinates
            I = p0.ist+1:p0.ien; J = p0.jst+1:p0.jen; K = p0.kst+1:p0.ken;
            p0.x = squeeze(b0.x(I,J,K)); p0.r = squeeze(b0.r(I,J,K));
            p0.rt = squeeze(b0.rt(I,J,K)); p0.t = p0.rt ./ p0.r;
            p0.ni = p0.ien - p0.ist; p0.nj = p0.jen - p0.jst; p0.nk = p0.ken - p0.kst;
            
            % Matching coordinates
            p1 = ts_get_match(p0,p1,b1);
            
            % Average meridional coordinates on both patches
            x = 0.5 * (p0.x + p1.x); r = 0.5 * (p0.r + p1.r);
            
            % Correct tangential coordinates by a whole number of pitches
            t0 = mean(p0.t(:)); t1 = mean(p1.t(:)); 
            dt = 2 * pi / double(g{bid+1}.bv.nblade);
            dt = round((t1 - t0) / dt) * dt;
            p1.t = p1.t - dt;
            t = 0.5 * (p0.t + p1.t);
            
            % Overwrite coordinates in first block with average values
            g{bid+1}.x(I,J,K) = x; g{bid+1}.r(I,J,K) = r; g{bid+1}.rt(I,J,K) = r .* t; 
            
            % Get indices of the matching block
            q = sub2ind([b1.attribute.ni b1.attribute.nj b1.attribute.nk],p1.I,p1.J,p1.K);
            
            % Overwrite coordinates in the matching block with average values
            g{p0.nxbid+1}.x(q) = x; g{p0.nxbid+1}.r(q) = r; g{p0.nxbid+1}.rt(q) = r .* (t + dt); 
                
        end
    end
end

end