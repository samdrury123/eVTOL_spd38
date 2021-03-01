function g = ts_downsample(g);
% TS_DOWNSAMPLE  Downsample a mesh and all patches in TURBOSTREAM

% Downsample all coordinates and match up patches
for bid = ts_get_bids(g)

    % Downsample block properties
    varnames = fieldnames(g{bid+1});
    for v = 1:length(varnames)
        if isnumeric(g{bid+1}.(varnames{v})) == 1
            g{bid+1}.(varnames{v}) = g{bid+1}.(varnames{v})(1:2:end,1:2:end,1:2:end);
        end
    end

    % Patch attributes
    varnames = {'ist' 'ien' 'jst' 'jen' 'kst' 'ken'};
    for pid = ts_get_pids(g,bid)
        for v = 1:length(varnames)
            g{bid+1}.patch{pid+1}.attribute.(varnames{v}) = ...
                round(g{bid+1}.patch{pid+1}.attribute.(varnames{v}) / 2);
        end
    end

    % Inlet patch
    varnames = {'pstag' 'tstag' 'yaw' 'pitch' 'ym0'};
    for pid = ts_get_pids(g,bid)
        if g{bid+1}.patch{pid+1}.attribute.kind == 0
            nj = g{bid+1}.attribute.nj; nk = g{bid+1}.attribute.nk;
            j = repmat((1:2:nj)',[1 round(nk/2)]); k = repmat((1:2:nk),[round(nj/2) 1]);
            n = reshape(sub2ind([nj nk],j,k),[],1);
            for v = 1:length(varnames)
                if isfield(g{bid+1}.patch{pid+1},varnames{v}) == 1
                    g{bid+1}.patch{pid+1}.(varnames{v}) = g{bid+1}.patch{pid+1}.(varnames{v})(n);
                end
            end
        end
    end

    % Block attributes
    varnames = {'ni' 'nj' 'nk'};
    for v = 1:length(varnames)
        g{bid+1}.attribute.(varnames{v}) = round(g{bid+1}.attribute.(varnames{v}) / 2);
    end
end
    

end