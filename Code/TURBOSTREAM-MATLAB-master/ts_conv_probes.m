function g = ts_conv_probes(g,seal_flag,mp_flag,xr_hub,xr_cas)
% Add steady probes to a data set to allow convergence of shroud wells and
% stages to be observed

% Where shroud wells join main passage, probe block just under hub line
if seal_flag == 1
    for bid = ts_get_bids(g)
        x = sum(sum(sum(g{bid+1}.x))) / numel(g{bid+1}.x);
        r = sum(sum(sum(g{bid+1}.r))) / numel(g{bid+1}.x);
        if r > interp1(xr_hub(:,1),xr_hub(:,2),x)
            for pid = ts_get_pids(g,bid)
                if g{bid+1}.patch{pid+1}.attribute.kind == 5
                    nxbid = g{bid+1}.patch{pid+1}.attribute.nxbid;
                    x = sum(sum(sum(g{nxbid+1}.x))) / numel(g{nxbid+1}.x);
                    r = sum(sum(sum(g{nxbid+1}.r))) / numel(g{nxbid+1}.x);
                    if r < interp1(xr_hub(:,1),xr_hub(:,2),x)
                        g{nxbid+1}.attribute.np = g{nxbid+1}.attribute.np + 1;
                        g{nxbid+1}.patch{end+1} = g{nxbid+1}.patch{end};
%                         g{nxbid+1}.patch{end}.probe_append = 0;
                        g{nxbid+1}.patch{end}.attribute.pid = length(g{nxbid+1}.patch) - 1;
                        g{nxbid+1}.patch{end}.attribute.kind = 20;
                        g{nxbid+1}.patch{end}.attribute.ist = 0;
                        g{nxbid+1}.patch{end}.attribute.ien = 1;
                        g{nxbid+1}.patch{end}.attribute.jst = 0;
                        g{nxbid+1}.patch{end}.attribute.jen = g{nxbid+1}.attribute.nj;
                        g{nxbid+1}.patch{end}.attribute.kst = 0;
                        g{nxbid+1}.patch{end}.attribute.ken = g{nxbid+1}.attribute.nk;
                        
                        break
                    end
                end
            end
        end
    end
end

% Probes on downstream sides of mixing planes, in the main passage only
if mp_flag == 1
    for bid = ts_get_bids(g)
        for pid = ts_get_pids(g,bid)
            if g{bid+1}.patch{pid+1}.attribute.kind == 2 && ...
                g{bid+1}.patch{pid+1}.attribute.ist == 0
                x = sum(sum(sum(g{bid+1}.x))) / numel(g{bid+1}.x);
                r = sum(sum(sum(g{bid+1}.r))) / numel(g{bid+1}.x);
                if r > interp1(xr_hub(:,1),xr_hub(:,2),x) && r < interp1(xr_cas(:,1),xr_cas(:,2),x)
                    g{bid+1}.attribute.np = g{bid+1}.attribute.np + 1;
                    g{bid+1}.patch{end+1} = g{bid+1}.patch{end};
%                     g{bid+1}.patch{end}.probe_append = 0;
                    g{bid+1}.patch{end}.attribute.pid = length(g{bid+1}.patch) - 1;
                    g{bid+1}.patch{end}.attribute.kind = 20;
                    g{bid+1}.patch{end}.attribute.ist = 0;
                    g{bid+1}.patch{end}.attribute.ien = 1;
                    g{bid+1}.patch{end}.attribute.jst = 0;
                    g{bid+1}.patch{end}.attribute.jen = g{bid+1}.attribute.nj;
                    g{bid+1}.patch{end}.attribute.kst = 0;
                    g{bid+1}.patch{end}.attribute.ken = g{bid+1}.attribute.nk;
                    
                    break
                end
            end
        end
    end
end


end