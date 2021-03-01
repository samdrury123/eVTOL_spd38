function [] = ts_rev_rot(filename)

g = ts_read_hdf5(filename);

varnames = fieldnames(g{1});

for bid = ts_get_bids(g)
    % Flip coordinates of geometry
    g{bid+1}.rorvt = - g{bid+1}.rorvt;
    g{bid+1}.rt = - g{bid+1}.rt;
    
    % Determine whether to flip i or k directions
    flip = 'k';
    for pid = ts_get_pids(g,bid)
        if g{bid+1}.patch{pid+1}.attribute.kind == 5;
            if g{bid+1}.patch{pid+1}.attribute.nxbid == bid;
                ist = g{bid+1}.patch{pid+1}.attribute.ist;
                ien = g{bid+1}.patch{pid+1}.attribute.ien;
                if abs(ist - ien) == 1
                    flip = 'i';
                end
            end
        end
    end
    
    % Flip dimensions to ensure left handed set
    for v = 1:length(varnames)
        varname = varnames{v};
        if size(g{bid+1}.(varname),3) > 1
            if strcmp(flip,'k') == 1
                g{bid+1}.(varname) = flipdim(g{bid+1}.(varname),3);
            else
                g{bid+1}.(varname) = flipdim(g{bid+1}.(varname),1);
            end
        end
    end
    
    % Change patch starts and ends if necessary
    if strcmp(flip,'k') == 1
        nk = g{bid+1}.attribute.nk;
        for pid = ts_get_pids(g,bid)
            kst = g{bid+1}.patch{pid+1}.attribute.kst;
            ken = g{bid+1}.patch{pid+1}.attribute.ken;

            if kst ~= 0 || ken ~= nk
                g{bid+1}.patch{pid+1}.attribute.kst = nk - ken;
                g{bid+1}.patch{pid+1}.attribute.ken = nk - kst;
            end
        end
    else
        ni = g{bid+1}.attribute.ni;
        for pid = ts_get_pids(g,bid)
            ist = g{bid+1}.patch{pid+1}.attribute.ist;
            ien = g{bid+1}.patch{pid+1}.attribute.ien;

            if ist ~= 0 || ien ~= ni
                g{bid+1}.patch{pid+1}.attribute.ist = ni - ien;
                g{bid+1}.patch{pid+1}.attribute.ien = ni - ist;
            end
        end        
    end
        
    % Flip rpms to spin the right way
    g{bid+1}.bv.rpm = - g{bid+1}.bv.rpm;
    g{bid+1}.bv.rpmi1 = - g{bid+1}.bv.rpmi1;
    g{bid+1}.bv.rpmi2 = - g{bid+1}.bv.rpmi2;
    g{bid+1}.bv.rpmj1 = - g{bid+1}.bv.rpmj1;
    g{bid+1}.bv.rpmj2 = - g{bid+1}.bv.rpmj2;
    g{bid+1}.bv.rpmk1 = - g{bid+1}.bv.rpmk1;
    g{bid+1}.bv.rpmk2 = - g{bid+1}.bv.rpmk2;    
    
    % Change whirl angles at inlet
    for pid = ts_get_pids(g,bid)
        if g{bid+1}.patch{pid+1}.attribute.kind == 0
            g{bid+1}.patch{pid+1}.yaw = - g{bid+1}.patch{pid+1}.yaw;
        end
    end
end

% Reset patch directions
g = ts_set_dir(g);

ts_check_patches(g,2)

ts_write_hdf5(g,filename);

end