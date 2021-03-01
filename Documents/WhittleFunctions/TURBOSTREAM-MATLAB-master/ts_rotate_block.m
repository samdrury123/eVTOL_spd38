function g = ts_rotate_block(g, bid, dimensions)
% Rotates the specified block in the direction required by dimensions,
% rearranges patches to match up too

% Get block and property names
b = ts_get_block(bid,g);
varnames = fieldnames(b);

% Loop over properties and flip and permute as required
for v = 1:length(varnames)
    varname = varnames{v};
    if ndims(b.(varname)) == 3
        b.(varname) = permute(b.(varname),abs(dimensions));
        
        if sum(dimensions == -1) == 1
            b.(varname) = flipdim(b.(varname),1);
            rpmi1 = b.bv.rpmi1;
            b.bv.rpmi1 = b.bv.rpmi2;
            b.bv.rpmi2 = rpmi1;
        end
        if sum(dimensions == -2) == 1
            b.(varname) = flipdim(b.(varname),2);
            rpmj1 = b.bv.rpmj1;
            b.bv.rpmj1 = b.bv.rpmj2;
            b.bv.rpmj2 = rpmj1;
        end
        if sum(dimensions == -3) == 1
            b.(varname) = flipdim(b.(varname),3);
            rpmk1 = b.bv.rpmk1;
            b.bv.rpmk1 = b.bv.rpmk2;
            b.bv.rpmk2 = rpmk1;            
        end
    end
end

% Take care of attributes and block variables
n = [b.attribute.ni b.attribute.nj b.attribute.nk];
n = n(abs(dimensions));
b.attribute.ni = n(1);
b.attribute.nj = n(2);
b.attribute.nk = n(3);
rpms = [b.bv.rpmi1 b.bv.rpmi2 ; b.bv.rpmj1 b.bv.rpmj2 ; b.bv.rpmk1 b.bv.rpmk2];
rpms = rpms(abs(dimensions),:);
b.bv.rpmi1 = rpms(1,1);
b.bv.rpmi2 = rpms(1,2);
b.bv.rpmj1 = rpms(2,1);
b.bv.rpmj2 = rpms(2,2);
b.bv.rpmk1 = rpms(3,1);
b.bv.rpmk2 = rpms(3,2);

if sum(dimensions == -1) == 1
    rpmi1 = b.bv.rpmi1;
    b.bv.rpmi1 = b.bv.rpmi2;
    b.bv.rpmi2 = rpmi1;
end
if sum(dimensions == -2) == 1
    rpmj1 = b.bv.rpmj1;
    b.bv.rpmj1 = b.bv.rpmj2;
    b.bv.rpmj2 = rpmj1;
end
if sum(dimensions == -3) == 1
    rpmk1 = b.bv.rpmk1;
    b.bv.rpmk1 = b.bv.rpmk2;
    b.bv.rpmk2 = rpmk1;            
end

% Loop over patches and change directions as required
for pid = ts_get_pids(b)
    if b.patch{pid+1}.attribute.kind == 5
        % Simple I-J swap
        if dimensions(1) == 2 && dimensions(2) == 1
            nxbid = b.patch{pid+1}.attribute.nxbid;
            nxpid = b.patch{pid+1}.attribute.nxpid;

            if nxbid ~= bid
                jdir = b.patch{pid+1}.attribute.jdir;
                idir = b.patch{pid+1}.attribute.idir;
                b.patch{pid+1}.attribute.idir = jdir;
                b.patch{pid+1}.attribute.jdir = idir;

                if g{nxbid+1}.patch{nxpid+1}.attribute.idir == 0
                    g{nxbid+1}.patch{nxpid+1}.attribute.idir = 1;
                elseif g{nxbid+1}.patch{nxpid+1}.attribute.jdir == 1
                    g{nxbid+1}.patch{nxpid+1}.attribute.jdir = 0;
                end
            end
            
            ist = b.patch{pid+1}.attribute.ist;
            ien = b.patch{pid+1}.attribute.ien;
            jst = b.patch{pid+1}.attribute.jst;
            jen = b.patch{pid+1}.attribute.jen;
            
            b.patch{pid+1}.attribute.ist = jst;
            b.patch{pid+1}.attribute.ien = jen;
            b.patch{pid+1}.attribute.jst = ist;
            b.patch{pid+1}.attribute.jen = ien;            
        end
        % Rotate I-J swap
        if dimensions(1) == -2 && dimensions(2) == 1
            nxbid = b.patch{pid+1}.attribute.nxbid;
            nxpid = b.patch{pid+1}.attribute.nxpid;

            if nxbid ~= bid
                jdir = b.patch{pid+1}.attribute.jdir;
                idir = b.patch{pid+1}.attribute.idir;
                b.patch{pid+1}.attribute.idir = jdir;
                b.patch{pid+1}.attribute.jdir = 3;

                if g{nxbid+1}.patch{nxpid+1}.attribute.idir == 0
                    g{nxbid+1}.patch{nxpid+1}.attribute.idir = 4;
                elseif g{nxbid+1}.patch{nxpid+1}.attribute.jdir == 1
                    g{nxbid+1}.patch{nxpid+1}.attribute.jdir = 0;
                end
            end
            
            ist = b.patch{pid+1}.attribute.ist;
            ien = b.patch{pid+1}.attribute.ien;
            jst = b.patch{pid+1}.attribute.jst;
            jen = b.patch{pid+1}.attribute.jen;
            
            ni = b.attribute.ni;
            nj = b.attribute.nj;
            
            b.patch{pid+1}.attribute.ist = jst;
            b.patch{pid+1}.attribute.ien = jen;
            b.patch{pid+1}.attribute.jst = nj - ien;
            b.patch{pid+1}.attribute.jen = nj - ist;            
        end
    end
end

g = ts_add_block(b,g);

end