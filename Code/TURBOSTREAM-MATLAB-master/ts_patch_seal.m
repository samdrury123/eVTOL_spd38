function [g] = ts_patch_seal(g,Tree,T,seal)
% Go through the entire seal mesh and patch it together, looking up
% coincident points with a kd-tree

% Patch all pitchwise periodic faces
for bid = seal.bids
    Q0 = [0 g{bid+1}.attribute.ni ; 0 g{bid+1}.attribute.nj ; 0 1];
    dir0 = [0 1 2];
    Q1 = [0 g{bid+1}.attribute.ni ; 0 g{bid+1}.attribute.nj ;...
        g{bid+1}.attribute.nk-1 g{bid+1}.attribute.nk];
    dir1 = [0 1 2];
    [ g ] = ts_add_periodic_patch(bid, bid, Q0, Q1, dir0, dir1, g, 5);
end

% Patch all surfaces that extend over entire faces
for bid = seal.bids
    % I1 - I = 0 faces
    XRRT = [g{bid+1}.x(1,1,1) ; g{bid+1}.r(1,1,1) ; g{bid}.rt(1,1,1)];
    index_vals = kdtree_range(Tree,double([XRRT - tolerance , XRRT + tolerance]));

    for index_val = index_vals
        if T(index_val,8) ~= bid && T(index_val,7) == 0
            IJK1a = T(index_val,1:3);
            bid1 = T(index_val,end-1);
        end
    end
    
    XRRT = [g{bid+1}.x(1,1,end) ; g{bid+1}.r(1,1,end) ; g{bid}.rt(1,1,end)];
    index_vals = kdtree_range(Tree,double([XRRT - tolerance , XRRT + tolerance]));

    for index_val = index_vals
        if T(index_val,8) ~= bid && T(index_val,7) == 0
            IJK1b = T(index_val,1:3);
        end
    end  
    
    XRRT = [g{bid+1}.x(1,end,1) ; g{bid+1}.r(1,end,1) ; g{bid}.rt(1,end,1)];
    index_vals = kdtree_range(Tree,double([XRRT - tolerance , XRRT + tolerance]));

    for index_val = index_vals
        if T(index_val,8) ~= bid && T(index_val,7) == 0
            IJK1c = T(index_val,1:3);
        end
    end
    
    XRRT = [g{bid+1}.x(1,end,end) ; g{bid+1}.r(1,end,end) ; g{bid}.rt(1,end,end)];
    index_vals = kdtree_range(Tree,double([XRRT - tolerance , XRRT + tolerance]));

    for index_val = index_vals
        if T(index_val,8) ~= bid && T(index_val,7) == 0
            IJK1d = T(index_val,1:3);
        end
    end      
    
    if isequal(IJK1a(1),IJK1a(1),IJK1a(1),IJK1a(1))
    
    Q0 = [0 1 ; 0 g{bid+1}.attribute.nj ; 0 g{bid+1}.attribute.nk];
    dir0 = [0 1 2];
    Q1 = [0 g{bid1+1}.attribute.ni ; 0 g{bid+1}.attribute.nj ;...
        g{bid+1}.attribute.nk-1 g{bid+1}.attribute.nk];
    dir1 = [0 1 2];
    [ g ] = ts_add_periodic_patch(bid, bid, Q0, Q1, dir0, dir1, g, 5);
    
    
    