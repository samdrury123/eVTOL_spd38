function [g] = ts_patch_row(g,Tree,T,blade_row)
% Patch a blade row together given data about block locations

tolerance = 1e-5;

% Check prescence of inner inlet & exit blocks
if isfield(blade_row,'inner_inlet')
    blade_row.outer_inlet = blade_row.inlet;
    blade_row.inlet = blade_row.inner_inlet;
end
if isfield(blade_row,'inner_exit')
    blade_row.outer_exit = blade_row.exit;
    blade_row.exit = blade_row.inner_exit;
end

% Periodic patches across passage
% Inlet
Q0 = [0 g{blade_row.inlet+1}.attribute.ni ; 0 g{blade_row.inlet+1}.attribute.nj ;...
    0 1];
dir0 = [0 1 2];
Q1 = [0 g{blade_row.inlet+1}.attribute.ni ; 0 g{blade_row.inlet+1}.attribute.nj ;...
    g{blade_row.inlet+1}.attribute.nk-1 g{blade_row.inlet+1}.attribute.nk];
dir1 = [0 1 2];
[ g ] = ts_add_periodic_patch(blade_row.inlet, blade_row.inlet, Q0, Q1, dir0, dir1, g, 5);

% Exit
Q0 = [0 g{blade_row.exit+1}.attribute.ni ; 0 g{blade_row.exit+1}.attribute.nj ;...
    0 1];
dir0 = [0 1 2];
Q1 = [0 g{blade_row.exit+1}.attribute.ni ; 0 g{blade_row.exit+1}.attribute.nj ;...
    g{blade_row.exit+1}.attribute.nk-1 g{blade_row.exit+1}.attribute.nk];
dir1 = [0 1 2];
[ g ] = ts_add_periodic_patch(blade_row.exit, blade_row.exit, Q0, Q1, dir0, dir1, g, 5);

% Upper-Lower
Q0 = [0 g{blade_row.lower+1}.attribute.ni ; 0 g{blade_row.lower+1}.attribute.nj ;...
    0 1];
dir0 = [0 1 2];
Q1 = [0 g{blade_row.upper+1}.attribute.ni ; 0 g{blade_row.upper+1}.attribute.nj ;...
    g{blade_row.upper+1}.attribute.nk-1 g{blade_row.upper+1}.attribute.nk];
dir1 = [0 1 2];
[ g ] = ts_add_periodic_patch(blade_row.lower, blade_row.upper, Q0, Q1, dir0, dir1, g, 5);

% Find corner points of omesh
% P1 - front of lower block
XRRT = [g{blade_row.lower+1}.x(1,1,end) ; g{blade_row.lower+1}.r(1,1,end) ; g{blade_row.lower+1}.rt(1,1,end)];
index_vals = kdtree_range(Tree,double([XRRT - tolerance , XRRT + tolerance]));

for index_val = index_vals
    if T(index_val,8) == blade_row.o
        C1 = T(index_val,1);
    end
    if T(index_val,8) == blade_row.inlet
        In1 = T(index_val,3);
    end
end

% P2 rear of lower block
XRRT = [g{blade_row.lower+1}.x(end,1,end) ; g{blade_row.lower+1}.r(end,1,end) ; g{blade_row.lower+1}.rt(end,1,end)];
index_vals = kdtree_range(Tree,double([XRRT - tolerance , XRRT + tolerance]));

for index_val = index_vals
    if T(index_val,8) == blade_row.o
        C2 = T(index_val,1);
    end
    if T(index_val,8) == blade_row.exit
        Ex1 = T(index_val,3);
    end
end

% P3 front of upper block
XRRT = [g{blade_row.upper+1}.x(1,1,1) ; g{blade_row.upper+1}.r(1,1,1) ; g{blade_row.upper+1}.rt(1,1,1)];
index_vals = kdtree_range(Tree,double([XRRT - tolerance , XRRT + tolerance]));

for index_val = index_vals
    if T(index_val,8) == blade_row.o
        C3 = T(index_val,1);
    end
    if T(index_val,8) == blade_row.inlet
        In2 = T(index_val,3);
    end
end

% P3 rear of upper block
XRRT = [g{blade_row.upper+1}.x(end,1,1) ; g{blade_row.upper+1}.r(end,1,1) ; g{blade_row.upper+1}.rt(end,1,1)];
index_vals = kdtree_range(Tree,double([XRRT - tolerance , XRRT + tolerance]));

for index_val = index_vals
    if T(index_val,8) == blade_row.o
        C4 = T(index_val,1);
    end
    if T(index_val,8) == blade_row.exit
        Ex2 = T(index_val,3);
    end
end

% Seam of O-mesh
XRRT = [g{blade_row.o+1}.x(1,1,end) ; g{blade_row.o+1}.r(1,1,end) ; g{blade_row.o+1}.rt(1,1,end)];
index_vals = kdtree_range(Tree,double([XRRT - tolerance , XRRT + tolerance]));

for index_val = index_vals
    if T(index_val,8) == blade_row.inlet
        In3 = T(index_val,3);
    end
end

% Patch O-mesh to

% Inlet 1
Q0 = [g{blade_row.inlet+1}.attribute.ni-1 g{blade_row.inlet+1}.attribute.ni ; 0 g{blade_row.inlet+1}.attribute.nj ;...
    In1 In3+1];
dir0 = [-1 1 0];
Q1 = [C1 g{blade_row.o+1}.attribute.ni ; 0 g{blade_row.o+1}.attribute.nj ;...
    g{blade_row.o+1}.attribute.nk-1 g{blade_row.o+1}.attribute.nk];
dir1 = [2 1 -1];
[ g ] = ts_add_periodic_patch(blade_row.inlet, blade_row.o, Q0, Q1, dir0, dir1, g, 5);

% Inlet 2
Q0 = [g{blade_row.inlet+1}.attribute.ni-1 g{blade_row.inlet+1}.attribute.ni ; 0 g{blade_row.inlet+1}.attribute.nj ;...
    In3 In2+1];
dir0 = [-1 1 0];
Q1 = [0 C3+1 ; 0 g{blade_row.o+1}.attribute.nj ;...
    g{blade_row.o+1}.attribute.nk-1 g{blade_row.o+1}.attribute.nk];
dir1 = [2 1 -1];
[ g ] = ts_add_periodic_patch(blade_row.inlet, blade_row.o, Q0, Q1, dir0, dir1, g, 5);

% Exit
Q0 = [0 1 ; 0 g{blade_row.exit+1}.attribute.nj ;...
    Ex1 Ex2+1];
dir0 = [-1 1 3];
Q1 = [C4 C2+1 ; 0 g{blade_row.o+1}.attribute.nj ;...
    g{blade_row.o+1}.attribute.nk-1 g{blade_row.o+1}.attribute.nk];
dir1 = [5 1 -1];
[ g ] = ts_add_periodic_patch(blade_row.exit, blade_row.o, Q0, Q1, dir0, dir1, g, 5);

% Lower
Q0 = [0 g{blade_row.lower+1}.attribute.ni ; 0 g{blade_row.lower+1}.attribute.nj ;...
    g{blade_row.lower+1}.attribute.nk-1 g{blade_row.lower+1}.attribute.nk];
dir0 = [3 1 -1];
Q1 = [C2 C1+1 ; 0 g{blade_row.o+1}.attribute.nj ;...
    g{blade_row.o+1}.attribute.nk-1 g{blade_row.o+1}.attribute.nk];
dir1 = [3 1 -1];
[ g ] = ts_add_periodic_patch(blade_row.lower, blade_row.o, Q0, Q1, dir0, dir1, g, 5);

% Upper
Q0 = [0 g{blade_row.upper+1}.attribute.ni ; 0 g{blade_row.upper+1}.attribute.nj ;...
    0 1];
dir0 = [0 1 -1];
Q1 = [C3 C4+1 ; 0 g{blade_row.o+1}.attribute.nj ;...
    g{blade_row.o+1}.attribute.nk-1 g{blade_row.o+1}.attribute.nk];
dir1 = [0 1 -1];
[ g ] = ts_add_periodic_patch(blade_row.upper, blade_row.o, Q0, Q1, dir0, dir1, g, 5);

% Itself
Q0 = [0 1 ; 0 g{blade_row.o+1}.attribute.nj ;...
    0 g{blade_row.o+1}.attribute.nk];
dir0 = [-1 1 2];
Q1 = [g{blade_row.o+1}.attribute.ni-1 g{blade_row.o+1}.attribute.ni ; 0 g{blade_row.o+1}.attribute.nj ;...
    0 g{blade_row.o+1}.attribute.nk];
dir1 = [-1 1 2];
[ g ] = ts_add_periodic_patch(blade_row.o, blade_row.o, Q0, Q1, dir0, dir1, g, 5);


% Patch inlet block to

% Lower
Q0 = [0 1 ; 0 g{blade_row.lower+1}.attribute.nj ;...
    0 g{blade_row.lower+1}.attribute.nk];
dir0 = [-1 1 2];
Q1 = [g{blade_row.inlet+1}.attribute.ni-1 g{blade_row.inlet+1}.attribute.ni ; 0 g{blade_row.inlet+1}.attribute.nj ;...
    0 In1+1];
dir1 = [-1 1 2];
[ g ] = ts_add_periodic_patch(blade_row.lower, blade_row.inlet, Q0, Q1, dir0, dir1, g, 5);

% Upper
Q0 = [0 1 ; 0 g{blade_row.upper+1}.attribute.nj ;...
    0 g{blade_row.upper+1}.attribute.nk];
dir0 = [-1 1 2];
Q1 = [g{blade_row.inlet+1}.attribute.ni-1 g{blade_row.inlet+1}.attribute.ni ; 0 g{blade_row.inlet+1}.attribute.nj ;...
    In2 g{blade_row.inlet+1}.attribute.nk];
dir1 = [-1 1 2];
[ g ] = ts_add_periodic_patch(blade_row.upper, blade_row.inlet, Q0, Q1, dir0, dir1, g, 5);


% Patch exit block to

% Lower
Q0 = [g{blade_row.lower+1}.attribute.ni-1 g{blade_row.lower+1}.attribute.ni ; 0 g{blade_row.lower+1}.attribute.nj ;...
    0 g{blade_row.lower+1}.attribute.nk];
dir0 = [-1 1 2];
Q1 = [0 1 ; 0 g{blade_row.exit+1}.attribute.nj ;...
    0 Ex1+1];
dir1 = [-1 1 2];
[ g ] = ts_add_periodic_patch(blade_row.lower, blade_row.exit, Q0, Q1, dir0, dir1, g, 5);

% Upper
Q0 = [g{blade_row.upper+1}.attribute.ni-1 g{blade_row.upper+1}.attribute.ni ; 0 g{blade_row.upper+1}.attribute.nj ;...
    0 g{blade_row.upper+1}.attribute.nk];
dir0 = [-1 1 2];
Q1 = [0 1 ; 0 g{blade_row.exit+1}.attribute.nj ;...
    Ex2 g{blade_row.exit+1}.attribute.nk];
dir1 = [-1 1 2];
[ g ] = ts_add_periodic_patch(blade_row.upper, blade_row.exit, Q0, Q1, dir0, dir1, g, 5);



% Deal with tip gap
if isfield(blade_row,'gap_o')
    % O to Tip-O
    Q0 = [0 g{blade_row.gap_o+1}.attribute.ni ; 0 g{blade_row.gap_o+1}.attribute.nj ;...
        g{blade_row.gap_o+1}.attribute.nk-1 g{blade_row.gap_o+1}.attribute.nk];
    dir0 = [0 1 -1];
    Q1 = [0 g{blade_row.o+1}.attribute.ni ; g{blade_row.o+1}.attribute.nj-g{blade_row.gap_o+1}.attribute.nj g{blade_row.o+1}.attribute.nj ;...
        0 1];
    dir1 = [0 1 -1];
    [ g ] = ts_add_periodic_patch(blade_row.gap_o, blade_row.o, Q0, Q1, dir0, dir1, g, 5);
    
    % Find corner points of tip gap o mesh
    % G1 i=0 k=0
    XRRT = [g{blade_row.gap_h+1}.x(1,1,1) ; g{blade_row.gap_h+1}.r(1,1,1) ; g{blade_row.gap_h+1}.rt(1,1,1)];
    index_vals = kdtree_range(Tree,double([XRRT - tolerance , XRRT + tolerance]));

    for index_val = index_vals
        if T(index_val,8) == blade_row.gap_o
            G1 = T(index_val,1);
        end
    end

    % G2 i=ni k=0
    XRRT = [g{blade_row.gap_h+1}.x(end,1,1) ; g{blade_row.gap_h+1}.r(end,1,1) ; g{blade_row.gap_h+1}.rt(end,1,1)];
    index_vals = kdtree_range(Tree,double([XRRT - tolerance , XRRT + tolerance]));

    for index_val = index_vals
        if T(index_val,8) == blade_row.gap_o
            G2 = T(index_val,1);
        end
    end
    
    % G3 i=0 k=nk
    XRRT = [g{blade_row.gap_h+1}.x(1,1,end) ; g{blade_row.gap_h+1}.r(1,1,end) ; g{blade_row.gap_h+1}.rt(1,1,end)];
    index_vals = kdtree_range(Tree,double([XRRT - tolerance , XRRT + tolerance]));

    for index_val = index_vals
        if T(index_val,8) == blade_row.gap_o
            G3 = T(index_val,1);
        end
    end
    
    % G4 i=ni k=ni
    XRRT = [g{blade_row.gap_h+1}.x(end,1,end) ; g{blade_row.gap_h+1}.r(end,1,end) ; g{blade_row.gap_h+1}.rt(end,1,end)];
    index_vals = kdtree_range(Tree,double([XRRT - tolerance , XRRT + tolerance]));

    for index_val = index_vals
        if T(index_val,8) == blade_row.gap_o
            G4 = T(index_val,1);
        end
    end
    
    % Seam of gap O-mesh - k coordinate on H-mesh
    XRRT = [g{blade_row.gap_o+1}.x(1,1,1) ; g{blade_row.gap_o+1}.r(1,1,1) ; g{blade_row.gap_o+1}.rt(1,1,1)];
    index_vals = kdtree_range(Tree,double([XRRT - tolerance , XRRT + tolerance]));

    for index_val = index_vals
        if T(index_val,8) == blade_row.gap_h
            G5 = T(index_val,3);
        end
    end
    
    % Patch gap H-mesh to gap O-mesh clockwise from i = 0
    Q0 = [0 1 ; 0 g{blade_row.gap_h+1}.attribute.nj ;...
        G5 g{blade_row.gap_h+1}.attribute.nk];
    dir0 = [-1 1 0];
    Q1 = [0 G3+1 ; 0 g{blade_row.gap_o+1}.attribute.nj ;...
        0 1];
    dir1 = [2 1 -1];
    [ g ] = ts_add_periodic_patch(blade_row.gap_h, blade_row.gap_o, Q0, Q1, dir0, dir1, g, 5);
    
    Q0 = [0 g{blade_row.gap_h+1}.attribute.ni ; 0 g{blade_row.gap_h+1}.attribute.nj ;...
        g{blade_row.gap_h+1}.attribute.nk-1 g{blade_row.gap_h+1}.attribute.nk];
    dir0 = [0 1 -1];
    Q1 = [G3 G4+1 ; 0 g{blade_row.gap_o+1}.attribute.nj ;...
        0 1];
    dir1 = [0 1 -1];
    [ g ] = ts_add_periodic_patch(blade_row.gap_h, blade_row.gap_o, Q0, Q1, dir0, dir1, g, 5);
    
    Q0 = [g{blade_row.gap_h+1}.attribute.ni-1 g{blade_row.gap_h+1}.attribute.ni ; 0 g{blade_row.gap_h+1}.attribute.nj ;...
        0 g{blade_row.gap_h+1}.attribute.nk];
    dir0 = [-1 1 3];
    Q1 = [G4 G2+1 ; 0 g{blade_row.gap_o+1}.attribute.nj ;...
        0 1];
    dir1 = [5 1 -1];
    [ g ] = ts_add_periodic_patch(blade_row.gap_h, blade_row.gap_o, Q0, Q1, dir0, dir1, g, 5);    
    
    Q0 = [0 g{blade_row.gap_h+1}.attribute.ni ; 0 g{blade_row.gap_h+1}.attribute.nj ;...
        0 1];
    dir0 = [3 1 -1];
    Q1 = [G2 G1+1 ; 0 g{blade_row.gap_o+1}.attribute.nj ;...
        0 1];
    dir1 = [3 1 -1];
    [ g ] = ts_add_periodic_patch(blade_row.gap_h, blade_row.gap_o, Q0, Q1, dir0, dir1, g, 5);
    
    Q0 = [0 1 ; 0 g{blade_row.gap_h+1}.attribute.nj ;...
        0 G5+1];
    dir0 = [-1 1 0];
    Q1 = [G1 g{blade_row.gap_o+1}.attribute.ni ; 0 g{blade_row.gap_o+1}.attribute.nj ;...
        0 1];
    dir1 = [2 1 -1];
    [ g ] = ts_add_periodic_patch(blade_row.gap_h, blade_row.gap_o, Q0, Q1, dir0, dir1, g, 5);
    
    % Seam of gap O-mesh
    Q0 = [0 1 ; 0 g{blade_row.gap_o+1}.attribute.nj ;...
        0 g{blade_row.gap_o+1}.attribute.nk];
    dir0 = [-1 1 2];
    Q1 = [g{blade_row.gap_o+1}.attribute.ni-1 g{blade_row.gap_o+1}.attribute.ni ; 0 g{blade_row.gap_o+1}.attribute.nj ;...
        0 g{blade_row.gap_o+1}.attribute.nk];
    dir1 = [-1 1 2];
    [ g ] = ts_add_periodic_patch(blade_row.gap_o, blade_row.gap_o, Q0, Q1, dir0, dir1, g, 5);
    
elseif isfield(blade_row,'gap_h') == 1 && isfield(blade_row,'gap_o') == 0
    % Case of no intermediate o mesh in tip gap
    
    % Find corner points of tip gap h mesh to blade o mesh
    % G1 i=0 k=0
    XRRT = [g{blade_row.gap_h+1}.x(1,1,1) ; g{blade_row.gap_h+1}.r(1,1,1) ; g{blade_row.gap_h+1}.rt(1,1,1)];
    index_vals = kdtree_range(Tree,double([XRRT - tolerance , XRRT + tolerance]));

    for index_val = index_vals
        if T(index_val,8) == blade_row.o
            G1 = T(index_val,1);
        end
    end

    % G2 i=ni k=0
    XRRT = [g{blade_row.gap_h+1}.x(end,1,1) ; g{blade_row.gap_h+1}.r(end,1,1) ; g{blade_row.gap_h+1}.rt(end,1,1)];
    index_vals = kdtree_range(Tree,double([XRRT - tolerance , XRRT + tolerance]));

    for index_val = index_vals
        if T(index_val,8) == blade_row.o
            G2 = T(index_val,1);
        end
    end
    
    % G3 i=0 k=nk
    XRRT = [g{blade_row.gap_h+1}.x(1,1,end) ; g{blade_row.gap_h+1}.r(1,1,end) ; g{blade_row.gap_h+1}.rt(1,1,end)];
    index_vals = kdtree_range(Tree,double([XRRT - tolerance , XRRT + tolerance]));

    for index_val = index_vals
        if T(index_val,8) == blade_row.o
            G3 = T(index_val,1);
        end
    end
    
    % G4 i=ni k=ni
    XRRT = [g{blade_row.gap_h+1}.x(end,1,end) ; g{blade_row.gap_h+1}.r(end,1,end) ; g{blade_row.gap_h+1}.rt(end,1,end)];
    index_vals = kdtree_range(Tree,double([XRRT - tolerance , XRRT + tolerance]));

    for index_val = index_vals
        if T(index_val,8) == blade_row.o
            G4 = T(index_val,1);
        end
    end
    
    % Seam of gap O-mesh - k coordinate on H-mesh
    XRRT = [g{blade_row.o+1}.x(1,end,1) ; g{blade_row.o+1}.r(1,end,1) ; g{blade_row.o+1}.rt(1,end,1)];
    index_vals = kdtree_range(Tree,double([XRRT - tolerance , XRRT + tolerance]));

    for index_val = index_vals
        if T(index_val,8) == blade_row.gap_h
            G5 = T(index_val,3);
        end
    end
    
    % Patch gap H-mesh to gap O-mesh clockwise from i = 0
    Q0 = [0 1 ; 0 g{blade_row.gap_h+1}.attribute.nj ;...
        G5 g{blade_row.gap_h+1}.attribute.nk];
    dir0 = [-1 1 0];
    Q1 = [0 G3+1 ; g{blade_row.o+1}.attribute.nj-g{blade_row.gap_h+1}.attribute.nj g{blade_row.o+1}.attribute.nj ;...
        0 1];
    dir1 = [2 1 -1];
    [ g ] = ts_add_periodic_patch(blade_row.gap_h, blade_row.o, Q0, Q1, dir0, dir1, g, 5);
    
    Q0 = [0 g{blade_row.gap_h+1}.attribute.ni ; 0 g{blade_row.gap_h+1}.attribute.nj ;...
        g{blade_row.gap_h+1}.attribute.nk-1 g{blade_row.gap_h+1}.attribute.nk];
    dir0 = [0 1 -1];
    Q1 = [G3 G4+1 ; g{blade_row.o+1}.attribute.nj-g{blade_row.gap_h+1}.attribute.nj g{blade_row.o+1}.attribute.nj ;...
        0 1];
    dir1 = [0 1 -1];
    [ g ] = ts_add_periodic_patch(blade_row.gap_h, blade_row.o, Q0, Q1, dir0, dir1, g, 5);
    
    Q0 = [g{blade_row.gap_h+1}.attribute.ni-1 g{blade_row.gap_h+1}.attribute.ni ; 0 g{blade_row.gap_h+1}.attribute.nj ;...
        0 g{blade_row.gap_h+1}.attribute.nk];
    dir0 = [-1 1 3];
    Q1 = [G4 G2+1 ; g{blade_row.o+1}.attribute.nj-g{blade_row.gap_h+1}.attribute.nj g{blade_row.o+1}.attribute.nj ;...
        0 1];
    dir1 = [5 1 -1];
    [ g ] = ts_add_periodic_patch(blade_row.gap_h, blade_row.o, Q0, Q1, dir0, dir1, g, 5);    
    
    Q0 = [0 g{blade_row.gap_h+1}.attribute.ni ; 0 g{blade_row.gap_h+1}.attribute.nj ;...
        0 1];
    dir0 = [3 1 -1];
    Q1 = [G2 G1+1 ; g{blade_row.o+1}.attribute.nj-g{blade_row.gap_h+1}.attribute.nj g{blade_row.o+1}.attribute.nj ;...
        0 1];
    dir1 = [3 1 -1];
    [ g ] = ts_add_periodic_patch(blade_row.gap_h, blade_row.o, Q0, Q1, dir0, dir1, g, 5);
    
    Q0 = [0 1 ; 0 g{blade_row.gap_h+1}.attribute.nj ;...
        0 G5+1];
    dir0 = [-1 1 0];
    Q1 = [G1 g{blade_row.o+1}.attribute.ni ; g{blade_row.o+1}.attribute.nj-g{blade_row.gap_h+1}.attribute.nj g{blade_row.o+1}.attribute.nj ;...
        0 1];
    dir1 = [2 1 -1];
    [ g ] = ts_add_periodic_patch(blade_row.gap_h, blade_row.o, Q0, Q1, dir0, dir1, g, 5);    
    
end
    
% Deal with outer inlet and exit blocks if they exist
if isfield(blade_row,'outer_inlet')
    % Across passage patch
    Q0 = [0 g{blade_row.outer_inlet+1}.attribute.ni ; 0 g{blade_row.outer_inlet+1}.attribute.nj ;...
        0 1];
    dir0 = [0 1 2];
    Q1 = [0 g{blade_row.outer_inlet+1}.attribute.ni ; 0 g{blade_row.outer_inlet+1}.attribute.nj ;...
        g{blade_row.outer_inlet+1}.attribute.nk-1 g{blade_row.outer_inlet+1}.attribute.nk];
    dir1 = [0 1 2];
    [ g ] = ts_add_periodic_patch(blade_row.outer_inlet, blade_row.outer_inlet, Q0, Q1, dir0, dir1, g, 5);
    
    % Patch to inner inlet
    Q0 = [g{blade_row.outer_inlet+1}.attribute.ni-1 g{blade_row.outer_inlet+1}.attribute.ni ; 0 g{blade_row.outer_inlet+1}.attribute.nj ;...
        0 g{blade_row.outer_inlet+1}.attribute.nk];
    dir0 = [0 1 2];
    Q1 = [0 1 ; 0 g{blade_row.inlet+1}.attribute.nj ;...
        0 g{blade_row.inlet+1}.attribute.nk];
    dir1 = [0 1 2];
    [ g ] = ts_add_periodic_patch(blade_row.outer_inlet, blade_row.inlet, Q0, Q1, dir0, dir1, g, 5);    
end

if isfield(blade_row,'outer_exit')
    % Across passage patch
    Q0 = [0 g{blade_row.outer_exit+1}.attribute.ni ; 0 g{blade_row.outer_exit+1}.attribute.nj ;...
        0 1];
    dir0 = [0 1 2];
    Q1 = [0 g{blade_row.outer_exit+1}.attribute.ni ; 0 g{blade_row.outer_exit+1}.attribute.nj ;...
        g{blade_row.outer_exit+1}.attribute.nk-1 g{blade_row.outer_exit+1}.attribute.nk];
    dir1 = [0 1 2];
    [ g ] = ts_add_periodic_patch(blade_row.outer_exit, blade_row.outer_exit, Q0, Q1, dir0, dir1, g, 5);
    
    % Patch to inner exit
    Q0 = [0 1 ; 0 g{blade_row.outer_exit+1}.attribute.nj ;...
        0 g{blade_row.outer_exit+1}.attribute.nk];
    dir0 = [0 1 2];
    Q1 = [g{blade_row.exit+1}.attribute.ni-1 g{blade_row.exit+1}.attribute.ni ; 0 g{blade_row.exit+1}.attribute.nj ;...
        0 g{blade_row.exit+1}.attribute.nk];
    dir1 = [0 1 2];
    [ g ] = ts_add_periodic_patch(blade_row.outer_exit, blade_row.exit, Q0, Q1, dir0, dir1, g, 5);    
end
    
    
end










