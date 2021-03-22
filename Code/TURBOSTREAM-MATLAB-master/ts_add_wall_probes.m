function g = ts_add_wall_probes(g,bids)
% Add probe patches everywhere mwall = 0

if exist('bids','var') == 0
    bids = ts_get_bids(g);
end

for bid = bids
    [I,J,K] = ind2sub(size(g{bid+1}.mwall), find(g{bid+1}.mwall == 0));
    I = I-1; J = J-1; K = K-1;

    % Plot i = 0 surface
    i0 = 0; i1 = 0;
    n = find(I == 0);
    if ~isempty(n)
        j0 = min(J(n)); j1 = max(J(n));
        k0 = min(K(n)); k1 = max(K(n));
        a = (k1 - k0 + 1); b = (j1 - j0 + 1);
        if a * b == length(n) && a > 1 && b > 1 
            g = add_probe(g,bid,i0,i1,j0,j1,k0,k1);
        end
    end
  
    % Plot i = ni surface
    ni = g{bid+1}.attribute.ni-1;
    i0 = ni; i1 = ni;
    n = find(I == ni);
    if ~isempty(n)
        j0 = min(J(n)); j1 = max(J(n));
        k0 = min(K(n)); k1 = max(K(n));
        a = (k1 - k0 + 1); b = (j1 - j0 + 1);
        if a * b == length(n) && a > 1 && b > 1
            g = add_probe(g,bid,i0,i1,j0,j1,k0,k1);
        end
    end

    % Plot j = 0 surface
    j0 = 0; j1 = 0;
    n = find(J == 0);
    if ~isempty(n)
        i0 = min(I(n)); i1 = max(I(n));
        k0 = min(K(n)); k1 = max(K(n));
        a = (k1 - k0 + 1); b = (i1 - i0 + 1);
        if a * b == length(n) && a > 1 && b > 1
            g = add_probe(g,bid,i0,i1,j0,j1,k0,k1);
        end
    end
    % Plot j = nj surface
    nj = g{bid+1}.attribute.nj-1;
    j0 = nj; j1 = nj;
    n = find(J == nj);
    if ~isempty(n)
        i0 = min(I(n)); i1 = max(I(n));
        k0 = min(K(n)); k1 = max(K(n));
        a = (k1 - k0 + 1); b = (i1 - i0 + 1);
        if a * b == length(n) && a > 1 && b > 1
            g = add_probe(g,bid,i0,i1,j0,j1,k0,k1);
        end
    end

    % Plot k = 0 surface
    k0 = 0; k1 = 0;
    n = find(K == 0);
    if ~isempty(n)
        % Allow for tip gaps
        j = unique(J(n)); 
        if max(diff(j)) > 1
            j1 = j(end-1);
            n = find(K == 0 & J~= j(end));
        else
            j1 = max(J(n));
        end
        i0 = min(I(n)); i1 = max(I(n));
        j0 = min(J(n));
        a = (j1 - j0 + 1); b = (i1 - i0 + 1);
        if a * b == length(n) && a > 1 && b > 1
            g = add_probe(g,bid,i0,i1,j0,j1,k0,k1);
        end
    end
    % Plot k = nk surface
    nk = g{bid+1}.attribute.nk-1;
    k0 = nk; k1 = nk;
    n = find(K == nk);
    if ~isempty(n)
        i0 = min(I(n)); i1 = max(I(n));
        j0 = min(J(n)); j1 = max(J(n));
        a = (j1 - j0 + 1); b = (i1 - i0 + 1);
        if a * b == length(n) && a > 1 && b > 1
            g = add_probe(g,bid,i0,i1,j0,j1,k0,k1);
        end
    end    
end

end

function g = add_probe(g,bid,i0,i1,j0,j1,k0,k1)

g{bid+1}.attribute.np = g{bid+1}.attribute.np + 1;

patch.attribute.bid = bid;
patch.attribute.pid = length(g{bid+1}.patch);
patch.attribute.kind = 8;
patch.attribute.ist = i0; patch.attribute.ien = i1+1;
patch.attribute.jst = j0; patch.attribute.jen = j1+1;
patch.attribute.kst = k0; patch.attribute.ken = k1+1;

g{bid+1}.patch{end+1} = patch;

end