function g_surf = ts_get_surface(g,bids)
% TS_GET_SURFACE  Extract cell array of all surfaces in a TURBOSTREAM mesh

% Default to all blocks
if exist('bids','var') == 0 || isempty(bids) == 1
    bids = ts_get_bids(g);
end

% Initialise surface array
g_surf = cell(0,1);

% Loop over all blocks
for bid = bids

    % Indices of wall patches
    [I,J,K] = ind2sub(size(g{bid+1}.mwall), find(g{bid+1}.mwall == 0));
    I = I-1; J = J-1; K = K-1;

    % Extract i = 0 surface
    i0 = 0; i1 = 1; n = find(I == 0);
    if ~isempty(n)
        j0 = min(J(n)); j1 = max(J(n));
        k0 = min(K(n)); k1 = max(K(n));
        a = (k1 - k0 + 1); b = (j1 - j0 + 1);
        if a * b == length(n) && a > 1 && b > 1 
            g_surf{end+1} = get_surface(g,bid,i0,i1,j0,j1+1,k0,k1+1); 
            g_surf{end}.attribute.bid = length(g_surf) - 1;
        end
    end
  
    % Extract i = ni surface
    ni = g{bid+1}.attribute.ni-1;
    i0 = ni-1; i1 = ni; n = find(I == ni);
    if ~isempty(n)
        j0 = min(J(n)); j1 = max(J(n));
        k0 = min(K(n)); k1 = max(K(n));
        a = (k1 - k0 + 1); b = (j1 - j0 + 1);
        if a * b == length(n) && a > 1 && b > 1
            g_surf{end+1} = get_surface(g,bid,i0,i1,j0,j1+1,k0,k1+1); 
            g_surf{end}.attribute.bid = length(g_surf) - 1;
        end
    end

    % Extract j = 0 surface
    j0 = 0; j1 = 1; n = find(J == 0);
    if ~isempty(n)
        i0 = min(I(n)); i1 = max(I(n));
        k0 = min(K(n)); k1 = max(K(n));
        a = (k1 - k0 + 1); b = (i1 - i0 + 1);
        if a * b == length(n) && a > 1 && b > 1
            g_surf{end+1} = get_surface(g,bid,i0,i1+1,j0,j1,k0,k1+1); 
            g_surf{end}.attribute.bid = length(g_surf) - 1;
        end
    end
    
    % Extract j = nj surface
    nj = g{bid+1}.attribute.nj-1;
    j0 = nj-1; j1 = nj; n = find(J == nj);
    if ~isempty(n)
        i0 = min(I(n)); i1 = max(I(n));
        k0 = min(K(n)); k1 = max(K(n));
        a = (k1 - k0 + 1); b = (i1 - i0 + 1);
        if a * b == length(n) && a > 1 && b > 1
            g_surf{end+1} = get_surface(g,bid,i0,i1+1,j0,j1,k0,k1+1); 
            g_surf{end}.attribute.bid = length(g_surf) - 1;
        end
    end

    % Extract k = 0 surface
    k0 = 0; k1 = 1; n = find(K == 0);
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
            g_surf{end+1} = get_surface(g,bid,i0,i1+1,j0,j1+1,k0,k1); 
            g_surf{end}.attribute.bid = length(g_surf) - 1;
        end
    end
    
    % Extract k = nk surface
    nk = g{bid+1}.attribute.nk-1;
    k0 = nk-1; k1 = nk; n = find(K == nk);
    if ~isempty(n)
        i0 = min(I(n)); i1 = max(I(n));
        j0 = min(J(n)); j1 = max(J(n));
        a = (j1 - j0 + 1); b = (i1 - i0 + 1);
        if a * b == length(n) && a > 1 && b > 1
            g_surf{end+1} = get_surface(g,bid,i0,i1+1,j0,j1+1,k0,k1); 
            g_surf{end}.attribute.bid = length(g_surf) - 1;
        end
    end    

end


end
    
function c = get_surface(g,bid,i0,i1,j0,j1,k0,k1)

% Properties to copy
varnames = {'x' 'r' 'rt' 'ro' 'rovx' 'rovr' 'rovr' 'rorvt' 'roe'};

% Collect surface properties and attributes into cut
for v = 1:length(varnames)
    c.(varnames{v}) = squeeze(g{bid+1}.(varnames{v})(i0+1:i1,j0+1:j1,k0+1:k1));
end
[c.attribute.ni,c.attribute.nj,c.attribute.nk] = size(c.x); 

% Record wall rotation for i walls
if i1 - i0 == 1
    if i0 == 0; c.bv.rpm = g{bid+1}.bv.rpmi1; else; c.bv.rpm = g{bid+1}.bv.rpmi2; end;
end

% Record wall rotation for j walls
if j1 - j0 == 1
    if j0 == 0; c.bv.rpm = g{bid+1}.bv.rpmj1; else; c.bv.rpm = g{bid+1}.bv.rpmj2; end;
end

% Record wall rotation for k walls
if k1 - k0 == 1
    if k0 == 0; c.bv.rpm = g{bid+1}.bv.rpmk1; else; c.bv.rpm = g{bid+1}.bv.rpmk2; end;
end

end