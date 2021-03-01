function normal = ts_calc_normals(g,coincident)
% Function to calculate normals to every node in given surfaces

normal = cell(size(g)); coinc_normal = cell(size(g)); type = cell(size(g));
I_norm = cell(size(g));  J_norm = cell(size(g));  K_norm = cell(size(g)); 
for bid = ts_get_bids(g)
    normal{bid+1} = nan([size(g{bid+1}.x) 3 12]);
    coinc_normal{bid+1} = nan([size(g{bid+1}.x) 3 12]);
    type{bid+1} = zeros([size(g{bid+1}.x) 3]);
end

% Calculate face normals to all faces in all cells
for bid = ts_get_bids(g)
    xrrt = cat(4,g{bid+1}.x,g{bid+1}.r,g{bid+1}.rt);
    
    % Cross diagonal vectors on all faces
    ai = xrrt(:,2:end,2:end,:) - xrrt(:,1:end-1,1:end-1,:);
    bi = xrrt(:,2:end,1:end-1,:) - xrrt(:,1:end-1,2:end,:);
    ni = cross(bi,ai);
    I_norm{bid+1} = ni ./ repmat(((ni(:,:,:,1).^2 + ni(:,:,:,2).^2 + ni(:,:,:,3).^2).^0.5),[1 1 1 3]);
    
    aj = xrrt(2:end,:,1:end-1,:) - xrrt(1:end-1,:,2:end,:);
    bj = xrrt(2:end,:,2:end,:) - xrrt(1:end-1,:,1:end-1,:);
    nj = cross(bj,aj);
    J_norm{bid+1} = nj ./ repmat(((nj(:,:,:,1).^2 + nj(:,:,:,2).^2 + nj(:,:,:,3).^2).^0.5),[1 1 1 3]);
    
    ak = xrrt(1:end-1,2:end,:,:) - xrrt(2:end,1:end-1,:,:);
    bk = xrrt(2:end,2:end,:,:) - xrrt(1:end-1,1:end-1,:,:);
    nk = cross(bk,ak);
    K_norm{bid+1} = nk ./ repmat(((nk(:,:,:,1).^2 + nk(:,:,:,2).^2 + nk(:,:,:,3).^2).^0.5),[1 1 1 3]);
end

% Assign correct face normals on edges of blocks
for bid = ts_get_bids(g)
    [I,J,K] = ind2sub(size(g{bid+1}.mwall), find(g{bid+1}.mwall == 0));
    I = I-1; J = J-1; K = K-1;
    
    % I = 0 face
    n = find(I == 0);
    if ~isempty(n)
        j0 = min(J(n)); j1 = max(J(n));
        k0 = min(K(n)); k1 = max(K(n));
        a = (k1 - k0 + 1); b = (j1 - j0 + 1);
        if a * b == length(n) && a > 1 && b > 1 
            normal{bid+1}(1,j0+1:j1,k0+1:k1,:,1)     =  I_norm{bid+1}(1,j0+1:j1,k0+1:k1,:);
            normal{bid+1}(1,j0+2:j1+1,k0+1:k1,:,2)   =  I_norm{bid+1}(1,j0+1:j1,k0+1:k1,:);
            normal{bid+1}(1,j0+2:j1+1,k0+2:k1+1,:,3) =  I_norm{bid+1}(1,j0+1:j1,k0+1:k1,:);
            normal{bid+1}(1,j0+1:j1,k0+2:k1+1,:,4)   =  I_norm{bid+1}(1,j0+1:j1,k0+1:k1,:);
        end
    end
    
    % I = ni face
    ni = g{bid+1}.attribute.ni-1;
    n = find(I == ni);
    if ~isempty(n)
        j0 = min(J(n)); j1 = max(J(n));
        k0 = min(K(n)); k1 = max(K(n));
        a = (k1 - k0 + 1); b = (j1 - j0 + 1);
        if a * b == length(n) && a > 1 && b > 1
            normal{bid+1}(end,j0+1:j1,k0+1:k1,:,1)     =  - I_norm{bid+1}(end,j0+1:j1,k0+1:k1,:);
            normal{bid+1}(end,j0+2:j1+1,k0+1:k1,:,2)   =  - I_norm{bid+1}(end,j0+1:j1,k0+1:k1,:);
            normal{bid+1}(end,j0+2:j1+1,k0+2:k1+1,:,3) =  - I_norm{bid+1}(end,j0+1:j1,k0+1:k1,:);
            normal{bid+1}(end,j0+1:j1,k0+2:k1+1,:,4)   =  - I_norm{bid+1}(end,j0+1:j1,k0+1:k1,:);
        end
    end
    
    % J = 0 face
    n = find(J == 0);
    if ~isempty(n)
        i0 = min(I(n)); i1 = max(I(n));
        k0 = min(K(n)); k1 = max(K(n));
        a = (k1 - k0 + 1); b = (i1 - i0 + 1);
        if a * b == length(n) && a > 1 && b > 1 
            normal{bid+1}(i0+1:i1,1,k0+1:k1,:,5)     =  J_norm{bid+1}(i0+1:i1,1,k0+1:k1,:);
            normal{bid+1}(i0+2:i1+1,1,k0+1:k1,:,6)   =  J_norm{bid+1}(i0+1:i1,1,k0+1:k1,:);
            normal{bid+1}(i0+2:i1+1,1,k0+2:k1+1,:,7) =  J_norm{bid+1}(i0+1:i1,1,k0+1:k1,:);
            normal{bid+1}(i0+1:i1,1,k0+2:k1+1,:,8)   =  J_norm{bid+1}(i0+1:i1,1,k0+1:k1,:);
        end
    end
    
    % J = nj face
    nj = g{bid+1}.attribute.nj-1;
    n = find(J == nj);
    if ~isempty(n)
        i0 = min(I(n)); i1 = max(I(n));
        k0 = min(K(n)); k1 = max(K(n));
        a = (k1 - k0 + 1); b = (i1 - i0 + 1);
        if a * b == length(n) && a > 1 && b > 1
            normal{bid+1}(i0+1:i1,end,k0+1:k1,:,5)     =  - J_norm{bid+1}(i0+1:i1,end,k0+1:k1,:);
            normal{bid+1}(i0+2:i1+1,end,k0+1:k1,:,6)   =  - J_norm{bid+1}(i0+1:i1,end,k0+1:k1,:);
            normal{bid+1}(i0+2:i1+1,end,k0+2:k1+1,:,7) =  - J_norm{bid+1}(i0+1:i1,end,k0+1:k1,:);
            normal{bid+1}(i0+1:i1,end,k0+2:k1+1,:,8)   =  - J_norm{bid+1}(i0+1:i1,end,k0+1:k1,:);
        end
    end
    
    % K = 0 face
    n = find(K == 0);
    if ~isempty(n)
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
            normal{bid+1}(i0+1:i1,j0+1:j1,1,:,9)      =  K_norm{bid+1}(i0+1:i1,j0+1:j1,1,:);
            normal{bid+1}(i0+1:i1,j0+2:j1+1,1,:,10)   =  K_norm{bid+1}(i0+1:i1,j0+1:j1,1,:);
            normal{bid+1}(i0+2:i1+1,j0+2:j1+1,1,:,11) =  K_norm{bid+1}(i0+1:i1,j0+1:j1,1,:);
            normal{bid+1}(i0+2:i1+1,j0+1:j1,1,:,12)   =  K_norm{bid+1}(i0+1:i1,j0+1:j1,1,:);
        end
    end
    
    % K = nk face
    nk = g{bid+1}.attribute.nk-1;
    n = find(K == nk);
    if ~isempty(n)
        i0 = min(I(n)); i1 = max(I(n));
        j0 = min(J(n)); j1 = max(J(n));
        a = (j1 - j0 + 1); b = (i1 - i0 + 1);
        if a * b == length(n) && a > 1 && b > 1
            normal{bid+1}(i0+1:i1,j0+1:j1,end,:,9)      =  - K_norm{bid+1}(i0+1:i1,j0+1:j1,end,:);
            normal{bid+1}(i0+1:i1,j0+2:j1+1,end,:,10)   =  - K_norm{bid+1}(i0+1:i1,j0+1:j1,end,:);
            normal{bid+1}(i0+2:i1+1,j0+2:j1+1,end,:,11) =  - K_norm{bid+1}(i0+1:i1,j0+1:j1,end,:);
            normal{bid+1}(i0+2:i1+1,j0+1:j1,end,:,12)   =  - K_norm{bid+1}(i0+1:i1,j0+1:j1,end,:);
        end
    end    
end
    
% Concatenate data from coincident points to get correct normals at block
% boundaries
for bid = ts_get_bids(g)
    [I,J,K] = ind2sub(size(g{bid+1}.mwall), find(g{bid+1}.mwall == 0));
    I = I-1; J = J-1; K = K-1;
    
    % Loop over all wall nodes
    for p = 1:length(I)
        % Get coincident points
        v = coincident{bid+1}{I(p)+1,J(p)+1,K(p)+1};
        q = 1;
        if isempty(v) == 0
            for m = 1:size(v,1)
                % Pull out normals only
                o = isnan(normal{v(m,1)+1}(v(m,2)+1,v(m,3)+1,v(m,4)+1,1,:)) == 0;
                v_normals = normal{v(m,1)+1}(v(m,2)+1,v(m,3)+1,v(m,4)+1,:,o);
                % Asign into new normal cell array
                coinc_normal{bid+1}(I(p)+1,J(p)+1,K(p)+1,:,q:(q + size(v_normals,5) - 1)) = v_normals;
                q = q + size(v_normals,5);
            end
        end
    end
end

% Update old normal array with new coincident normals
for bid = ts_get_bids(g)
    o = isnan(coinc_normal{bid+1}(:,:,:,1,1)) == 0;
    o = repmat(o,[1 1 1 3 12]);
    normal{bid+1}(o) = coinc_normal{bid+1}(o);
end
    
% Calculate single node normal as an average of those face normals stored
for bid = ts_get_bids(g)
    normal{bid+1}(isnan(normal{bid+1})) = 0;
    
    sum_normal = sum(normal{bid+1},5);
    mag_normal = (sum_normal(:,:,:,1).^2 + sum_normal(:,:,:,2).^2 + sum_normal(:,:,:,3).^2).^0.5;
    normal{bid+1} = sum_normal ./ repmat(mag_normal, [1 1 1 3]);
end


end