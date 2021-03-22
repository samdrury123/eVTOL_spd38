function g = ts_merge_radial(g,bids,scale)
% Match up radial distributions at a sliding plane

for n = 1:size(bids,1);
    r_1 = g{bids(n,1)+1}.r(1,:,:);
    x_1 = g{bids(n,1)+1}.x(1,:,:);
    r_3 = g{bids(n,2)+1}.r(end,:,:);
    x_3 = g{bids(n,2)+1}.x(end,:,:);

    r = scale(n) *mean(r_1,3) + (1-scale(n)) * mean(r_3,3);

    x_21 = mean(g{bids(n,1)+1}.x(end,:,:),3);
    r_21 = mean(g{bids(n,1)+1}.r(end,:,:),3);

    g{bids(n,1)+1}.r = (repmat(r,[g{bids(n,1)+1}.attribute.ni 1 g{bids(n,1)+1}.attribute.nk]) - ...
        repmat(r_1,[g{bids(n,1)+1}.attribute.ni 1 1])) .* ...
        (g{bids(n,1)+1}.x - repmat(x_1,[g{bids(n,1)+1}.attribute.ni 1 1])) ./ ...
        (repmat(x_21,[g{bids(n,1)+1}.attribute.ni 1 g{bids(n,1)+1}.attribute.nk]) - ...
        repmat(x_1,[g{bids(n,1)+1}.attribute.ni 1 1])) +...
        repmat(r_1,[g{bids(n,1)+1}.attribute.ni 1 1]); 

    dx = zeros(1,g{bids(n,1)+1}.attribute.nj,g{bids(n,1)+1}.attribute.nk);
    for k = 1:g{bids(n,1)+1}.attribute.nk
        dx(1,:,k) = interp1(r_21,x_21,g{bids(n,1)+1}.r(end,:,k),'pchip') -g{bids(n,1)+1}.x(end,:,k);
    end

    g{bids(n,1)+1}.x = g{bids(n,1)+1}.x + repmat(dx,[g{bids(n,1)+1}.attribute.ni 1 1]) .* ...
    (g{bids(n,1)+1}.x - repmat(x_1,[g{bids(n,1)+1}.attribute.ni 1 1])) ./ ...
        (repmat(x_21,[g{bids(n,1)+1}.attribute.ni 1 g{bids(n,1)+1}.attribute.nk]) - ...
        repmat(x_1,[g{bids(n,1)+1}.attribute.ni 1 1]));

    dx = zeros(1,g{bids(n,1)+1}.attribute.nj,g{bids(n,1)+1}.attribute.nk);
    for k = 1:g{bids(n,1)+1}.attribute.nk
        dx(1,:,k) = interp1(r_21,x_21,g{bids(n,1)+1}.r(end,:,k),'pchip') -g{bids(n,1)+1}.x(end,:,k);
    end

    g{bids(n,1)+1}.x = g{bids(n,1)+1}.x + repmat(dx,[g{bids(n,1)+1}.attribute.ni 1 1]) .* ...
    (g{bids(n,1)+1}.x - repmat(x_1,[g{bids(n,1)+1}.attribute.ni 1 1])) ./ ...
        (repmat(x_21,[g{bids(n,1)+1}.attribute.ni 1 g{bids(n,1)+1}.attribute.nk]) - ...
        repmat(x_1,[g{bids(n,1)+1}.attribute.ni 1 1]));



    x_22 = mean(g{bids(n,2)+1}.x(1,:,:),3);
    r_22 = mean(g{bids(n,2)+1}.r(1,:,:),3);

    g{bids(n,2)+1}.r = (repmat(r,[g{bids(n,2)+1}.attribute.ni 1 g{bids(n,2)+1}.attribute.nk]) - ...
        repmat(r_3,[g{bids(n,2)+1}.attribute.ni 1 1])) .* ...
        (g{bids(n,2)+1}.x - repmat(x_3,[g{bids(n,2)+1}.attribute.ni 1 1])) ./ ...
        (repmat(x_22,[g{bids(n,2)+1}.attribute.ni 1 g{bids(n,2)+1}.attribute.nk]) - ...
        repmat(x_3,[g{bids(n,2)+1}.attribute.ni 1 1])) +...
        repmat(r_3,[g{bids(n,2)+1}.attribute.ni 1 1]);

    dx = zeros(1,g{bids(n,2)+1}.attribute.nj,g{bids(n,2)+1}.attribute.nk);
    for k = 1:g{bids(n,2)+1}.attribute.nk
        dx(1,:,k) = interp1(r_22,x_22,g{bids(n,2)+1}.r(1,:,k),'pchip') -g{bids(n,2)+1}.x(1,:,k);
    end

    g{bids(n,2)+1}.x = g{bids(n,2)+1}.x + repmat(dx,[g{bids(n,2)+1}.attribute.ni 1 1]) .* ...
    (g{bids(n,2)+1}.x - repmat(x_3,[g{bids(n,2)+1}.attribute.ni 1 1])) ./ ...
        (repmat(x_21,[g{bids(n,2)+1}.attribute.ni 1 g{bids(n,2)+1}.attribute.nk]) - ...
        repmat(x_3,[g{bids(n,2)+1}.attribute.ni 1 1]));

    dx = zeros(1,g{bids(n,2)+1}.attribute.nj,g{bids(n,2)+1}.attribute.nk);
    for k = 1:g{bids(n,2)+1}.attribute.nk
        dx(1,:,k) = interp1(r_22,x_22,g{bids(n,2)+1}.r(1,:,k),'pchip') -g{bids(n,2)+1}.x(1,:,k);
    end

    g{bids(n,2)+1}.x = g{bids(n,2)+1}.x + repmat(dx,[g{bids(n,2)+1}.attribute.ni 1 1]) .* ...
    (g{bids(n,2)+1}.x - repmat(x_3,[g{bids(n,2)+1}.attribute.ni 1 1])) ./ ...
        (repmat(x_21,[g{bids(n,2)+1}.attribute.ni 1 g{bids(n,2)+1}.attribute.nk]) - ...
        repmat(x_3,[g{bids(n,2)+1}.attribute.ni 1 1]));


    x = 0.5 * (mean(g{bids(n,2)+1}.x(1,:,:),3) + mean(g{bids(n,1)+1}.x(end,:,:),3));
    r = 0.5 * (mean(g{bids(n,2)+1}.r(1,:,:),3) + mean(g{bids(n,1)+1}.r(end,:,:),3));

    g{bids(n,2)+1}.x(1,:,:) = repmat(x,[1 1 g{bids(n,2)+1}.attribute.nk]);
    g{bids(n,1)+1}.x(end,:,:) = repmat(x,[1 1 g{bids(n,1)+1}.attribute.nk]);
    g{bids(n,2)+1}.r(1,:,:) = repmat(r,[1 1 g{bids(n,2)+1}.attribute.nk]);
    g{bids(n,1)+1}.r(end,:,:) = repmat(r,[1 1 g{bids(n,1)+1}.attribute.nk]);
end


end