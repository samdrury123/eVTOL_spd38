function g = ts_set_vs(g,filename)
% TS_SET_VS  Read application and block variables from file and apply to grid

% Load in file
load(filename)

% Set application variables
avnames = fieldnames(av);
for v = 1:length(avnames)
    g{1}.av.(avnames{v}) = av.(avnames{v});
end

% Set block variables
for bid = ts_get_bids(g)
    bvnames = fieldnames(bv{bid+1});
    for v = 1:length(bvnames)
        g{bid+1}.bv.(bvnames{v}) = bv{bid+1}.(bvnames{v});
    end
end

end