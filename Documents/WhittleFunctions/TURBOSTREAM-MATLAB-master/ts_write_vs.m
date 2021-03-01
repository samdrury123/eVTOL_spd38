function [av,bv] = ts_write_vs(g,filename)
% TS_WRITE_VS  Write application and block variables to file

% Get all application variables names and values
avnames = fieldnames(g{1}.av);
for v = 1:length(avnames)
    av.(avnames{v}) = g{1}.av.(avnames{v});
end

% Get all block application variables names and values
bv = cell(size(g));
for bid = ts_get_bids(g)
    bvnames = fieldnames(g{bid+1}.bv);
    for v = 1:length(bvnames); bv{bid+1}.(bvnames{v}) = g{bid+1}.bv.(bvnames{v}); end;
end

% Save variables
save(filename,'av','bv');


end

