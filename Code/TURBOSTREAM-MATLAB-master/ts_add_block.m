function [ Data ] = ts_add_block(b, Data)
% Add a block to a matlab data_set

% See if block number is already specified
if b.attribute.bid == -1
    try
        b.attribute.bid = Data{end}.attribute.bid + 1;
    catch ME 
        b.attribute.bid = 0;
    end
end

b.attribute.np = length(ts_get_pids(b));

bids = ts_get_bids(Data);
block_location = bids == b.attribute.bid;

if sum(block_location) > 0
    Data{block_location} = b;
else
    Data = [Data ; b];
end

end

