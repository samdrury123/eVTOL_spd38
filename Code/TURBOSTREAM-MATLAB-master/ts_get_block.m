function [ b ] = ts_get_block(bid, Data)
% Pull a block from a matlab data cell array

% Try indexed blocks
b_temp = Data{bid+1};

% If it fails loop over all blocks
if b_temp.attribute.bid ~= bid
    for n = 1:length(Data)
        b_temp = Data{n};
        if b_temp.attribute.bid == bid
            b = b_temp;
            break
        end
    end
else
    b = b_temp;
end


end