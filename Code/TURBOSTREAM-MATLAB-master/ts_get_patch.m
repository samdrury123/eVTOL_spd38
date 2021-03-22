function [ p ] = ts_get_patch(bid, pid, Data)
% Pull a patch from a matlab data cell array

b = ts_get_block(bid, Data);

% Try indexed patches
p_temp = b.patch{pid+1};

% If it fails loop over all patches
if p_temp.attribute.pid ~= pid
    for n = 1:length(b.patch)
        p_temp = b.patch{n};
        if p_temp.attribute.pid == pid
            p = p_temp;
            break
        end
    end
else
    p = p_temp;
end


end