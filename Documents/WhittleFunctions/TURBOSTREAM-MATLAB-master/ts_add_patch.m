function [ Data ] = ts_add_patch(bid, p, Data)
% Add a patch to a block in the data set

b = ts_get_block(bid, Data);

% Check if any previous patches exist and increment pid
if p.attribute.pid == -1
    try 
        p.attribute.pid = b.patch{end}.attribute.pid + 1;
    catch ME
        p.attribute.pid = 0;
    end
end

pids = ts_get_pids(b);
patch_location = pids == p.attribute.pid;

if sum(patch_location) > 0
    b.patch{patch_location} = p;
else
    b.patch = [b.patch ; p];
    b.attribute.np = b.attribute.np + 1;
end

Data = ts_add_block(b, Data);

end

