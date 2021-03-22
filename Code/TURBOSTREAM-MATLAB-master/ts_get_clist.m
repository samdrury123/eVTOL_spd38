function [ cl ] = ts_get_clist(bid, clid, Data)
% Pull a patch from a matlab data cell array

b = ts_get_block(bid, Data);

% Try indexed patches
cl_temp = b.clist{clid+1};

% If it fails loop over all patches
if cl_temp.attribute.clid ~= clid
    for n = 1:length(b.clist)
        cl_temp = b.clist{n};
        if cl_temp.attribute.clid == clid
            cl = cl_temp;
            break
        end
    end
else
    cl = cl_temp;
end


end