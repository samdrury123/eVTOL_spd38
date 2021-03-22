function g = ts_remove_patch(g,bid,pid)
% Delete a patch from TS array and update the rest of the mesh

% Delete the patch
g{bid+1}.patch(pid+1) = [];

% Update the patches on the rest of the mesh
for pid = pid:length(g{bid+1}.patch)-1
    
    % Update the patch id
    g{bid+1}.patch{pid+1}.attribute.pid = pid;
    
    % Update attributes on connecting blocks
    if sum(g{bid+1}.patch{pid+1}.attribute.kind == [2 3 5]) > 0
        nxpid = g{bid+1}.patch{pid+1}.attribute.nxpid;
        nxbid = g{bid+1}.patch{pid+1}.attribute.nxbid;
        if length(g{nxbid+1}.patch) > nxpid
            g{nxbid+1}.patch{nxpid+1}.attribute.nxpid = pid;
        end
    end
    
end

% Update the number of patches
g{bid+1}.attribute.np = length(g{bid+1}.patch);

end