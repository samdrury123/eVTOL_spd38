function g = ts_renumber_patches(g,bid,pids)
% Change the order or delete patches from a TURBOSTREAM mesh

% Get all connecting patch data
nxbids = zeros(size(g{bid+1}.patch)); nxpids = zeros(size(g{bid+1}.patch)); 
for pid = ts_get_pids(g,bid)
    nxbids(pid+1) = g{bid+1}.patch{pid+1}.attribute.nxbid;
    nxpids(pid+1) = g{bid+1}.patch{pid+1}.attribute.nxpid;
end
    
% Renumber nxpids on other blocks
% nxbids = nxbids(pids+1); nxpids = nxpids(pids+1);
for n = 1:length(nxbids)
    g{nxbids(n)+1}.patch{nxpids(n)+1}.attribute.nxpid = find(pids == (n-1)) - 1;
end

% Reorder all patches
g{bid+1}.patch = g{bid+1}.patch(pids+1);

% Renumber pids on all patches
for pid = 0:length(g{bid+1}.patch)-1
    g{bid+1}.patch{pid+1}.attribute.pid = pid;
end

% % Move old sliding planes to the end of the list
% for n = 1:length(bid_slide)
%     bid = bid_slide(n);
%     for pid = pid_slide(n)+1:g{bid+1}.attribute.np-1
%         g{bid+1}.patch{pid+1}.attribute.pid = pid-1;
%         nxbid = g{bid+1}.patch{pid+1}.attribute.nxbid; 
%         nxpid = g{bid+1}.patch{pid+1}.attribute.nxpid;
%         g{nxbid+1}.patch{nxpid+1}.attribute.nxpid = pid-1;
%     end
%     g{bid+1}.patch = [g{bid+1}.patch(1:pid_slide(n)) ; g{bid+1}.patch(pid_slide(n)+2:end) ; ...
%         g{bid+1}.patch(pid_slide(n)+1)];
%     g{bid+1}.patch{end}.attribute.pid = length(g{bid+1}.patch)-1;
%     nxbid = g{bid+1}.patch{end}.attribute.nxbid;
%     g{bid+1}.patch{end}.attribute.nxpid = length(g{nxbid+1}.patch)-1;
% end

end