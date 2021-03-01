function g = ts_renumber_blocks(g,bids)
% Extract a list of blocks from a data structure and renumber blocks and
% patches

% Delete all patches associated with blocks that will not exist in new mesh
for bid = bids
    pid_remove = [];
    for pid = 0:length(g{bid+1}.patch) - 1
        kind = g{bid+1}.patch{pid+1}.attribute.kind;
        if kind == 5 || kind == 15 || kind == 2
            nxbid = g{bid+1}.patch{pid+1}.attribute.nxbid;
            
            if sum(nxbid == bids) == 0
                % Renumber all other patches below it
                if pid ~= length(g{bid+1}.patch) - 1
                    for pid2 = pid+1:length(g{bid+1}.patch) - 1
                        g{bid+1}.patch{pid2+1}.attribute.pid = pid2 - 1;
                        kind = g{bid+1}.patch{pid2+1}.attribute.kind;
                        if kind == 5 || kind == 15 || kind == 2
                            nxbid = g{bid+1}.patch{pid2+1}.attribute.nxbid;
                            nxpid = g{bid+1}.patch{pid2+1}.attribute.nxpid;
                            
                            g{nxbid+1}.patch{nxpid+1}.attribute.nxpid = ...
                                g{nxbid+1}.patch{nxpid+1}.attribute.nxpid - 1;
                        end
                    end
                end
                
                g{bid+1}.attribute.np = g{bid+1}.attribute.np -1;
                
                pid_remove = [pid_remove pid];
            end
        end
    end
    % Delete the patches
    g{bid+1}.patch(pid_remove+1) = [];
end

% Extract blocks of interest
if exist('bids','var')
    g_new = cell(length(bids),1);
    for new_bid = 0:length(bids)-1
        bid = bids(new_bid+1);
        g_new{new_bid+1} = g{bid+1};
    end
    g = g_new;
else
    bids = ts_get_bids(g);
end

% Create bid reference list
ref_bids = -ones(max(bids)+1,1);
for bid = 0:length(g)-1
    ref_bids(bids(bid+1)+1) = bid;
end

% Renumber block and patch attributes
for bid = 0:length(g)-1
    g{bid+1}.attribute.bid = bid;
    for pid = 0:length(g{bid+1}.patch) - 1
        g{bid+1}.patch{pid+1}.attribute.bid = bid;
        kind = g{bid+1}.patch{pid+1}.attribute.kind;
        if kind == 5 || kind == 15 || kind == 2
            nxbid = g{bid+1}.patch{pid+1}.attribute.nxbid;
            if sum(nxbid == bids) > 0
                g{bid+1}.patch{pid+1}.attribute.nxbid = ref_bids(nxbid+1);
            end
        end
    end
end

% Catch any misnumbered pids
for bid = 0:length(g)-1
    for pid = 0:length(g{bid+1}.patch)-1
        g{bid+1}.patch{pid+1}.attribute.pid = pid;
    end
end

end