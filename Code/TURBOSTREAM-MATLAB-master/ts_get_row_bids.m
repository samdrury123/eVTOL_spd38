function [B,bid_os] = ts_get_row_bids(g,xr_hub,xr_cas)
% Return a cell array with bids in each row determined by mixing planes

% Default hub and casing radii if not specified
if exist('xr_hub','var') == 0 || isempty(xr_hub) == 1
    xr_hub = [-100 min(g{1}.r(:)) ; 100 min(g{1}.r(:))];
end
if exist('xr_cas','var') == 0 || isempty(xr_cas) == 1
    xr_cas = [-100 max(g{1}.r(:)) ; 100 max(g{1}.r(:))];
end

% Get o-mesh bids
bid_os = ts_find_os(g,[mean(xr_hub(:,2)) mean(xr_cas(:,2))]);
B = cell(length(bid_os),1);

% Loop over all blocks and patches to find blocks connected by periodic patches
for n = 1:length(bid_os)
    bids = bid_os(n);
    old_length = 0;
    while length(bids) ~= old_length
        old_length = length(bids);
        for bid = bids
            for pid = ts_get_pids(g,bid)
                if g{bid+1}.patch{pid+1}.attribute.kind == 5
                    bids = [bids g{bid+1}.patch{pid+1}.attribute.nxbid];
                end
            end
        end
        bids = unique(bids);
    end
    B{n} = bids;
end

% Return only unique cell arrays for each row
bids_used = zeros(length(g),1); B_new = cell(0,1);
for n = 1:length(B)
    if all(bids_used(B{n}+1) == 0)
        B_new = [B_new ; B(n)];
        bids_used(B{n}+1) = 1;
    end
end
B = B_new;

end