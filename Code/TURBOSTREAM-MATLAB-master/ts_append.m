function g = ts_append(g,g_end)
% TS_APPEND  Join two TURBOSTREAM cell arrays together

% Starting bid of new blocks
bid_st = length(g);

% Update block numbers in attributes and patch attributes
for bid = 0:length(g_end)-1
    g_end{bid+1}.attribute.bid = g_end{bid+1}.attribute.bid + bid_st;
    for pid = 0:length(g_end{bid+1}.patch)-1
        g_end{bid+1}.patch{pid+1}.attribute.bid = ...
            g_end{bid+1}.patch{pid+1}.attribute.bid + bid_st;
        g_end{bid+1}.patch{pid+1}.attribute.nxbid = ...
            g_end{bid+1}.patch{pid+1}.attribute.nxbid + bid_st;
    end
end

% Append end blocks to main mesh
g = [g ; g_end];

end

