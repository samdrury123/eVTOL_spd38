function [] = ts_check_mwall(g)
% Checks that the mesh has been properly patched by looking at the
% locations of mwall = 0

for bid = ts_get_bids(g)
    
    fprintf('\nBlock %i \n\n', bid)
    fprintf('\t ni = %i; nj = %i; nk = %i; \n', g{bid+1}.attribute.ni, g{bid+1}.attribute.nj, g{bid+1}.attribute.nk)
    
    mwall_flag = sum(sum(sum(g{bid+1}.mwall(2:end-1,2:end-1, 2:end-1) == 0)));
    
    fprintf('\t Number of wall nodes not on a wall = %i \n\n', mwall_flag)
    
end

end