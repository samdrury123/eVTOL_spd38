function g = ts_set_IJK(g)
% Scans through and sets indexing of arrays to ensure that each block is a
% right handed set

bid_flip = [];
for bid = ts_get_bids(g)
    QI = [g{bid+1}.x(2,1,1)-g{bid+1}.x(1,1,1) ; 
        g{bid+1}.r(2,1,1)-g{bid+1}.r(1,1,1) ;
        g{bid+1}.rt(2,1,1)-g{bid+1}.rt(1,1,1)];
    QJ = [g{bid+1}.x(1,2,1)-g{bid+1}.x(1,1,1) ; 
        g{bid+1}.r(1,2,1)-g{bid+1}.r(1,1,1) ;
        g{bid+1}.rt(1,2,1)-g{bid+1}.rt(1,1,1)];
    QK = [g{bid+1}.x(1,1,2)-g{bid+1}.x(1,1,1) ; 
        g{bid+1}.r(1,1,2)-g{bid+1}.r(1,1,1) ;
        g{bid+1}.rt(1,1,2)-g{bid+1}.rt(1,1,1)];
    
    % Check for righthanded set and flip i-direction
    if dot(QK,cross(QI,QJ)) < 0
%         bid_flip = [bid_flip ; bid];
        g{bid+1}.x = flipdim(g{bid+1}.x,1);
        g{bid+1}.r = flipdim(g{bid+1}.r,1);
        g{bid+1}.rt = flipdim(g{bid+1}.rt,1);
        
        % Change patch starts and ends in i direction
        ni = g{bid+1}.attribute.ni;
        for pid = ts_get_pids(g,bid)
            ist = g{bid+1}.patch{pid+1}.attribute.ist;
            ien = g{bid+1}.patch{pid+1}.attribute.ien;

            if ist ~= 0 || ien ~= ni
                g{bid+1}.patch{pid+1}.attribute.ist = ni - ien;
                g{bid+1}.patch{pid+1}.attribute.ien = ni - ist;
            end
        end        
    end   
    
end
% size(bid_flip)
% bid_flip
end