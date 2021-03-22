function C = ts_cut_rows(g,xr_hub,xr_cas,cut_ss,cut_tip,x)
% TS_CUT_ROWS  Take structured cuts on individual rows and assemble into a cell array

% Default to extract suction side block
if exist('cut_ss','var') == 0
    cut_ss = 1;
end

% Default to not extract tip cuts
if exist('cut_tip','var') == 0
    cut_tip = 0;
end

% Default to not take any unstructured axial cuts
if exist('x','var') == 0
    x = [];
end

% Get row and o-mesh bids
[B,bid_os] = ts_get_row_bids(g,xr_hub,xr_cas);

% Prellocate cell array
nn = length(B); C = cell(nn,1);

% Take o-mesh cuts
for n = 1:nn
    C{n}.blade = ts_structured_cut(g,bid_os(n),1,'en',1,'en',1,1);
end

% Find and take cuts one cell away from the mixing planes
for n = 1:nn
    for bid = B{n}
        for pid = ts_get_pids(g,bid)
            
            % Row inlet
            if any(g{bid+1}.patch{pid+1}.attribute.kind == [2 3]) && ...
                    g{bid+1}.patch{pid+1}.attribute.ist == 0
                C{n}.inlet = ts_structured_cut(g,bid,2,2,1,'en',1,'en');
            end            
            
            % Row outlet
            if any(g{bid+1}.patch{pid+1}.attribute.kind == [2 3]) && ...
                    g{bid+1}.patch{pid+1}.attribute.ien == g{bid+1}.attribute.ni
                i = g{bid+1}.attribute.ni - 1;
                C{n}.outlet = ts_structured_cut(g,bid,i,i,1,'en',1,'en');
            end
        end
    end
end

% Find inlet and exit planes to domain
for bid = ts_get_bids(g)
    for pid = ts_get_pids(g,bid)
        if g{bid+1}.patch{pid+1}.attribute.kind == 0
            C{1}.inlet = ts_structured_cut(g,bid,2,2,1,'en',1,'en');
        elseif g{bid+1}.patch{pid+1}.attribute.kind == 1
            i = g{bid+1}.attribute.ni - 1;
            C{end}.outlet = ts_structured_cut(g,bid,i,i,1,'en',1,'en');
        end
    end
end

% Find suction surface h-blocks
if cut_ss == 1
    for n = 1:nn

        % Measure pressure on both sides of the blade
        c = ts_secondary(C{n}.blade); j_mid = round(double(c.attribute.nj)/2); 
        [~,i_le] = min(c.x(:,j_mid));
        p_1 = mean(mean(c.P(1:i_le,:))); p_2 = mean(mean(c.P(i_le:end,:)));

        % Determine an index in the middle of the suction surface
        if p_1 < p_2; i_ss = round(0.5 * i_le); else i_ss = round(1.5 * i_le); end;

        % Find the patch in the middle of the suction side
        bid_o = c.attribute.bid;
        for pid_o = ts_get_pids(c)
            i_o = [c.patch{pid_o+1}.attribute.ist+1 c.patch{pid_o+1}.attribute.ien];
            if i_o(1) < i_ss && i_o(2) > i_ss && c.patch{pid_o+1}.attribute.kst == c.attribute.nk-1
                break
            end
        end

        % Extract parameters describing the join between o and h-block
        bid_ss = c.patch{pid_o+1}.attribute.nxbid; pid_ss = c.patch{pid_o+1}.attribute.nxpid;
        i_ss = [g{bid_ss+1}.patch{pid_ss+1}.attribute.ist+1 g{bid_ss+1}.patch{pid_ss+1}.attribute.ien];
        j_ss = [g{bid_ss+1}.patch{pid_ss+1}.attribute.jst+1 g{bid_ss+1}.patch{pid_ss+1}.attribute.jen];

        % Join blocks together
        varnames = fieldnames(g{bid_o+1});
        b = g{bid_o+1};
        for v = 1:length(varnames)
            if isnumeric(b.(varnames{v})) == 1
                b.(varnames{v}) = double(cat(3,g{bid_o+1}.(varnames{v})(i_o(1):i_o(2),:,:),...
                    g{bid_ss+1}.(varnames{v})(i_ss(1):i_ss(2),j_ss(1):j_ss(2),2:end)));
            end
        end

        % Record new block and its size
        [b.attribute.ni,b.attribute.nj,b.attribute.nk] = size(b.x); 
        C{n}.ss = b;

    end
end

% Extract cuts around the tip gap
if cut_tip == 1
    for n = 1:nn
        
        % Find a block joined to the blade o-mesh in the right orientation
        for bid = B{n}
            for pid = ts_get_pids(g,bid)
                if g{bid+1}.patch{pid+1}.attribute.nxbid == bid_os(n) && ...
                        g{bid+1}.patch{pid+1}.attribute.idir == 0 && ...
                        g{bid+1}.patch{pid+1}.attribute.jdir == 1 && ...
                        g{bid+1}.patch{pid+1}.attribute.ken == g{bid+1}.attribute.nk
                    
                    % Take the cut at the join with the o-mesh
                    C{n}.tip = ts_structured_cut(g,bid,1,'en',1,'en','en','en');
                    
                end
            end
        end
        
    end
end

% Extract unstructured axial cuts
if isempty(x) == 0
    for n = 1:nn
        
        % Find axial positions in the current row
        q = x > min(C{n}.inlet.x(:)) & x < max(C{n}.outlet.x(:));
        
        % Check if any cuts lie in the current row
        if sum(q) > 0
            
            % Loop over axial positions and take cuts
            C{n}.axial = cell(sum(q),1);
            is = find(q);
            for i = 1:length(is)
                C{n}.axial{i} = ts_unstructured_cut(g,x(is(i)));
            end
        end
    end
end

% Record meridional lines
for n = 1:nn
    C{n}.xr_hub = xr_hub; C{n}.xr_cas = xr_cas;
end


end