function [coincident neighbour] = ts_get_connectivity(g)
% A function to loop over all patches and connectivity lists and write data
% into a large indexed cell array for speed of recalling
% This function only returns data that lies in the surfaces of the passage

% Initialise cell arrays
coincident = cell(length(g),1);
neighbour.diagonal = cell(length(g),1);
% neighbour.square = cell(length(g),1);

for b = 1:length(g)
    coincident{b} = cell(g{b}.attribute.ni, g{b}.attribute.nj,...
        g{b}.attribute.nk);
    neighbour.diagonal{b} = cell(g{b}.attribute.ni, g{b}.attribute.nj,...
        g{b}.attribute.nk);
%     neighbour.square{b} = cell(g{b}.attribute.ni, g{b}.attribute.nj,...
%         g{b}.attribute.nk);
end

% Loop over all clists and patches and write connectivity data into indexed array
% This represents only coincident data, although it may be spaced by a
% pitch
disp('Defining Coincidents')
for bid = ts_get_bids(g)
    b = ts_get_block(bid,g);
    
    % clists
    for clid = ts_get_clids(b)
        cl = ts_get_clist(bid,clid,g);
        
        for n = 1:cl.attribute.n
            I = cl.ilist(n); J = cl.jlist(n); K = cl.klist(n);
            
            if g{bid+1}.mwall(I+1,J+1,K+1) == 0
                nxI = cl.nxilist(n); nxJ = cl.nxjlist(n); nxK = cl.nxklist(n);
                nxbid = cl.attribute.nxbid;

                % Write connectivity arrays
                coincident{bid+1}{I+1,J+1,K+1} = ...
                    [coincident{bid+1}{I+1,J+1,K+1};...
                    nxbid nxI nxJ nxK];
                coincident{bid+1}{I+1,J+1,K+1} = unique(...
                    coincident{bid+1}{I+1,J+1,K+1} , 'rows');
            end

        end
    end

    % Periodic patches
    for pid = ts_get_pids(b)
        p0 = ts_get_patch(bid,pid,g);
        att = p0.attribute;
        if att.kind == 5
            % Get attributes
            ien0 = att.ien;  ist0 = att.ist;
            jen0 = att.jen;  jst0 = att.jst;
            ken0 = att.ken;  kst0 = att.kst;
            bid0 = att.bid; nxbid0 = att.nxbid; nxpid0 = att.nxpid;
            idir0 = att.idir; jdir0 = att.jdir; kdir0 = att.kdir;
            
            p1 = ts_get_patch(nxbid0, nxpid0, g);
            att = p1.attribute;
            ien1 = att.ien;  ist1 = att.ist;
            jen1 = att.jen;  jst1 = att.jst;
            ken1 = att.ken;  kst1 = att.kst;
   
            
            % Define array of absolute indices accross the patch
            index0 = zeros(ien0-ist0-1, jen0-jst0-1, ken0-kst0-1, 3, 'int32');
            index1 = zeros(ien1-ist1-1, jen1-jst1-1, ken1-kst1-1, 3, 'int32');

            for i = ist0:ien0-1
                for j = jst0:jen0-1
                    for k = kst0:ken0-1
                        index0(i-ist0+1,j-jst0+1,k-kst0+1,:) = [i j k];
                    end
                end
            end
            
            for i = ist1:ien1-1
                for j = jst1:jen1-1
                    for k = kst1:ken1-1
                        index1(i-ist1+1,j-jst1+1,k-kst1+1,:) = [i j k];
                    end
                end
            end


            % Deal with idir
            if idir0 == 1
                index1 = permute(index1,[2 1 3 4]);
            elseif idir0 == 2
                index1 = permute(index1,[3 2 1 4]);
            elseif idir0 == 3
                index1 = flipdim(index1,1);
            elseif idir0 == 4
                index1 = flipdim(index1,2);
                index1 = permute(index1,[2 1 3 4]);
            elseif idir0 == 5
                index1 = flipdim(index1,3);
                index1 = permute(index1,[3 2 1 4]);
            end
            
            % Deal with jdir
            if jdir0 == 0
                index1 = permute(index1,[2 1 3 4]);
            elseif jdir0 == 2
                index1 = permute(index1,[1 3 2 4]);
            elseif jdir0 == 3
                index1 = flipdim(index1,1);
                index1 = permute(index1,[2 1 3 4]);
            elseif jdir0 == 4
                index1 = flipdim(index1,2);
            elseif jdir0 == 5
                index1 = flipdim(index1,3);
                index1 = permute(index1,[1 3 2 4]);
            end
            
            % Deal with kdir
            if kdir0 == 0
                index1 = permute(index1,[3 2 1 4]);
            elseif kdir0 == 1
                index1 = permute(index1,[1 3 2 4]);
            elseif kdir0 == 3
                index1 = flipdim(index1,1);
                index1 = permute(index1,[3 2 1 4]);
            elseif kdir0 == 4
                index1 = flipdim(index1,2);
                index1 = permute(index1,[1 3 2 4]);
            elseif kdir0 == 5
                index1 = flipdim(index1,3);
            end            
            
            % Loop over index0 and look up corresponding points in index1
            for i = 1:size(index0,1)
                for j = 1:size(index0,2)
                    for k = 1:size(index0,3)
                        I = index0(i,j,k,1); J = index0(i,j,k,2); K = index0(i,j,k,3); 
                        if g{bid+1}.mwall(I+1,J+1,K+1) == 0
                            nxI = index1(i,j,k,1); nxJ = index1(i,j,k,2); nxK = index1(i,j,k,3); 
                            nxbid = nxbid0;

                            coincident{bid0+1}{I+1,J+1,K+1} = ...
                                [coincident{bid0+1}{I+1,J+1,K+1};...
                                int16([nxbid nxI nxJ nxK])];

                            coincident{bid0+1}{I+1,J+1,K+1} = unique(...
                                coincident{bid0+1}{I+1,J+1,K+1} , 'rows');
                        end
                    end
                end
            end
                                              
        end
    end
    
end

disp('Adding edges to coincidents')

% Add in coincident points that are left out due to a wall lying along a
% block edge
for bid = ts_get_bids(g)
    emp = cellfun(@isempty,coincident{bid+1});
    [I,J,K] = ind2sub(size(coincident{bid+1}), find(emp == 0));
    for m = 1:length(I)
        i = I(m)-1; j = J(m)-1; k = K(m)-1;
        for n = 1:size(coincident{bid+1}{i+1,j+1,k+1},1)
            v = coincident{bid+1}{i+1,j+1,k+1}(n,:);
            coincident{v(1)+1}{v(2)+1,v(3)+1,v(4)+1} =unique([coincident{v(1)+1}{v(2)+1,v(3)+1,v(4)+1}...
                ; [bid i j k]],'rows');
        end
    end
end

for bid = ts_get_bids(g)
    emp = cellfun(@isempty,coincident{bid+1});
    [I,J,K] = ind2sub(size(coincident{bid+1}), find(emp == 0));
    for m = 1:length(I)
        i = I(m)-1; j = J(m)-1; k = K(m)-1;
        for n = 1:size(coincident{bid+1}{i+1,j+1,k+1},1)
            v = coincident{bid+1}{i+1,j+1,k+1}(n,:);
            coincident{v(1)+1}{v(2)+1,v(3)+1,v(4)+1} =unique([coincident{v(1)+1}{v(2)+1,v(3)+1,v(4)+1}...
                ; coincident{bid+1}{i+1,j+1,k+1}],'rows');
        end
    end
end

% Change mwall to update edges of blocks that lie on an edge of the wall
for bid = ts_get_bids(g)
    emp = cellfun(@isempty,coincident{bid+1});
    [I,J,K] = ind2sub(size(coincident{bid+1}), find(emp == 0));
    for m = 1:length(I)
        i = I(m)-1; j = J(m)-1; k = K(m)-1;
        for n = 1:size(coincident{bid+1}{i+1,j+1,k+1},1)
            v = coincident{bid+1}{i+1,j+1,k+1}(n,:);
            if g{v(1)+1}.mwall(v(2)+1,v(3)+1,v(4)+1) == 0
                g{bid+1}.mwall(i+1,j+1,k+1) = 0;
            end
        end
    end
end
            
mwall = cell(size(g));
for bid = ts_get_bids(g)
    % Create new mwall array to remove results outside the bounds
    mwall{bid+1} = ones(size(g{bid+1}.mwall)+2);
    mwall{bid+1}(2:end-1,2:end-1,2:end-1) = g{bid+1}.mwall;
end

% Create neighbours array for every point including coincident data
disp('Initialising Neighbours')
for bid = ts_get_bids(g)
    i = 0;
    for j = 0:g{bid+1}.attribute.nj-1
        for k = 0:g{bid+1}.attribute.nk-1
            if g{bid+1}.mwall(i+1,j+1,k+1) == 0
                neighbour.diagonal{bid+1}{i+1,j+1,k+1} = [neighbour.diagonal{bid+1}{i+1,j+1,k+1} ; ...
                    int16([bid i j+1 k; bid i j+1 k+1; bid i j k+1; bid i j-1 k+1; ... 
                    bid i j-1 k; bid i j-1 k-1; bid i j k-1; bid i j+1 k-1])];
%                 neighbour.square{bid+1}{i+1,j+1,k+1} = [neighbour.square{bid+1}{i+1,j+1,k+1} ; ...
%                     int16([bid i j+1 k; bid i j k+1; bid i j-1 k; bid i j k-1])];
            end
        end
    end
    i = g{bid+1}.attribute.ni-1;
    for j = 0:g{bid+1}.attribute.nj-1
        for k = 0:g{bid+1}.attribute.nk-1
            if g{bid+1}.mwall(i+1,j+1,k+1) == 0
                neighbour.diagonal{bid+1}{i+1,j+1,k+1} = [neighbour.diagonal{bid+1}{i+1,j+1,k+1} ; ...
                    int16([bid i j+1 k; bid i j+1 k+1; bid i j k+1; bid i j-1 k+1; ... 
                    bid i j-1 k; bid i j-1 k-1; bid i j k-1; bid i j+1 k-1])];
%                 neighbour.square{bid+1}{i+1,j+1,k+1} = [neighbour.square{bid+1}{i+1,j+1,k+1} ; ...
%                     int16([bid i j+1 k; bid i j k+1; bid i j-1 k; bid i j k-1])];
            end
        end
    end    
    j = 0;
    for i = 0:g{bid+1}.attribute.ni-1
        for k = 0:g{bid+1}.attribute.nk-1
            if g{bid+1}.mwall(i+1,j+1,k+1) == 0
                neighbour.diagonal{bid+1}{i+1,j+1,k+1} = [neighbour.diagonal{bid+1}{i+1,j+1,k+1} ; ...
                    int16([bid i+1 j k; bid i+1 j k+1; bid i j k+1; bid i-1 j k+1; ... 
                    bid i-1 j k; bid i-1 j k-1; bid i j k-1; bid i+1 j k-1])];
%                 neighbour.square{bid+1}{i+1,j+1,k+1} = [neighbour.square{bid+1}{i+1,j+1,k+1} ; ...
%                     int16([bid i+1 j k; bid i j k+1; bid i-1 j k; bid i j k-1])];
            end
        end
    end
    j = g{bid+1}.attribute.nj-1;
    for i = 0:g{bid+1}.attribute.ni-1
        for k = 0:g{bid+1}.attribute.nk-1
            if g{bid+1}.mwall(i+1,j+1,k+1) == 0
                neighbour.diagonal{bid+1}{i+1,j+1,k+1} = [neighbour.diagonal{bid+1}{i+1,j+1,k+1} ; ...
                    int16([bid i+1 j k; bid i+1 j k+1; bid i j k+1; bid i-1 j k+1; ... 
                    bid i-1 j k; bid i-1 j k-1; bid i j k-1; bid i+1 j k-1])];
%                 neighbour.square{bid+1}{i+1,j+1,k+1} = [neighbour.square{bid+1}{i+1,j+1,k+1} ; ...
%                     int16([bid i+1 j k; bid i j k+1; bid i-1 j k; bid i j k-1])];
            end
        end
    end    
    k = 0;
    for i = 0:g{bid+1}.attribute.ni-1
        for j = 0:g{bid+1}.attribute.nj-1
            if g{bid+1}.mwall(i+1,j+1,k+1) == 0
                neighbour.diagonal{bid+1}{i+1,j+1,k+1} = [neighbour.diagonal{bid+1}{i+1,j+1,k+1} ; ...
                    int16([bid i+1 j k; bid i+1 j+1 k; bid i j+1 k; bid i-1 j+1 k; ... 
                    bid i-1 j k; bid i-1 j-1 k; bid i j-1 k; bid i+1 j-1 k])];
%                 neighbour.square{bid+1}{i+1,j+1,k+1} = [neighbour.square{bid+1}{i+1,j+1,k+1} ; ...
%                     int16([bid i+1 j k; bid i j+1 k; bid i-1 j k; bid i j-1 k])];
            end
        end
    end
    k = g{bid+1}.attribute.nk-1;
    for i = 0:g{bid+1}.attribute.ni-1
        for j = 0:g{bid+1}.attribute.nj-1
            if g{bid+1}.mwall(i+1,j+1,k+1) == 0
                neighbour.diagonal{bid+1}{i+1,j+1,k+1} = [neighbour.diagonal{bid+1}{i+1,j+1,k+1} ; ...
                    int16([bid i+1 j k; bid i+1 j+1 k; bid i j+1 k; bid i-1 j+1 k; ... 
                    bid i-1 j k; bid i-1 j-1 k; bid i j-1 k; bid i+1 j-1 k])];
%                 neighbour.square{bid+1}{i+1,j+1,k+1} = [neighbour.square{bid+1}{i+1,j+1,k+1} ; ...
%                     int16([bid i+1 j k; bid i j+1 k; bid i-1 j k; bid i j-1 k])];
            end
        end
    end

end

disp('Trimming neighbours')

for bid = ts_get_bids(g)
    [I,J,K] = ind2sub(size(g{bid+1}.mwall), find(g{bid+1}.mwall == 0));
    for m = 1:length(I)
        i = I(m)-1; j = J(m)-1; k = K(m)-1;
                    
        % Put in points surrounding coincident points
        if isempty(coincident{bid+1}{i+1,j+1,k+1}) == 0
            for n = 1:size(coincident{bid+1}{i+1,j+1,k+1},1)
                v = coincident{bid+1}{i+1,j+1,k+1}(n,:);
                neighbour.diagonal{bid+1}{i+1,j+1,k+1} = [neighbour.diagonal{bid+1}{i+1,j+1,k+1} ;...
                    neighbour.diagonal{v(1)+1}{v(2)+1,v(3)+1,v(4)+1}];
%                 neighbour.square{bid+1}{i+1,j+1,k+1} = [neighbour.square{bid+1}{i+1,j+1,k+1} ;...
%                     neighbour.square{v(1)+1}{v(2)+1,v(3)+1,v(4)+1}];                
            end
        end
        
        % Remove points outside of bounds
        for n = size(neighbour.diagonal{bid+1}{i+1,j+1,k+1},1):-1:1
            v = neighbour.diagonal{bid+1}{i+1,j+1,k+1}(n,:);
            if mwall{v(1)+1}(v(2)+2,v(3)+2,v(4)+2) ~= 0
                neighbour.diagonal{bid+1}{i+1,j+1,k+1}(n,:) = [];
            end
        end
        
%         for n = size(neighbour.square{bid+1}{i+1,j+1,k+1},1):-1:1
%             v = neighbour.square{bid+1}{i+1,j+1,k+1}(n,:);
%             if mwall{v(1)+1}(v(2)+2,v(3)+2,v(4)+2) ~= 0
%                 neighbour.square{bid+1}{i+1,j+1,k+1}(n,:) = [];
%             end
%         end        

        % Put in points that are coincident with the neighbours
        for n = 1:size(neighbour.diagonal{bid+1}{i+1,j+1,k+1},1)
            v = neighbour.diagonal{bid+1}{i+1,j+1,k+1}(n,:);
            neighbour.diagonal{bid+1}{i+1,j+1,k+1} = [neighbour.diagonal{bid+1}{i+1,j+1,k+1} ;...
                coincident{v(1)+1}{v(2)+1,v(3)+1,v(4)+1}];
        end
        
%         for n = 1:size(neighbour.square{bid+1}{i+1,j+1,k+1},1)
%             v = neighbour.square{bid+1}{i+1,j+1,k+1}(n,:);
%             neighbour.square{bid+1}{i+1,j+1,k+1} = [neighbour.square{bid+1}{i+1,j+1,k+1} ;...
%                 coincident{v(1)+1}{v(2)+1,v(3)+1,v(4)+1}];
%         end        
        
        % Take only unique points
        neighbour.diagonal{bid+1}{i+1,j+1,k+1} = unique(neighbour.diagonal{bid+1}{i+1,j+1,k+1},'rows');
%         neighbour.square{bid+1}{i+1,j+1,k+1} = unique(neighbour.square{bid+1}{i+1,j+1,k+1},'rows');
    end
end

neighbour = neighbour.diagonal;

end

