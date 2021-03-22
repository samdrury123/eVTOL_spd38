function cut = ts_structunstruct_cut(g,xr,bids,IJ,recomb, splits)
% Take an unstructured cut along a plane from (x0, r0) to (x1,r1)

% Specify block list to speed up code
if exist('bids','var') == 0
    bids = ts_get_bids(g);
elseif isempty(bids)
    bids = ts_get_bids(g);
end

% IJswitch changes between i constant lines (1) blade2blade or j constant
% lines (2) traverse plane
if exist('IJ','var') == 0
    IJ = 2;
end

% Change output format, ifrecomb == 1 then cell arrays will be recombined
% in the passage, if not then output will be in  cell array format
if exist('recomb','var') == 0
    recomb = 0;
end

tol = 1e-5;

% Get actual variable names of interest
varnames = fieldnames(g{1});
varnames_new = {};
interpnames = {};
for v = 1:length(varnames)
    varname = varnames{v};
    if ndims(g{1}.(varname)) == 3;
        varnames_new = [varnames_new ; varname];
        if strcmp(varname,'x') ~= 1 && strcmp(varname,'r') ~= 1 && strcmp(varname,'rt') ~= 1
            interpnames = [interpnames ; varname];
        end
    end
end
interpnames = [interpnames ; 'rt'];
varnames = varnames_new;

% Create any additional splits as specified
if exist('splits','var')
    for n = 1:size(splits,1)
        bid = splits(n,1);
        i = splits(n,2);
        for v = 1:length(varnames)
            varname = varnames{v};
            b.(varname) = g{bid+1}.(varname)(1:i+1,:,:);
            g{bid+1}.(varname) = g{bid+1}.(varname)(i+1:end,:,:);
        end
        
        g{bid+1}.attribute.ni = size(g{bid+1}.x,1);
        g{bid+1}.attribute.nj = size(g{bid+1}.x,2);
        g{bid+1}.attribute.nk = size(g{bid+1}.x,3);
        
        b.attribute.ni = size(b.x,1);
        b.attribute.nj = size(b.x,2);
        b.attribute.nk = size(b.x,3);
        b.attribute.bid = size(g,1);
        b.bv = g{bid+1}.bv;
        
        g{end+1} = b;
        
        bids = [bids size(g,1)-1];
    end
end        

% Split o-meshes into 4 and re-order LE and TE pieces
% Find o-mesh bids and corner points
n = 1;
bid_os = [];
for bid = bids
    for pid = ts_get_pids(bid,g)
        if g{bid+1}.patch{pid+1}.attribute.kind == 5
            if g{bid+1}.patch{pid+1}.attribute.nxbid == g{bid+1}.patch{pid+1}.attribute.bid
                if g{bid+1}.patch{pid+1}.attribute.ist == g{bid+1}.attribute.ni - 1 && g{bid+1}.patch{pid+1}.attribute.ien == g{bid+1}.attribute.ni
                    bid_os(n) = bid;
                    Qtemp = [];
                    % Get corner points
%                     for pid2 = ts_get_pids(bid,g)
%                         if g{bid+1}.patch{pid2+1}.attribute.nxbid ~= g{bid+1}.patch{pid2+1}.attribute.bid
%                             if g{bid+1}.patch{pid2+1}.attribute.ien ~= g{bid+1}.attribute.ni
%                                 Qtemp = [Qtemp g{bid+1}.patch{pid2+1}.attribute.ien + 1];
%                             end
%                         end
%                     end
%                     Q(n,:) = sort(Qtemp);
                    % Get corner points based on which points are closest
                    % to corners of a large square
                    nj = round(g{bid+1}.attribute.nj/2);
%                     q = [reshape(g{bid+1}.x(:,nj,end) - g{bid+1}.x(:,nj,1),[],1,1),...
%                         reshape(g{bid+1}.r(:,nj,end) - g{bid+1}.r(:,nj,1),[],1,1),...
%                         reshape(g{bid+1}.rt(:,nj,end) - g{bid+1}.rt(:,nj,1),[],1,1)];
%                     norm_q = (q(:,1).^2 + q(:,2).^2 + q(:,3).^2).^0.5;
%                     q = q./repmat(norm_q,1,3);
%                     [~,I1] = max(dot(q,repmat([1 0 1],size(q,1),1),2));
%                     [~,I2] = max(dot(q,repmat([-1 0 1],size(q,1),1),2));
%                     [~,I3] = max(dot(q,repmat([-1 0 -1],size(q,1),1),2));
%                     [~,I4] = max(dot(q,repmat([1 0 -1],size(q,1),1),2));
                    q = [reshape(g{bid+1}.x(:,nj,end),[],1,1),...
                        reshape(g{bid+1}.r(:,nj,end),[],1,1),...
                        reshape(g{bid+1}.rt(:,nj,end),[],1,1)];
                    q1 = q - 100*[ones(size(q,1),1) zeros(size(q,1),1) ones(size(q,1),1)];
                    [~,I1] = min(q1(:,1).^2 + q1(:,2).^2 + q1(:,3).^2);
                    q2 = q - 100*[-ones(size(q,1),1) zeros(size(q,1),1) ones(size(q,1),1)];
                    [~,I2] = min(q2(:,1).^2 + q2(:,2).^2 + q2(:,3).^2);
%                     [~, I2] = min(q(:,1)); % Min x coordinate
                    q3 = q - 100*[-ones(size(q,1),1) zeros(size(q,1),1) -ones(size(q,1),1)];
                    [~,I3] = min(q3(:,1).^2 + q3(:,2).^2 + q3(:,3).^2);
%                     [~, I3] = min(q(:,2)); % Min r coordinate
                    q4 = q - 100*[ones(size(q,1),1) zeros(size(q,1),1) -ones(size(q,1),1)];
                    [~,I4] = min(q4(:,1).^2 + q4(:,2).^2 + q4(:,3).^2);
                    Q(n,:) = sort([I1 I2 I3 I4]);

                    n = n+1;
                end
            end
        end
    end
end

% Split up into 4 and cat onto end of g
n = 1;
for bid = bid_os
    ni = g{bid+1}.attribute.ni;
    
    b1 = g{bid+1};
    b2 = g{bid+1};
    b3 = g{bid+1};
    b4 = g{bid+1};
    
    % Use corner points to get variables
    for v = 1:length(varnames)
        varname = varnames{v};
        b1.(varname) = b1.(varname)(Q(n,1):Q(n,2),:,:);
        b2.(varname) = b2.(varname)(Q(n,2):Q(n,3),:,:);
        b3.(varname) = b3.(varname)(Q(n,3):Q(n,4),:,:);
        b4.(varname) = b4.(varname)([Q(n,4):ni 2:Q(n,1)],:,:);
    end
    
    % Adjust ni indices
    b1.attribute.ni = size(b1.(varname),1);
    b2.attribute.ni = size(b2.(varname),1);
    b3.attribute.ni = size(b3.(varname),1);
    b4.attribute.ni = size(b4.(varname),1);
    
    % Flip i and k directions on b2 and b4
    for v = 1:length(varnames)
        varname = varnames{v};
        b2.(varname) = permute(b2.(varname),[3 2 1]);
        b4.(varname) = permute(b4.(varname),[3 2 1]);
    end
    
    b2.attribute.ni = size(b2.(varname),1);
    b2.attribute.nk = size(b2.(varname),3);
    b4.attribute.ni = size(b4.(varname),1);
    b4.attribute.nk = size(b4.(varname),3);
    
    % Reassign block numbers
    b2.attribute.bid = size(g,1);
    b3.attribute.bid = size(g,1)+1;
    b4.attribute.bid = size(g,1)+2;
    
    % Reassign into g cell array
    g{bid+1} = b1;
    g = [g ; b2 ; b3 ; b4];
    bids = [bids b2.attribute.bid b3.attribute.bid b4.attribute.bid];
    
    n = n+1;
end

for bid = bids
    % Determine which indices are in the vicinity of the cut permute and
    % flip block so x is increasing in the i direction, r in the j and rt
    % in the k
    ni = g{bid+1}.attribute.ni; nj = g{bid+1}.attribute.nj; nk = g{bid+1}.attribute.nk;
    x = repmat(mean(xr(:,1)),[ni nj nk]); r = repmat(mean(xr(:,2)),[ni nj nk]); rt = repmat(mean(mean(mean(g{bid+1}.rt))),[ni nj nk]);
    dist = (g{bid+1}.x-x).^2 + (g{bid+1}.r-r).^2 + (g{bid+1}.rt-rt).^2; 
    [~, n] = min(reshape(dist,[],1,1));
    [I,J,K] = ind2sub([ni nj nk],n); I = I-1; J = J-1; K = K-1;
    
    I = [I-1 I I+1];
    if min(I) == -1
        I = I+1;
    elseif max(I) == ni;
        I = I-1;
    end
    
    J = [J-1 J J+1];
    if min(J) == -1
        J = J+1;
    elseif max(J) == nj;
        J = J-1;
    end
    
    K = [K-1 K K+1];
    if min(K) == -1
        K = K+1;
    elseif max(K) == nk;
        K = K-1;
    end    

    % Get vector in grid directions
    i = [g{bid+1}.x(I(3)+1,J(2)+1,K(2)+1) - g{bid+1}.x(I(1)+1,J(2)+1,K(2)+1);
        g{bid+1}.r(I(3)+1,J(2)+1,K(2)+1) - g{bid+1}.r(I(1)+1,J(2)+1,K(2)+1);
        g{bid+1}.rt(I(3)+1,J(2)+1,K(2)+1) - g{bid+1}.rt(I(1)+1,J(2)+1,K(2)+1)];
    i = i / norm(i);
    j = [g{bid+1}.x(I(2)+1,J(3)+1,K(2)+1) - g{bid+1}.x(I(2)+1,J(1)+1,K(2)+1);
        g{bid+1}.r(I(2)+1,J(3)+1,K(2)+1) - g{bid+1}.r(I(2)+1,J(1)+1,K(2)+1);
        g{bid+1}.rt(I(2)+1,J(3)+1,K(2)+1) - g{bid+1}.rt(I(2)+1,J(1)+1,K(2)+1)];
    j = j / norm(j);
    k = [g{bid+1}.x(I(2)+1,J(2)+1,K(3)+1) - g{bid+1}.x(I(2)+1+1,J(2)+1,K(1)+1);
        g{bid+1}.r(I(2)+1,J(2)+1,K(3)+1) - g{bid+1}.r(I(2)+1+1,J(2)+1,K(1)+1);
        g{bid+1}.rt(I(2)+1,J(2)+1,K(3)+1) - g{bid+1}.rt(I(2)+1+1,J(2)+1,K(1)+1)];
    k = k / norm(k);
    
    % Find vector most aligned with x direction and set i grid direction
    % along that dimension
%     xi = abs(dot(i,[1;0;0])); xj = abs(dot(j,[1;0;0])); xk = abs(dot(k,[1;0;0]));
%     if max([xi xj xk]) == xj
%         for v = 1:length(varnames)
%             varname = varnames{v};
%             g{bid+1}.(varname) = permute(g{bid+1}.(varname),[2 1 3]);
%         end
%         temp = i; i = j; j = temp;
%     elseif max([xi xj xk]) == xk
%         for v = 1:length(varnames)
%             varname = varnames{v};
%             g{bid+1}.(varname) = permute(g{bid+1}.(varname),[3 2 1]);
%         end
%         temp = i; i = k; k = temp;
%     end
%     % Find remaining vector most aligned with the r direction and set j
%     rj = abs(dot(j,[0;1;0])); rk = abs(dot(k,[0;0;1]));
%     if max([rj rk]) == rk
%         for v = 1:length(varnames)
%             varname = varnames{v};
%             g{bid+1}.(varname) = permute(g{bid+1}.(varname),[1 3 2]);
%         end
%         temp = j; j = k; k = temp;
%     end
    
    % Flip blocks to increase in respective directions    
    
    % Flip all blocks in i direction if not in the increasing x direction
    if dot(i,[1;0;0]) < 0
        for v = 1:length(varnames)
            varname = varnames{v};
            g{bid+1}.(varname) = flipdim(g{bid+1}.(varname),1);
        end
    end
    
    % Flip all blocks in j direction if not in the increasing r direction
    if dot(j,[0;1;0]) < 0
        for v = 1:length(varnames)
            varname = varnames{v};
            g{bid+1}.(varname) = flipdim(g{bid+1}.(varname),2);
        end
    end
    
    % Flip all blocks in k direction if not in the increasing rt direction
    if dot(k,[0;0;1]) < 0
        for v = 1:length(varnames)
            varname = varnames{v};
            g{bid+1}.(varname) = flipdim(g{bid+1}.(varname),3);
        end
    end
    
    % Reset ni, nj, nk attributes
    g{bid+1}.attribute.ni = size(g{bid+1}.x,1);
    g{bid+1}.attribute.nj = size(g{bid+1}.x,2);
    g{bid+1}.attribute.nk = size(g{bid+1}.x,3);
    
    % Turn everything to doubles
    for v = 1:length(varnames)
        varname = varnames{v};
        g{bid+1}.(varname) = double(g{bid+1}.(varname));
    end
end

% Initialise cut cell array
cut = cell(size(g));

if IJ == 2
    % Loop over each block and check intersection of the line of the
    % unstructured cut with the constant lines in the block
    for bid = bids
        % Check if the line lies in the block or not
        [flag n] = inblock(g,bid,xr,2);
        if flag == 1
            
            % Use matrix operations to determine all intersecting
            % coordinates
            [cut{bid+1}.xi cut{bid+1}.ri cut{bid+1}.Ind] = findintersection2(g{bid+1}.x,g{bid+1}.r,xr(n,:),1);
            if isfield(g{1},'av')
                cut{bid+1}.av = g{1}.av;
            end
            cut{bid+1}.bv = g{bid+1}.bv;

            % Remove coordinates of arrays that do not intersect
            J_remove = sum(cut{bid+1}.xi,2)==0;
            cut{bid+1}.xi(J_remove,:) = nan;
            cut{bid+1}.ri(J_remove,:) = nan;
            cut{bid+1}.Ind(1:4,J_remove,:) = 1;            

            % Select 4 bounding elements
            for v = 1:length(varnames)
                varname = varnames{v};
                var = cat(1,2*g{bid+1}.(varname)(1,:,:)-g{bid+1}.(varname)(2,:,:)...
                    ,g{bid+1}.(varname)(:,:,:), ...
                    2*g{bid+1}.(varname)(end,:,:)-g{bid+1}.(varname)(end-1,:,:));                
                cut{bid+1}.(varname) = reshape(var(cut{bid+1}.Ind==1),4,[],g{bid+1}.attribute.nk);
                cut{bid+1}.(varname)(:,J_remove,:) = nan;
                cut{bid+1}.(varname) = permute(cut{bid+1}.(varname),[2 3 1]);
            end
            
            % Interpolate using a spline
            for v = 1:length(interpnames)
                varname = interpnames{v};            
                cut{bid+1}.(varname) = spline1(cut{bid+1}.x, cut{bid+1}.(varname), cut{bid+1}.xi);
            end

            cut{bid+1}.x = cut{bid+1}.xi;
            cut{bid+1}.r = cut{bid+1}.ri;

            % Remove data
            for v = 1:length(varnames)
                varname = varnames{v};
                cut{bid+1}.(varname) = cut{bid+1}.(varname)(:,:,1);
            end            
            
            cut{bid+1} = rmfield(cut{bid+1},{'ri', 'Ind', 'xi'});

        end
    end

    % Remove empty arrays from cut
    c_remove = [];
    for c = 1:length(cut)
        if isempty(cut{c}) == 1
            c_remove = [c_remove c];
        elseif size(cut{c}.x,1) == 0
            c_remove = [c_remove c];
        elseif sum(sum(isnan(cut{c}.x))) == numel(cut{c}.x)
            c_remove = [c_remove c];
        end
    end
    cut(c_remove) = [];

    if recomb == 1
        % Reconstruct a single mesh in the passage from array pieces
        % Determine whether recombining accross a blade or mixing plane.
        % Recombining in r or rt direction.
        rmax = zeros(size(cut));
        for c = 1:size(cut,1)
            rmax(c) = max(max(cut{c}.r));
        end

        if length(cut) == 1
            recomb_type = 2; % single array
        elseif isalmost(repmat(rmax,1,size(rmax,1)),repmat(rmax.',size(rmax,1),1),tol)
            recomb_type = 0; % rt
        else
            recomb_type = 1; % r
        end

        if recomb_type == 0
            % Shouldn't be any NaNs in the array and so recombination is
            % simple to do in the positive rt direction
            rtmin = [];
            for c = 1:size(cut,1)
                rtmin = [rtmin min(cut{c}.rt(1,:))];
            end
            
            % Sort and initialise data structure
            [~, c_order] = sort(rtmin);
            varnames = fieldnames(cut{1});
            cut_2 = cell2struct(cell(size(varnames)),varnames,1);

            for c = c_order
                if cut{c}.rt(1) > cut{c}.rt(end)
                    flip_t = 1;
                else
                    flip_t = 0;
                end
                
                for v = 1:length(varnames)
                    varname = varnames{v};
                    if flip_t == 1
                        cut_2.(varname) = [cut_2.(varname) flipdim(cut{c}.(varname),2)];
                    else
                        cut_2.(varname) = [cut_2.(varname) cut{c}.(varname)];
                    end
                end
            end
            cut = cut_2; cut.av = cut.av(1); cut.bv = cut.bv(1); 
            
            % Split and recombine in opposite direction if intersecting a blade
            rt_nondim = (cut.rt - repmat(min(cut.rt,[],2),1,size(cut.rt,2))) ./...
               (repmat(max(cut.rt,[],2),1,size(cut.rt,2)) - repmat(min(cut.rt,[],2),1,size(cut.rt,2)));
            blade_present = 0;
            rt_grad = zeros(size(rt_nondim(1,:)));
            for k = 2:size(rt_nondim,2)
                rt_grad(k) = rt_nondim(1,k) - rt_nondim(1,k-1);
            end

            % Detect whether a blade intersection is present
            [max_grad ksplit] = max(rt_grad);
            if max_grad > 5*mean(rt_grad)
                blade_present = 1;
            end
            
            % Recalculate theta coords and reorder flow properties
            if blade_present == 1
                nk = size(cut.r,2);
                nj = size(cut.r,1);     
                cut.t = cut.rt ./ cut.r; 
                t_new = zeros(nj,nk);
                t_new(:,1:nk-ksplit+1) = cut.t(:,ksplit:end);
                t_new(:,nk-ksplit+2:end) = cut.t(:,1:ksplit-1) + ...
                    2 * pi / double(cut.bv.nblade);
%                 t_new = t_new - mean(mean(t_new));
                cut.t = t_new; cut.rt = cut.t .* cut.r;

                cut.x(:,1:nk-ksplit+1) = cut.x(:,ksplit:end);
                cut.x(:,nk-ksplit+2:end) = cut.x(:,1:ksplit-1);
                
                for v = 1:length(interpnames)
                    interpname = interpnames{v};
                    if strcmp(interpname,'rt') == 0
                        varnew = zeros(nj,nk);
                        varnew(:,1:nk-ksplit+1) = cut.(interpname)(:,ksplit:end);
                        varnew(:,nk-ksplit+2:end) = cut.(interpname)(:,1:ksplit-1);
                        cut.(interpname) = varnew;
                    end
                end
            end
            
            % Remove lines of coincident nodes
            drt = diff(cut.rt(1,:));
            k_coinc = find(drt < drt(1) / 1000);
            for v = 1:length(varnames)
                varname = varnames{v};
                if strcmp(varname,'rt') == 0 && size(cut.(varname),3) > 1
                    var_mean = zeros(nj,length(k_coinc));
                    for p = 1:length(k_coinc)
                        k = k_coinc(p);
                        var_mean(:,p) = mean(cut.(varname)(:,k:k+1),2);
                    end
                    cut.(varname)(:,k_coinc) = var_mean;
                    cut.(varname)(:,k_coinc+1) = [];
                end
            end
            
        end

        if recomb_type == 1;
            % Split arrays into chunks in the case of a curved mixing plane
            % and multiple intersections
            rmin = []; nk = []; n = 1;
            
            for c = 1:size(cut,1)
%                 varnames = fieldnames(cut{c});
                J = sum(isnan([nan(1,size(cut{c}.r,2)) ; cut{c}.r ; nan(1,size(cut{c}.r,2))]),2) == 0;
                J = J - circshift(J,1);
                J_starts = find(J==1)-1;
                J_ends = find(J==-1)-2;
                
                for j = 1:length(J_starts)
                    for v = 1:length(varnames)
                        varname = varnames{v};
                        cut_2{n,1}.(varname) = cut{c}.(varname)(J_starts(j):J_ends(j),:);
                    end

                    bv = cut{c}.bv;
                    av = cut{1}.av;

                    rmin = [rmin min(min(cut_2{n,1}.r))];
                    nk = [nk size(cut_2{n,1}.r,2)];
                    n = n+1;
                end
            end
            
            % Remove unncessarary fields
%             unnecessary = {'dwallsq', 'misc', 'mwall', 'phi'};
%             for c = 1:size(cut_2,1)
%                 cut_2{c} = rmfield(cut_2{c},unnecessary);
%             end
            
            % If it does lie accross a mixing plane use an area average
            % on everything and return a single vector of a profile
            if ~isequal(repmat(nk.',1,size(nk,2)),repmat(nk,size(nk,2),1))
                for c = 1:size(cut_2,1)
                    for v = 1:length(interpnames)-1
                        varname = interpnames{v};
                        cut_2{c}.(varname) = ts_area_average(cut_2{c},varname,2).';
                    end
                    cut_2{c}.x = mean(cut_2{c}.x,2);
                    cut_2{c}.r = mean(cut_2{c}.r,2);
                    cut_2{c}.rt = max(cut_2{c}.rt,[],2) - min(cut_2{c}.rt,[],2);
                end
            end

            % Sort in terms of r coordinate and recombine
            [~, c_order] = sort(rmin);
            % Initialise new data structure
            cut = cell2struct(cell(size(varnames)),varnames,1);
            cut.av = av;
            cut.bv = bv;
            
            for c = c_order
                for v = 1:length(varnames)
                    varname = varnames{v};
                    cut.(varname) = [cut.(varname) ; cut_2{c}.(varname)];
                end
            end
            
        end
        
        if recomb_type == 2
            % Remove any stripes of nans in the case of a single array
            % which is bounded by the r coordinates
            J_remove = sum(isnan(cut{1}.x),2) ~= 0;
            for v = 1:length(varnames)
                varname = varnames{v};
                cut{1}.(varname)(J_remove,:) = [];
            end
        end                    
            
        
%         for c = 1:size(cut,1)
%             bw = cut{c}.x;
%             bw(isnan(bw)) = 0;
%             wb = cut{c}.x;
%             wb(~isnan(wb)) = 0;
%             wb(isnan(wb)) = 1;
%             cut{c}.Bf = bwboundaries(bw);
%             [cut{c}.Bn, ~, ~, cut{c}.An] = bwboundaries(wb);
%         end
        
        
%         nan_info = zeros(size(cut,1),3);
%         fl_info = zeros(size(cut,1),3);
%         for c = 1:size(cut,1)
%             if isempty(cut{c}) == 0
%                 nan_info(c,1) = sum(sum(isnan(cut{c}.x)));
%                 nan_info(c,2) = mean(cut{c}.x(isnan(cut{c}.x)));
%                 nan_info(c,3) = mean(cut{c}.r(isnan(cut{c}.r)));
% 
%                 fl_info(c,1) = sum(sum(~isnan(cut{c}.x)));
%                 fl_info(c,2) = mean(cut{c}.x(~isnan(cut{c}.x)));
%                 fl_info(c,3) = mean(cut{c}.r(~isnan(cut{c}.r)));
%             end
%         end
%         nan_info
%         fl_info

    end

end

if IJ == 1
    % Loop over each block and check intersection of the line of the
    % unstructured cut with the constant lines in the block
    for bid = bids

        % Check if the line lies in the block or not
        [flag n] = inblock(g,bid,xr,1);
        if flag == 1
            % Use matrix operations to determine all intersectin
            % coordinates
            [cut{bid+1}.ri cut{bid+1}.xi cut{bid+1}.Ind] = findintersection2(g{bid+1}.r,g{bid+1}.x,xr(n,[2 1]),2);

            % Select 4 bounding elements
            for v = 1:length(varnames)
                varname = varnames{v};
                var = cat(2,2*g{bid+1}.(varname)(:,1,:)-g{bid+1}.(varname)(:,2,:)...
                    ,g{bid+1}.(varname)(:,:,:), ...
                    2*g{bid+1}.(varname)(:,end,:)-g{bid+1}.(varname)(:,end-1,:));

                var = permute(var,[2 1 3]);
                cut{bid+1}.(varname) = reshape(var(cut{bid+1}.Ind==1),[4,g{bid+1}.attribute.ni,g{bid+1}.attribute.nk]);
                cut{bid+1}.(varname) = permute(cut{bid+1}.(varname),[2 3 1]);
            end

            % Interpolate using a spline
            for v = 1:length(interpnames)
                varname = interpnames{v};            
                cut{bid+1}.(varname) = spline1(cut{bid+1}.r, cut{bid+1}.(varname), cut{bid+1}.ri);
            end

            cut{bid+1}.x(:,:,1) = cut{bid+1}.xi;
            cut{bid+1}.r(:,:,1) = cut{bid+1}.ri;
            
            % Remove data
            for v = 1:length(varnames)
                varname = varnames{v};
                cut{bid+1}.(varname) = cut{bid+1}.(varname)(:,:,1);
            end            
            
            cut{bid+1} = rmfield(cut{bid+1},{'ri', 'Ind', 'xi'});
            
        end
    end
    
    % Remove empty arrays from cut
    for c = size(cut,1):-1:1
        if isempty(cut{c}) == 1
            cut(c) = [];
        end
    end
end

if size(cut,1) == 1 && iscell(cut)
    cut = cut{1};
end
    

% Turn everything back to singles
if iscell(cut) == 1
    for c = 1:size(cut,1)
        for v = 1:length(varnames)
            varname = varnames{v};
            cut{c}.(varname) = single(cut{c}.(varname));
        end
    end
else
    for v = 1:length(varnames)
        varname = varnames{v};
        if isfloat(cut.(varname))
            cut.(varname) = single(cut.(varname));
        end
    end
end
            
end

function [flag n] = inblock(g,bid,xr,type)
% Check if  line intersects block geometry and return intersecting coords
    flag = 0;
    
    xmin = min(min(min(g{bid+1}.x))); xmax = max(max(max(g{bid+1}.x)));
    rmin = min(min(min(g{bid+1}.r))); rmax = max(max(max(g{bid+1}.r)));
    
    % type 1 is a x1 -> x2 line - blade to blade slice
    if type == 1
        if size(xr,1) > 2
            % Reduce xr data based on what actual lies in the block
            n = find(((xr(:,1) < xmin) | (xr(:,1) > xmax)) == 0);
            if sum(xr(n,2) > rmin & xr(n,2) < rmax) > 0
                flag = 1;
                n = [n(1)-2 ; n(1)-1 ; n ; n(end)+1 ; n(end)+2];
                n(n<1) = [];
                n(n>size(xr,1)) = [];
            else
                n = [];
            end
        else
            % Don't reduce xr data as it is a straight line intersecting
            if max(xr(:,2)) > rmin && min(xr(:,2)) < rmax
                flag = 1;
                n = 1:size(xr,1);
            else
                n = [];
            end
        end
    end
    
    % type 2 is a r1 -> r2 line - traverse plane
    if type == 2
        if max(xr(:,1)) > xmin && min(xr(:,1)) < xmax
            flag = 1;
            n = 1:size(xr,1);
        else
            n = [];
        end
%         n = find(((xr(:,2) < rmin) | (xr(:,2) > rmax)) == 0);
%         if sum(xr(n,1) > xmin & xr(n,1) < xmax) > 0
%             flag = 1;
%             n = [n(1)-2 ; n(1)-1 ; n ; n(end)+1 ; n(end)+2];
%         else
%             n = [];
%         end
    end

end

function [x y ind] = findintersection2(x1,y1,xy2,dim)
% Always operate along 1st dimension, x1, y1 are permuted by dim if
% necessary
if dim == 2
    x1 = permute(x1,[2 1 3]);
    y1 = permute(y1,[2 1 3]);
end

x2 = xy2(:,1); y2 = xy2(:,2);
ind = zeros(size(x1,1)+2, size(x1,2), size(x1,3));

% Define absolute coordinates and differences for first lines
xp1 = repmat(x1(1:end-1,:,:),[1 1 1 size(x2,1)-1]);
yp1 = repmat(y1(1:end-1,:,:),[1 1 1 size(x2,1)-1]);

xs1 = repmat(x1(2:end,:,:)-x1(1:end-1,:,:),[1 1 1 size(x2,1)-1]);
ys1 = repmat(y1(2:end,:,:)-y1(1:end-1,:,:),[1 1 1 size(x2,1)-1]);

% Define absolute coordinates and differences for second repeated lines
xp2 = repmat(reshape(x2(1:end-1),1,1,1,[]),[size(xp1,1) size(x1,2) size(x1,3) 1]);
yp2 = repmat(reshape(y2(1:end-1),1,1,1,[]),[size(yp1,1) size(y1,2) size(y1,3) 1]);

xs2 = repmat(reshape(x2(2:end)-x2(1:end-1),1,1,1,[]),[size(xp1,1) size(x1,2) size(x1,3) 1]);
ys2 = repmat(reshape(y2(2:end)-y2(1:end-1),1,1,1,[]),[size(yp1,1) size(y1,2) size(y1,3) 1]);

% Determine scalar constants
t1 = ((xp2 - xp1).*ys2 - (yp2 - yp1).*xs2) ./ (xs1.*ys2 - ys1.*xs2);
t2 = ((xp1 - xp2).*ys1 - (yp1 - yp2).*xs1) ./ (xs2.*ys1 - ys2.*xs1);

% Remove all those that do not correspond to a mutual intersection
t1(t1<0) = 0; t1(t1>1) = 0;
t2(t2<0) = 0; t2(t2>1) = 0;

% Calcuate x,y and build indexing matrix
x = xp1 + t1.*xs1; y = yp2 + t2.*ys2;

t1(t1~=0) = 0.5;
t2(t2~=0) = 0.5;

t = t1+t2;
t(t~=1) = 0;

x = squeeze(sum(sum(x.*t,1),4));
y = squeeze(sum(sum(y.*t,1),4));

ind(1:end-3,:,:) = sum(t,4);

% Set index matrix to be a vector of length 4
I = find(ind);
I = reshape([I I+1 I+2 I+3],[],1);
ind(I) = 1;

end


function [yi] = spline1(x,y,xi)
% Take an array 4 deep and compute the function array y using first order
% accurate spline interpolation
A = (x(:,:,3) - xi) ./ (x(:,:,3) - x(:,:,2));
B = 1 - A;

C = 1/6 * (A.*A.*A - A) .* ((x(:,:,3) - x(:,:,2)).^2);
D = 1/6 * (B.*B.*B - B) .* ((x(:,:,3) - x(:,:,2)).^2);

y1 = y(:,:,2); y2 = y(:,:,3);

top = (x(:,:,2) - x(:,:,1)).*y(:,:,3) - (x(:,:,3) - x(:,:,1)).*y(:,:,2) + (x(:,:,3) - x(:,:,2)).*y(:,:,1);
bottom = 0.5*(x(:,:,3) - x(:,:,2)).*(x(:,:,2)-x(:,:,1)).*(x(:,:,3)-x(:,:,1));
d2y1 = top./bottom;

top = (x(:,:,3) - x(:,:,2)).*y(:,:,4) - (x(:,:,4) - x(:,:,2)).*y(:,:,3) + (x(:,:,4) - x(:,:,3)).*y(:,:,2);
bottom = 0.5*(x(:,:,4) - x(:,:,3)).*(x(:,:,3)-x(:,:,2)).*(x(:,:,4)-x(:,:,2));
d2y2 = top./bottom;

yi = A.*y1 + B.*y2 + C.*d2y1 + D.*d2y2;

end
     
        
function [x r ii] = findintersection(x_block,r_block,x0,x1,r0,r1)
    % Compute the location of the intersection between two lines, one
    % defined by a number of points the other straight between two points
    r_line = r_block;
    x_line = x0 + (r_line - r0) / (r1 - r0) * (x1 - x0);

    x_diff = x_block-x_line;
    if x_diff(end) < 0 || x_diff(1) > 0
        ii = [];
        x = [];
        r = [];
        return
    end
        
    [~, ii(1)] = min(abs(x_diff));
    
    if x_diff(ii(1)) > 0
        ii(2) = ii(1)-1;
    else
        ii(2) = ii(1)+1;
    end

    ii = sort(ii);

    xv = x1;
    rv = r1;
    xu = x0;
    ru = r0;
    xb = x_block(ii(2));
    rb = r_block(ii(2));
    xa = x_block(ii(1));
    ra = r_block(ii(1));    
    
    f = ((xv-xu)*(ra-ru) - (rv-ru)*(xa-xu)) / ((rv-ru)*(xb-xa) - (xv-xu)*(rb-ra));
    x = xa + f*(xb-xa);
    r = ra + f*(rb-ra);
%     figure()
%     hold on
%     plot([xu xv],[ru rv],'k-')
%     plot([xa xb],[ra rb],'r-')
%     plot(x,r,'k.')
%     [xa xb ra rb]
%     axis([xa xb ra rb])
%     
%     gtg
end