function zi = interp_2D(X,Y,Zs,xi,yi)
% Interpolate flow parameters from a general non-uniform grid onto a new 2D grid

% Grid orientation
i = [X(2,1) Y(2,1)] - [X(1,1) Y(1,1)];
j = [X(1,2) Y(1,2)] - [X(1,1) Y(1,1)];
z = i(1) * j(2) - j(1) * i(2);
if z < 0
    temp = X; X = Y; Y = temp;
    temp = xi; xi = yi; yi = temp;
end

% Sizes of grid
ni_d = size(X,1); nj_d = size(X,2); 
ni_t = size(xi,1); nj_t = size(xi,2); 

% Preallocate primary variables for target
zi = zeros(ni_t,nj_t,size(Zs,3));

% Construct line segment vectors for domain cells
xy_d = cat(3,X,Y);
l = cat(4,xy_d(2:end,1:end-1,:) - xy_d(1:end-1,1:end-1,:),  xy_d(2:end,2:end,:) - xy_d(2:end,1:end-1,:),...
    xy_d(1:end-1,2:end,:) - xy_d(2:end,2:end,:),  xy_d(1:end-1,1:end-1,:) - xy_d(1:end-1,2:end,:));

xy_t = cat(3,xi,yi);
c = zeros(ni_d-1,nj_d-1,2,4); z_lc = zeros(ni_d-1,nj_d-1,1,4);

for i = 1:ni_t
    for j = 1:nj_t
        % Vectors from points to cell corners
        xy = repmat(xy_t(i,j,:),[ni_d-1 nj_d-1]);
        c(:,:,:,1) = xy - xy_d(1:end-1,1:end-1,:);
        c(:,:,:,2) = xy - xy_d(2:end,1:end-1,:);
        c(:,:,:,3) = xy - xy_d(2:end,2:end,:);
        c(:,:,:,4) = xy - xy_d(1:end-1,2:end,:);
        
        % Z component of cross product between a and l
        z_lc(:,:,1,:) = l(:,:,1,:) .* c(:,:,2,:) - l(:,:,2,:) .* c(:,:,1,:);
        q1 = sum(sign(z_lc),4) == 4;
%         [iq, jq] = find(sum(sign(z_lc),4) == 4);
        
        if sum(sum(q1)) == 1 %length(iq) == 1
            % Extract cell coordinates
            q = zeros(ni_d,nj_d);
            q(1:end-1,1:end-1) = q(1:end-1,1:end-1) + q1;
            q(2:end,1:end-1) = q(2:end,1:end-1) + q1;
            q(2:end,2:end) = q(2:end,2:end) + q1;
            q(1:end-1,2:end) = q(1:end-1,2:end) + q1;

            % Bilinear mapping function
            xy_q = reshape(xy_d(repmat(q,[1 1 2]) > 0),[4 2]); xy_q = xy_q([1 2 4 3],:);
%             xy_q = xy_d(iq:iq+1,jq:jq+1,:); xy_q = xy_q([1 2 4 3],:);

            A = [1 0 0 0 ; -1 1 0 0 ; -1 0 0 1 ; 1 -1 1 -1];
            a = A * xy_q(:,1);
            b = A * xy_q(:,2);

            % Logical coordinates, s,r describe position in square
            aa = a(4)*b(3) - a(3)*b(4);
            bb = a(4)*b(1) - a(1)*b(4) + a(2)*b(3) - a(3)*b(2) + xy_t(i,j,1)*b(4) - xy_t(i,j,2)*a(4);
            cc = a(2)*b(1) - a(1)*b(2) + xy_t(i,j,1)*b(2) - xy_t(i,j,2)*a(2);
            if aa ~= 0 
                % Solve quadratic equation
                det = sqrt(bb*bb - 4*aa*cc);
                s = (-bb+det)/(2*aa);

                r = (xy_t(i,j,1)-a(1)-a(3)*s)/(a(2)+a(4)*s);
            else
                % Solve linear equation
                s = -cc / bb;
                r = (xy_t(i,j,1)-a(1)-a(3)*s)/(a(2)+a(4)*s);
            end
            
            % Cell corner flow values
            p_q = reshape(Zs(repmat(q,[1 1 size(Zs,3)]) > 0),[4 size(Zs,3)]); 
            p_q = p_q([1 2 4 3],:);
            
            % Bilinear interpolation on target point
            zi(i,j,:) = (1-r)*(1-s)*p_q(1,:) + r*(1-s)*p_q(2,:) + r*s*p_q(3,:) + (1-r)*s*p_q(4,:);
        else
            % Use NaNs for points outside of the domain.
            zi(i,j,:) = nan;
        end
        
    end
end


end