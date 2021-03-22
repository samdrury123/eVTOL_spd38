function g_target = ts_interp_2D(g_domain,g_target,x_name,y_name)
% Interpolate flow parameters onto a new 2D grid

if exist('x_name','var') == 0
    x_name = 'x';
end
if exist('y_name','var') == 0
    y_name = 'rt';
end

varnames = {'x' 'r' 'rt' 'ro' 'rovx' 'rovr' 'rorvt' 'roe'};

% Double precision
for v = 1:length(varnames)
    varname = varnames{v};
    g_domain.(varname) = double(g_domain.(varname));
%     g_target.(varname) = double(g_target.(varname));
end

% Repeat domain by +/- a pitch if traverse plane is required.
if strcmp(x_name,'r') == 1 && strcmp(y_name,'rt') == 1
    g_domain.rt = [g_domain.rt(:,1:end-1) - (2*pi*g_domain.r(:,1:end-1)) / double(g_domain.bv.nblade) ...
        g_domain.rt g_domain.rt(:,2:end) + (2*pi*g_domain.r(:,2:end)) / double(g_domain.bv.nblade)];
    for v = 1:length(varnames)
        varname = varnames{v};
        if strcmp(varname,'rt') == 0
            g_domain.(varname) = [g_domain.(varname)(:,1:end-1) g_domain.(varname)...
                g_domain.(varname)(:,2:end)];

        end
    end
end

% Grid orientation
i = [g_domain.(x_name)(2,1) g_domain.(y_name)(2,1)] - [g_domain.(x_name)(1,1) g_domain.(y_name)(1,1)];
j = [g_domain.(x_name)(1,2) g_domain.(y_name)(1,2)] - [g_domain.(x_name)(1,1) g_domain.(y_name)(1,1)];
z = i(1) * j(2) - j(1) * i(2);
if z < 0
    temp = x_name; x_name = y_name; y_name = temp;
end

ni_d = size(g_domain.(x_name),1); nj_d = size(g_domain.(x_name),2); 
ni_t = size(g_target.(x_name),1); nj_t = size(g_target.(x_name),2); 

% Pull out primary variables for domain
p_d = zeros(ni_d,nj_d,length(varnames));
for v = 1:length(varnames)
    varname = varnames{v};
    p_d(:,:,v) = g_domain.(varname);
end

% Preallocate primary variables for target
o = zeros(ni_t,nj_t,length(varnames));

% Construct line segment vectors for domain cells
xy_d = cat(3,g_domain.(x_name),g_domain.(y_name));
l = cat(4,xy_d(2:end,1:end-1,:) - xy_d(1:end-1,1:end-1,:),  xy_d(2:end,2:end,:) - xy_d(2:end,1:end-1,:),...
    xy_d(1:end-1,2:end,:) - xy_d(2:end,2:end,:),  xy_d(1:end-1,1:end-1,:) - xy_d(1:end-1,2:end,:));

xy_t = cat(3,g_target.(x_name),g_target.(y_name));
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
            p_q = reshape(p_d(repmat(q,[1 1 length(varnames)]) > 0),[4 length(varnames)]); 
            p_q = p_q([1 2 4 3],:);
            
            % Bilinear interpolation on target point
            o(i,j,:) = (1-r)*(1-s)*p_q(1,:) + r*(1-s)*p_q(2,:) + r*s*p_q(3,:) + (1-r)*s*p_q(4,:);
        else
            % Use NaNs for points outside of the domain.
            o(i,j,:) = nan;
        end
        
    end
end

% Put interpolated flow back into data structure
for v = 1:length(varnames)
    varname = varnames{v};
    if strcmp(varname,x_name) == 0 && strcmp(varname,y_name) == 0
        g_target.(varname) = o(:,:,v);
    end
end

end