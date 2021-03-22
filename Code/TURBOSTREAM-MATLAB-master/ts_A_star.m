function [A path] = ts_A_star(g, u_bump, u_target, neighbour, A, max_dist, plot_path)
% A* algorithm for calculating path through a surface

if exist('plot_path','var') == 0
    plot_path = 0;
end

nblade = double(g{u_bump(1)+1}.bv.nblade);

%% A* ALGORITHM FOR SHORTEST PATH DISTANCE

% First Point variables
u = u_bump;
xrrt_u_target = reshape(g{u_target(1)+1}.xrrt(u_target(2)+1,u_target(3)+1,u_target(4)+1,:),3,1);
Q = 0; e = 0; h = pythag(reshape(g{u(1)+1}.xrrt(u(2)+1,u(3)+1,u(4)+1,:),3,1), xrrt_u_target);
f = e + h; bid_prev = u(1); I_prev = u(2); J_prev = u(3); K_prev = u(4);

A{u(1)+1}(u(2)+1,u(3)+1,u(4)+1,:) = [Q e h f bid_prev I_prev J_prev K_prev];

% Initialise queue [ bid I J K f ] as those that have been visited already
queue = [];
for bid = ts_get_bids(g)
    n = find(A{bid+1}(:,:,:,1) == 0);
    [I,J,K] = ind2sub(size(A{bid+1}(:,:,:,1)), n);
    if isempty(I) == 0
        F = A{bid+1}(:,:,:,4);
        queue = [queue ; repmat(bid,length(I),1) I-1 J-1 K-1 F(n)];
    end
end

% Check queue is not empty
while isempty(queue) == 0
    % Sort queue and choose node with lowest f
    queue = sortrows(queue,size(queue,2));
    u = int16(queue(1,1:4));
    
    % Check if at target point and quit function
    if isequal(u, u_target) == 1
        A{u(1)+1}(u(2)+1,u(3)+1,u(4)+1,1) = 0;
        %disp(A{u(1),u(2)}(u(3),u(4),2));
        break
    elseif A{u(1)+1}(u(2)+1,u(3)+1,u(4)+1,2) > max_dist
        A{u_target(1)+1}(u_target(2)+1,u_target(3)+1,u_target(4)+1,2) = inf;
        path = [];
        return
    end
    
    % Get neighbouring nodes for inspection if point isn't already visited
    [v_neighbours] = neighbour{u(1)+1}{u(2)+1,u(3)+1,u(4)+1};
    xrrt_u = reshape(g{u(1)+1}.xrrt(u(2)+1,u(3)+1,u(4)+1,:),3,1);

    % Loop over all given neighbouring co-ordinates
    for n = 1:size(v_neighbours,1)
        v = v_neighbours(n,:);

        if A{v(1)+1}(v(2)+1,v(3)+1,v(4)+1,1) == 1 % If Q is unvisited
            % Calculate distance between neighbour and current visited node
            % Use cylindrical polar coords for correct surface distances
            % Correct for the pitch distance to allow crossing of
            % periodic walls
            xrrt_v = reshape(g{v(1)+1}.xrrt(v(2)+1,v(3)+1,v(4)+1,:),3,1);

            if abs(xrrt_v(3) - xrrt_u(3)) > xrrt_u(2) * pi / nblade
                if xrrt_v(3) > xrrt_u(3)
                    xrrt_v(3) = xrrt_v(3) - (2*pi)/nblade * xrrt_v(2);
                else
                    xrrt_v(3) = xrrt_v(3) + (2*pi)/nblade * xrrt_v(2);
                end
            end

            alt = A{u(1)+1}(u(2)+1,u(3)+1,u(4)+1,2) + pythag(xrrt_u, xrrt_v);

            % Check if e, distance function is less than that already stored
            if alt < A{v(1)+1}(v(2)+1,v(3)+1,v(4)+1,2)
                A{v(1)+1}(v(2)+1,v(3)+1,v(4)+1,2) = alt;

                % Calculate heuristic function - shortest straight distance
                % squared
                A{v(1)+1}(v(2)+1,v(3)+1,v(4)+1,3) = sum((xrrt_v - xrrt_u_target).^2)^0.5;

                % Calculate total function and add neighbour to queue
                A{v(1)+1}(v(2)+1,v(3)+1,v(4)+1,4) = A{v(1)+1}(v(2)+1,v(3)+1,v(4)+1,2) + A{v(1)+1}(v(2)+1,v(3)+1,v(4)+1,3);
                queue = [queue ; double(v) A{v(1)+1}(v(2)+1,v(3)+1,v(4)+1,4)];

                % Record previous node indices
                A{v(1)+1}(v(2)+1,v(3)+1,v(4)+1,5:8) = u;
            end
        end
    end
   
    % Set node to visited and remove from the queue
    A{u(1)+1}(u(2)+1,u(3)+1,u(4)+1,1) = 0;
    queue(1,:) = [];
end

%% BACKTRACK FROM TARGET TO FIND AND PLOT PATH
if plot_path == 1
    % Create first point (now work in cartesian for accurate plotting)
    u = u_target;
    path = [reshape(g{u(1)+1}.xrrt(u(2)+1,u(3)+1,u(4)+1,:),3,1).' u];

    % Track back through previous datasets to create list of jumps
    while isequal(u,u_bump) == 0
        u = reshape(A{u(1)+1}(u(2)+1,u(3)+1,u(4)+1,5:8), 1, 4);
        path = [path ; reshape(g{u(1)+1}.xrrt(u(2)+1,u(3)+1,u(4)+1,:),3,1).' u];
    end
    
    C = cell(size(g));
    for bid = ts_get_bids(g)
        % Highlight visited nodes in colour matrix
        Q = A{bid+1}(:,:,:,1);

        C{bid+1} = zeros([size(g{bid+1}.x) 3]);
        [I,J,K] = ind2sub(size(Q), find(Q == 0));
        if isempty(I) == 0
            for n = 1:length(I)
                C{bid+1}(I(n),J(n),K(n),:) = [0.7 0.5 0.5];
            end
        end
        
    end
    ts_plot_surface(g, C, 1, 0)
        
    % Plot on path
    d = plot3(path(:,1), path(:,2), path(:,3), 'r-');
    set(d,'LineWidth',2);

else
    path = 'Path has not be created for reduced runtime';
end

end