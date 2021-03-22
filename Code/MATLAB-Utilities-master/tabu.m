function [x,Q] = tabu(g,x,d,lb,ub,options)
% TABU  Optimise an objective function with many degrees of freedom
%   
%   [x,Q] = TABU_SEARCH(f,x,lb,ub,options)

% Reset random number lists
rng('default');

% Initialise counters
i = 1; n_calls = 1; n_improve = 0; d_scale = 1; i_diversify = 1;

% Record size of input matrix and reshape into a list
sx = size(x); nx = numel(x); x = reshape(x,nx,1); 
nj = size(d,ndims(d)); d = reshape(d,nx,nj);
lb = reshape(lb,nx,1); ub = reshape(ub,nx,1); 

% Initialise banned move arrays
n_ban = zeros(1,nj); q_ban = zeros(1,nj);

% Plot contours of the function if requested
if options.plot_func == 1
    
    % Define range to plot over
    np = 100; 
    x1 = repmat(linspace(lb(1),ub(1),np),[np 1]); 
    x2 = repmat(linspace(lb(2),ub(2),np)',[1 np]); 
    x_func = [reshape(x1,1,[]) ; reshape(x2,1,[])];
    
    % Evaluate function
    f = reshape(g(x_func),[np np]);
    
    % Evaluate contour levels
    v = [min(f(:)) quantile(f(:),18) max(f(:))];
    
    % Plot contours of objective function
    h.func = figure('units','normalized','position',[0.5000 0.1 0.4000 0.7]); 
    hold on; grid on; box on; axis equal;
    xlabel('x_1'); ylabel('x_2');
    contour(x1,x2,f,v);
    
end

% Function value at the start point
f = g(reshape(x,sx)); 

% Start recording coordinates and function values
f_iter(i) = f; x_iter(:,i) = x; n_calls_iter(i) = n_calls;

% Open a figure window to plot progress of optimiser
if options.plot_stuff == 1
    
    % Open window
    h.progress = figure('units','normalized','position',[0.02 0.4 0.45 0.4]);
    
    % Plot input values
    cols = jet(nx);
    subplot(3,1,1); hold on; grid on; box on; ylabel('Input');
    for a = 1:nx; l.x{a} = plot(x_iter(a,:)','.-','color',cols(a,:)); end;
    
    % Plot output values
    subplot(3,1,2); hold on; grid on; box on; ylabel('Output');
    l.f = plot(f_iter,'k.-');
    
    % Plot function call count
    subplot(3,1,3); hold on; grid on; box on; ylabel('Function Count'); xlabel('Iteration');
    l.m = plot(n_calls_iter,'k.-');
    
end

% Plot the path of the function so far
if options.plot_func == 1
    figure(h.func); l.p = plot(x_iter(1,:),x_iter(2,:),'k.-');
end

% Initialise matrix to store all move data
Q = nan(1,nx+2,nj);

% Start iterating on the input parameters
d_max = max(abs(d(:)));
while i < options.n_max && d_scale * d_max >= options.d_lim

    % Get the last point
    x = x_iter(:,end); f = f_iter(end);

    % Apply the actions
    x_batch = nan(nx,nj); f_batch = inf(1,nj);
    for j = 1:nj

        % Apply the move
        x_temp = x + d_scale * d(:,j);

        % Check the limits
        if any(x_temp < lb) || any(x_temp > ub)
            continue
        end

        % Check the point is not near previous ones
        q = pdist2(x_iter',x_temp');
        if min(q) < 0.4 * d_scale * d_max
            continue
        end

        % Check the current move is not banned
        if q_ban(j) > 0
            continue
        end
        
        % Record the move
        x_batch(:,j) = x_temp;
        
    end

    % Check all the function values
    j = any(isnan(x_batch) == 0,1); ni_batch = sum(j);
    f_batch(j) = g(reshape(x_batch(:,j),[sx ni_batch])); 
%     j = find(j); for i = 1:ni_batch; f_batch(j(i)) = g(x_batch(:,j(i))); end;
    n_calls = n_calls + ni_batch;

    % Sort the action list by objective function value
    [~,j] = sort(f_batch); 
    d = d(:,j); q_ban = q_ban(j); n_ban = n_ban(j);
    x_batch = x_batch(:,j); f_batch = f_batch(j);

    % Sort the names and print the list if specified
    if isfield(options,'d_name') == 1
        for j = 1:nj; fprintf('%s\t%g\n',options.d_name{j},f_batch(j)); end;
        fprintf('\n');
    end
    
    % Save the sorted batch data in master array
    Q(i,1,:) = i; Q(i,2:nx+1,:) = x_batch; Q(i,end,:) = f_batch;
    
    % Modify improvement counter
    if all(isinf(f_batch) == 1) == 1
        % Set to intensify if no valid moves exist
        n_improve = options.improve;
    elseif f_batch(1) > f
        % Increase intensification counter if no moves are better than the current one
        n_improve = n_improve + 1;
    else
        % Reset if found a better move
        n_improve = 0;
    end
    
    % Ban moves if they repeatedly produce bad results
    q_ban = max(q_ban - 1,0);
    q = f_batch > f & q_ban == 0; n_ban(q) = n_ban(q) + 1;
    n_ban(q == 0) = 0; q_ban(n_ban >= options.improve) = options.ban;

    % Remove bans from the top half of the move list
    q_ban(1:round(nj/2)) = 0;
    
    % Apply the best (or least worst) move
    x = x_batch(:,1); f = f_batch(1);

    % Increase iteration count
    i = i + 1;

    % Record the current point
    f_iter(i) = f; x_iter(:,i) = x; n_calls_iter(i) = n_calls;
    
    % Find the best point since the last diversification
    [~,i_min] = min(f_iter(i_diversify(end):end)); i_min = i_min + i_diversify(end) - 1;
    n_abandon = (i - i_min);
    
    % Check status of current point relative to best so far and decide whether to diversify
    if n_improve >= options.improve || n_abandon >= options.abandon
        
        % Increase iteration count
        i = i + 1; 
        
        % Diversify or intensify depending on desired number of diversifications
        if length(i_diversify) <= options.n_diversify

            % Evaluate a new random point
            x_iter(:,i) = d_max * round((rand(size(lb)) .* (ub - lb) + lb) / d_max); 
            f_iter(i) = g(reshape(x_iter(:,end),sx)); 

            % Modify counters
            n_calls_iter(i) = n_calls + 1; i_diversify = [i_diversify i];
        else
            
            % Find the global best point
            [~,i_min] = min(f_iter);
            
            % Half step size and jump back to the best point so far
            d_scale = 0.5 * d_scale;
            x_iter(:,i) = x_iter(:,i_min); 
            
            % Record the loss for the best point so far
            n_calls_iter(i) = n_calls; f_iter(i) = 0.999999999 * f_iter(i_min);
        end
        
        % Save the new point in the master array
        Q(i,1,:) = i; Q(i,2:nx+1,1) = x_iter(:,i); Q(i,end,1) = f_iter(i);
        
        % Remove all moves from ban list
        q_ban(:) = 0;
        
        % Reset improvement counter
        n_improve = 0;
    end

    % Display counters showing tabu memory
    fprintf('Iteration %i\tImprove %i\tAbandon %i\tScale %f\n\n',[i n_improve n_abandon d_scale]); 
    
    % Update the convergence plots
    if options.plot_stuff == 1
        figure(h.progress);
        for a = 1:nx; set(l.x{a},'xdata',1:i,'ydata',x_iter(a,:)); end;
        set(l.f,'xdata',1:i,'ydata',f_iter); set(l.m,'xdata',1:i,'ydata',n_calls_iter);
        drawnow;
    end

    % Update the path and plot evaluated points
    if options.plot_func == 1
        q = isinf(f_batch) == 0;
        figure(h.func); plot(x_batch(1,q),x_batch(2,q),'r.');
        set(l.p,'xdata',x_iter(1,:),'ydata',x_iter(2,:))
        drawnow;
    end
    
end

% Return the inputs for the best result found in the original matrix size
[~,i_min] = min(f_iter); x = reshape(x_iter(:,i_min),sx);

% Plot the best result
if options.plot_stuff == 1
    figure(h.func); plot(x(1),x(2),'ko','markerfacecolor','k');
end

end

