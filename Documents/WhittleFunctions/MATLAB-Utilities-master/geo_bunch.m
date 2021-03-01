function [s,r] = geo_bunch(ni,s1,s3,ds1,ds3,f)
% GEO_BUNCH  Geometric point bunching to a specified spacing and value

% Limits on common ratios
r_min = 1e-5; r_max = 1e5; r_old = 1; r_tol = 1e-9;

% Determine type of spacings to use
if exist('ds3','var') == 0 || isempty(ds3) == 1
    % Generate single sided geometric progression
    
    % Bisect common ratio until converged
    r = 0.5 * (r_min + r_max);
    while abs(r - r_old) > r_tol
        
        % Evaluate the sum of the series
        S1 = ds1 * (1 - r^(ni-1)) / (1 - r);
        
        % Check the sum and move the limits
        if S1 > abs(s3 - s1)
            r_max = r;
        else
            r_min = r;
        end
        
        % Update current ratio
        r_old = r; r = 0.5 * (r_min + r_max);
        
    end
    
    % Generate series
    ds = zeros(1,ni-1); ds(1) = ds1; for i = 2:ni-1; ds(i) = ds(i-1) * r; end;
    s = [0 cumsum(ds)] + s1; 
    
elseif exist('f','var') == 0 || isempty(f) == 1
    % Generate two sided geometric progression with minimum common ratio
    
    % Divide number of points between both parts of the series
    ni1 = round(ni/2); ni3 = ni - ni1 + 1;
    
    % Bisect common ratio until converged
    r = 0.5 * (r_min + r_max);
    while abs(r - r_old) > r_tol
        
        % Evaluate the sum of the series
        S1 = ds1 * (1 - r^(ni1-1)) / (1 - r);
        S3 = ds3 * (1 - r^(ni3-1)) / (1 - r);
        
        % Check the sum and move the limits
        if S1 + S3 > abs(s3 - s1)
            r_max = r;
        else
            r_min = r;
        end
        
        % Update current ratio
        r_old = r; r = 0.5 * (r_min + r_max);
        
    end
    
    % Generate geometric differences
    ds1 = [ds1 zeros(1,ni1-2)]; for i = 2:ni1-1; ds1(i) = ds1(i-1) * r; end;
    ds3 = [ds3 zeros(1,ni3-2)]; for i = 2:ni3-1; ds3(i) = ds3(i-1) * r; end;
    s = [0 cumsum([ds1 flip(ds3,2)])] + s1; 
    
else
    % Generated two sided geometric points with specified nubmer of linearly spaced points in middle
    
    % Limits on common ratios
    r_min = 1e-5; r_max = 1e5; r_old = 1; r_tol = 1e-9;

    % Set number of constant difference points
    ni2 = round(f * ni); 

    % Bisect common ratio until converged
    r = 0.5 * (r_min + r_max);
    while abs(r - r_old) > r_tol

        % Calculate number of points in geometric regions
        ni3 = round((ni - ni2 + 2) / (1 + log(r * ds3) / log(r * ds1)));
        ni1 = ni - ni2 - ni3 + 2;

        % Evaluate the sum of the geometric series
        S1 = ds1 * (1 - r^ni1) / (1 - r);
        S3 = ds3 * (1 - r^ni3) / (1 - r);

        % Evaluate the sum of the linear series
        ds2a = ds1 * r^ni1; ds2b = ds3 * r^ni3;
        ds2 = linspace(ds2a,ds2b,ni2-1); S2 = sum(ds2(2:end-1));

        % Check the sum and move the limits
        if S1 + S2 + S3 > abs(s3 - s1)
            r_max = r;
        else
            r_min = r;
        end

        % Update current ratio
        r_old = r; r = 0.5 * (r_min + r_max);

    end

    % Generate geometric differences
    ds1 = [ds1 zeros(1,ni1-1)]; for i = 2:ni1; ds1(i) = ds1(i-1) * r; end;
    ds3 = [ds3 zeros(1,ni3-1)]; for i = 2:ni3; ds3(i) = ds3(i-1) * r; end;

    % Generate linear differences
    ds2 = linspace(ds1(end)*r,ds3(end)*r,ni2-1); ds2 = ds2(2:end-1);

    % Coordinates of points
    s = [0 cumsum([ds1 ds2 flip(ds3,2)])] + s1; 

end

end

