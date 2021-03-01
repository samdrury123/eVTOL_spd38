function y_new = diff_smooth(x,y,y_lim)
% Smooth a function by downsampling based upon differences in value

% Size of vectors
ni = length(x);

% Downsample based upon limit value
x_set = x(1); y_set = y(1);
i_last = 2;
while i_last < ni
    
    % Calculate difference between last value and the remainder of the list
    dy = abs(y_set(end) - y(i_last:end));
    
    % Find the next value to cross the limit
    i_next = find(dy > y_lim,1) + i_last - 1;
    
    % Record the next values
    x_set = [x_set ; x(i_next)];
    y_set = [y_set ; y(i_next)];
    
    % Update index
    i_last = i_next;
end

% Record the final values as well as the first
x_set = [x_set ; x(end)]; y_set = [y_set ; y(end)];

% Interpolate new values from the set
[x_set,i] = unique(x_set); y_set = y_set(i);
y_new = interp1(x_set,y_set,x,'pchip');

end