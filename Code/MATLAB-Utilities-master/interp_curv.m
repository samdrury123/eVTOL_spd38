function xyi = interp_curv(s,xy,ds,lim,seg_type,plot_stuff)
% INTERP_CURV  Interpolate a function in segments based upon curvature limits
%
%   xyi = INTERP_CURV(s,xy,ds,lim,seg_type,plot_stuff)
%
%   s - distance through curve points
%   xy - coordinates of curve
%   ds - desired resolution in curve distance, or actual desired values
%   lim - curvature limit, set just below the highest values
%   seg_type - 'cont' or 'seg' for different type of evaluation
%   plot_stuff - 0 or 1 to show working

% Default to continuous interpolation through segments
if exist('seg_type','var') == 0
    seg_type = 'cont';
end

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Calculate curvature on polyline
k = curvature(xy(:,1),xy(:,2)); k(isnan(k)) = 0;

% Plot curvature
if plot_stuff == 1
    figure(); hold on; grid on; box on;
    plot(s,k,'.-')
end

% Find segment bounding indices
i = [find((k(1:end-1)-lim) .* (k(2:end) - lim) < 0) ; length(k)]; 
i_start = find(k(2:end) >= lim & k(1:end-1) < lim);
i_end = find(k(1:end-1) >= lim & k(2:end) < lim) + 1;

% Remove points close to each other to deal with sharp corners
for n = 1:length(i_start)
    for m = length(i_end):-1:1
        if i_end(m) == i_start(n) + 2
            i_end(m) = [];
            i_start(n) = i_start(n) + 1;
        end
    end
end

% Combine bounding indices and append start an end points
i = sort([1 ; i_start ; i_end ; length(s)]);

% Plot bounding indices
if plot_stuff == 1
    plot(s(i),k(i),'o');
end

% Define resolution to interpolate to
if strcmp(seg_type,'cont') == 1
    
    % Spacing continuous across all segments
    if numel(ds) == 1
        si = linspace(0,s(end),round(s(end)/ds))';
    else
        si = ds;
    end
    
elseif strcmp(seg_type,'seg') == 1
    
    % Loop over all segments and define partial interpolation functions
    si = [];
    for n = 1:length(i)-1
        ni = round((s(i(n+1)) - s(i(n))) / ds);
        si = [si ; linspace(s(i(n)),s(i(n+1)),ni)'];
    end
    
    % Take unique points only
    si = unique(si);
end

% Loop over all segments and interpolate
xyi = zeros(length(si),2);
for n = 1:length(i)-1

    % Determine indices in current segment
    q = si <= s(i(n+1)) & si >= s(i(n));

    % Interpolate function
    for m = 1:2
        xyi(q,m) = interp1(s(i(n):i(n+1)),xy(i(n):i(n+1),m),si(q),'spline');
    end
end

% Recalculate and plot curvature
if plot_stuff == 1
    ki = curvature(xyi(:,1),xyi(:,2));
    plot(si,ki,'.-');
end

end