function t = bl_guess_camber(xy,plot_stuff,ilete)
% BL_GUESS_CAMBER  Determine a quick camberline from a blade section
%
%   t = BL_GUESS_CAMBER(xy,plot_stuff)
%
%   xy - 2D array of section coordinates [x1 y1 ; x2 xy2 ; ...
%   plot_stuff - 0 or 1 for showing working
%   ilete - optional argument to specify leading and trailing edge indices
%   t - data structure of output parameters

% Default plotting options
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Open figure window
if plot_stuff == 1
    figure(); hold on; axis equal; grid on; box on;
end

% Record original geometry
xy_orig = xy;

% Find le and te by max and min x
[~,i_le] = min(xy(:,1)); [~,i_te] = max(xy(:,1)); 

% Flip to point nose down
if xy(i_te,2) < xy(i_le,2)
    xy(:,2) = - xy(:,2);
    flip_blade = 1;
else
    flip_blade = 0;
end

% Plot blade section
if plot_stuff == 1
    plot(xy(:,1),xy(:,2),'b-');
end

% Determine whether to find leading and trailing edge locations or not
if exist('ilete','var') ~= 0
    i_le = ilete(1); i_te = ilete(2);
end

% Rotate blade by first guess of chord line
gamma = atand((xy(i_te,2) - xy(i_le,2)) / (xy(i_te,1) - xy(i_le,1)));
xy = [xy(:,1) * cosd(-gamma) - xy(:,2) * sind(-gamma) ...
    xy(:,1) * sind(-gamma) + xy(:,2) * cosd(-gamma)];

% Plot blade section
if plot_stuff == 1
    plot(xy(:,1),xy(:,2),'r-');
end

% Find new estimates of leading and trailing edges
if exist('ilete','var') == 0
    [~,i_le] = min(xy(:,1)); [~,i_te] = max(xy(:,1)); 
end

% Max and min coordinates
x_min = xy(i_le,1); x_max = xy(i_te,1);

% Split blade surfaces
ni = size(xy,1);
if i_le > i_te
    i_1 = (i_le:-1:i_te)';
    i_2 = [i_le:ni 2:i_te]';
else
    i_1 = [i_le:-1:2 ni:-1:i_te]';
    i_2 = (i_le:i_te)';        
end

% Calculate non-dimensional blade distances
xy_1 = xy(i_1,:); xy_2 = xy(i_2,:);
[~,i] = unique(xy_1,'rows'); xy_1 = xy_1(i,:);
[~,i] = unique(xy_2,'rows'); xy_2 = xy_2(i,:);
s_1 = [0 ; cumsum(sum(diff(xy_1,1,1).^2,2).^0.5)]; s_1 = s_1 / max(s_1);
s_2 = [0 ; cumsum(sum(diff(xy_2,1,1).^2,2).^0.5)]; s_2 = s_2 / max(s_2);

% Interpolate a few points on each surface
s_fit = hyperbolic_bunch(12,0.03,0.01)';
xy_1 = [interp1(s_1,xy_1(:,1),s_fit) interp1(s_1,xy_1(:,2),s_fit)];
xy_2 = [interp1(s_2,xy_2(:,1),s_fit) interp1(s_2,xy_2(:,2),s_fit)];

% Find midpoints between two surfaces
xy_cl = 0.5 * (xy_1 + xy_2); x = linspace(x_min,x_max,100)';

% Fit a spline camberline
spline_cl = spap2(2,4,xy_cl(:,1),xy_cl(:,2));
xy_cl = [x fnval(fnxtr(spline_cl,4),x)];

% Fit a polynomial camberline
% p = polyfit(xy_cl(2:end-1,1),xy_cl(2:end-1,2),4);
% xy_cl = [x polyval(p,x)];

% Fit a circular camberline
% p = CircleFitByTaubin(xy_cl);
% th = atand((xy_cl(:,1) - p(1)) ./ (xy_cl(:,2) - p(2)));
% th = interp1(xy_cl(:,1),th,x);
% xy_cl = [p(1) + p(3) * sind(th) p(2) + p(3) * cosd(th)];

% Plot guessed camberline
if plot_stuff == 1
    plot(xy_cl(:,1),xy_cl(:,2),'b-');
end

% Find new estimates of leading and trailing edges
[~,i_le] = min(abs(sum((xy - repmat(xy_cl(1,:),[ni 1])).^2,2)));
[~,i_te] = min(abs(sum((xy - repmat(xy_cl(end,:),[ni 1])).^2,2)));

% Split blade surfaces
if i_le > i_te
    i_1 = (i_le:-1:i_te)';
    i_2 = [i_le:ni 2:i_te]';
else
    i_1 = [i_le:-1:2 ni:-1:i_te]';
    i_2 = (i_le:i_te)';        
end

% Calculate non-dimensional blade distances
xy_1 = xy_orig(i_1,:); xy_2 = xy_orig(i_2,:);
s_1 = [0 ; cumsum(sum(diff(xy_1,1,1).^2,2).^0.5)]; s_1 = s_1 / max(s_1);
s_2 = [0 ; cumsum(sum(diff(xy_2,1,1).^2,2).^0.5)]; s_2 = s_2 / max(s_2);
s_cl = [0 ; cumsum(sum(diff(xy_cl,1,1).^2,2).^0.5)]; s_cl = s_cl / max(s_cl);

% Rotate and flip camberline back
xy_cl = [xy_cl(:,1) * cosd(gamma) - xy_cl(:,2) * sind(gamma) ...
    xy_cl(:,1) * sind(gamma) + xy_cl(:,2) * cosd(gamma)];
if flip_blade == 1
    xy_cl(:,2) = -xy_cl(:,2);
end

% Record surfaces and camberline
t.xy_1 = xy_1; t.xy_2 = xy_2; 
t.s_1 = s_1; t.s_2 = s_2;
t.i_1 = i_1; t.i_2 = i_2;
t.i_le = i_le; t.i_te = i_te; 
t.xy_cl = xy_cl; t.s_cl = s_cl;
t.chord = sum((xy_cl(1,:) - xy_cl(end,:)).^2)^0.5;
t.chi_le = 360*atan2(xy_cl(2,2)-xy_cl(1,2),xy_cl(2,1)-xy_cl(1,1))/(2*pi);
t.chi_te = 360*atan2(xy_cl(end,2)-xy_cl(end-1,2),xy_cl(end,1)-xy_cl(end-1,1))/(2*pi);

% Plot split surfaces
if plot_stuff == 1
    plot(xy_1(:,1),xy_1(:,2),'k-')
    plot(xy_2(:,1),xy_2(:,2),'r-')
    plot(xy_cl(:,1),xy_cl(:,2),'b-')
    plot(xy_orig(i_le,1),xy_orig(i_le,2),'b.')
    plot(xy_orig(i_te,1),xy_orig(i_te,2),'b.')
end

end