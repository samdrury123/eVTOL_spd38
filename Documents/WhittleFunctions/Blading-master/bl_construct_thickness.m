function t = bl_construct_thickness(c,ni,ote,plot_stuff)
% BL_CONSTRUCT_THICKNESS  Generate a distribution in shape space that matches physical parameters
%
%   t = BL_CONSTRUCT_THICKNESS(c,ni,ote,plot_stuff)
%
%   c - input struct containing thickness parameters
%   ni - number of points to evaluate thickness distribution on
%   ote - 0 or 1 to evaluate with an open trailing edge, -1 to smooth discontinuity
%   plot_stuff - 0 or 1 for showing working
%
%   t - output struct containing thickness distribution
%
%   c must contain fields:
%       rad_le - radius of curvature at leading edge
%       s_thick_max - chordwise position of max thickness
%       rad_thick_max - radius of curvature at max thickness location
%       wedge_te - trailing edge wedge angle
%       thick_te - trailing edge thickness
%
%   In the case of a closed trailing edge c must also contain:
%       thick_max - maximum thickness value
%       tchord - true chord measured along the camber line
%
%   t contains fields:
%       thick - non-dimensional thickness distribution
%       s_cl - non-dimensional distance along the camber line
%       S - evaluated shape space function for reference

% Default number of points
if exist('ni','var') == 0
    ni = 301;
end

% Default to closed trailing edge
if exist('ote','var') == 0
    ote = 0;
end

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Define weighting of points on different parts of the surface
if ote ~= 1
    ni_le = round(0.2 * ni); ni_si = round(0.6 * ni) + 2; ni_te = round(0.2 * ni); ni_sm = 9;
else
    ni_le = round(0.4 * ni); ni_si = round(0.6 * ni) + 1;
end

% Open figure window
if plot_stuff == 1
    figure();
    subplot(2,1,1); hold on; grid on; box on; xlabel('Chord'); ylabel('Thickness');
    subplot(2,1,2); hold on; grid on; box on; xlabel('Chord'); ylabel('Shape Function');
end

% Input parameters
rad_le = c.rad_le; s_thick_max = c.s_thick_max; rad_thick_max = c.rad_thick_max; 
wedge_te = c.wedge_te; thick_te = c.thick_te;

% Extra parameters for closed trailing edge
% if ote ~= 1
    thick_max = c.thick_max; tchord = c.tchord;
% end


%% Determine point spacing with a guess of the shape space function

% Guess a first thickness distribution
q_0 = [0 0 0 0 polyfit([0 1],[2 * (2 * rad_le)^0.5 ...
    tand(wedge_te) + thick_te],1)];
s_cl = hyperbolic_bunch(ni,ni/6e5,ni/3e5).';
S = polyval(q_0,s_cl);
thick = S .* (s_cl.^0.5 .* (1 - s_cl)) + s_cl * thick_te;

% Refine points distribution based on total distance
d = [0 ; cumsum(sum(diff([s_cl thick],1,1).^2,2).^0.5)];
d = (d - min(d)) / (max(d) - min(d));
s_cl = interp1(d,s_cl,s_cl,'pchip');


%% Construct shape space function from parameters

% Shape space leading edge point
S_1 = (2 * rad_le)^0.5;

% Shape space value at max thickness point
s_2 = s_thick_max;
S_2 = (1 - s_2 * thick_te) ./ (s_2.^0.5 .* (1 - s_2)); 

% Shape space gradient at max thickness point
a = -thick_te / (s_2^0.5 - s_2^1.5);
b = (1 - s_2 * thick_te) *  ((1/(2*s_2^0.5)) - 1.5*s_2^0.5) / ( (s_2^0.5 -s_2^1.5)^2 );
dSds_2 = a - b;

% Shape space curvature at max thickness point
d2tds2_2 = -1 / rad_thick_max;
a = 2 * ((1/(2*s_2^0.5)) - 1.5*s_2^0.5) * (-thick_te) / ( (s_2^0.5 -s_2^1.5)^2 );
b = 2 * (((1/(2*s_2^0.5)) - 1.5*s_2^0.5)^2) / ( (s_2^0.5 -s_2^1.5)^3 );
d = ( (-1 / (4*s_2^1.5)) - (0.75 / (s_2^0.5)) ) / ( (s_2^0.5 -s_2^1.5)^2 );
d2Sds2_2 = -a + (b-d) * (1 - s_2* thick_te) + d2tds2_2 / (s_2^0.5 -s_2^1.5);

% Shape space trailing edge point
S_3 = tand(wedge_te) + thick_te;

% Construct shape space spline from two cubics, rear section first
b = [S_2 ; dSds_2 ; d2Sds2_2 ; S_3];
A = [s_2^3 s_2^2 s_2 1 ; 3*s_2^2 2*s_2 1 0 ; 6*s_2 2 0 0 ; 1 1 1 1];
x_r = A\b;

% Calculate value, gradient and curvature at join
s_j = 0.3;
S_j = polyval(x_r,s_j); 
dSds_j = polyval(polyder(x_r),s_j); 
d2Sds2_j = polyval(polyder(polyder(x_r)),s_j);

% Construct cubic for front section
s_split = 0.11; s_stretch = 0.08;
b = [S_1 ; S_j ; dSds_j ; d2Sds2_j];
A = [(-s_stretch)^3 (-s_stretch)^2 -s_stretch 1 ; s_j^3 s_j^2 s_j 1 ; 3*s_j^2 2*s_j 1 0 ; 6*s_j 2 0 0];
x_f = A\b;

% Construct stretching cubic
A = [0 0 0 1 ; s_split^3 s_split^2 s_split 1 ; 3*s_split^2 2*s_split 1 0 ; 6*s_split 2 0 0];
b = [s_stretch ; 0 ; 0 ; 0];
x = A\b;

% Calculate thickness from shape space polynomials
[thick,S] = calc_thick(s_cl,x,x_f,x_r,s_j,s_split,thick_te);
    
% Plot shape space function
if plot_stuff == 1
    subplot(2,1,2); plot(s_cl,S,'b-'); v = axis;
end


%% Add a trailing edge circle

% Check if case is an open trailing edge
if ote ~= 1
    
    % Dimensionalise thickness distribution
    thick_dim = 0.5 * thick * thick_max; s_dim = s_cl * tchord;

    % Plot thickness distribution
    if plot_stuff == 1
        subplot(2,1,1); plot(s_cl,thick_dim,'b-')
    end
    
    % Display options and tight tolerancing
    if plot_stuff == 1
        options = optimoptions('fsolve','Display','iter','TolFun',1e-10); 
    else
        options = optimoptions('fsolve','Display','off','TolFun',1e-10);
    end

    % Initial guess of angle extent and circle centre position
    q = [s_dim(end) - 1.2 * thick_dim(end) 0.6 * pi];

    % Drive over camberwise length and angle at end of arc
    f = @(q) fun_te(q,s_dim,thick_dim,tchord,plot_stuff);
    q = fsolve(f,q,options);

    % Circle centre and radius
    xy_cen = [q(1) 0]; rad = s_dim(end) - q(1);
    
    % Evaluate circle geometry
    psi = linspace(q(2),0,ni_te)';
    xy_circ = repmat(xy_cen,[ni_te 1]) + rad * [cos(psi) sin(psi)];
    
    % Insert into thickness distribution
    q = s_dim < xy_circ(1,1);
    s_dim = [s_dim(q) ; xy_circ(:,1)]; thick_dim = [thick_dim(q) ; xy_circ(:,2)];
    
    % Non-dimensionalise thickness
    thick = 2 * thick_dim / thick_max; s_cl = s_dim / tchord;
else
    % Plot thickness distribution
    if plot_stuff == 1
        subplot(2,1,1); plot(s_cl,thick,'b-')
    end
end


%% Refine points spacing on the distribution

% Calculate total distance through the distribution
d = [0 ; cumsum(sum(diff([s_cl thick],1,1).^2,2).^0.5)]; d = d / d(end);

% Interpolate spacing at leading edge by distance
d_le = interp1(s_cl,d,0.05,'pchip');
s_le = interp1(d,s_cl,hyperbolic_bunch(ni_le,1.5/ni_le,0.5/ni_le)' * d_le,'pchip');
% s_le_2 = interp1(d,s_cl,linspace(0,d_le,ni_le)','pchip');

% % Interpolate spacing at leading edge by angle
% p = unit([diff(s_cl * tchord) diff(thick * thick_max * 0.5)]);
% dth = acosd(min(dot(p(1:end-1,:),p(2:end,:),2),1));
% dth(dth == 0) = 1e-5;
% th = [0 ; cumsum(dth)]; th_le = interp1(s_cl(1:end-1),th,0.05); 
% i = find(th > th_le,1);
% s_le_1 = interp1(th(1:i),s_cl(1:i),linspace(0,th_le,ni_le)','pchip');

% Take combination of leading edge spacings
% s_le = smooth(sort([s_le_1 ; s_le_2]),9); s_le = s_le(1:2:end);

% Spacing at trailing edge
if ote ~= 1
    ds_te = (xy_circ(2,1) - xy_circ(1,1)) / tchord;
    s_te = xy_circ(1,1) / tchord;
else
    ds_te = 0.5 / ni_si;
    s_te = 1;
end

% Spacing on surface
ds_le = s_le(end) - s_le(end-1); s_tot = (s_te - s_le(end));
s_si = hyperbolic_bunch(ni_si,ds_te/s_tot,ds_le/s_tot)' * s_tot + s_le(end);

% Assemble spacings together and recalculate thickness
s_cl = [s_le(1:end-1) ; s_si];
i = ni_le-13:ni_le+13; s_temp = smooth(s_cl(i),7,'sgolay');
s_cl(i(4:end-4)) = s_temp(4:end-4);
thick = calc_thick(s_cl,x,x_f,x_r,s_j,s_split,thick_te);

% Add in trailing edge circle if required
if ote ~= 1
    s_cl = [s_cl ; xy_circ(2:end,1) / tchord]; thick = [thick ; 2 * xy_circ(2:end,2) / thick_max];
end

% Smooth the thickness distribution at the trailing edge join
if ote == -1
    
    % Calculate boundary conditions on surface and trailing edge circle
    i_join = ni - ni_te + 1; i_1 = i_join - ni_sm; i_2 = i_join + ni_sm*2;
    [dydx,d2ydx2] = grad_mg(s_cl,thick);
    x_1 = s_cl(i_1); y_1 = thick(i_1); dydx_1 = dydx(i_1); d2ydx2_1 = d2ydx2(i_1);
    x_2 = s_cl(i_2); y_2 = thick(i_2); dydx_2 = dydx(i_2); d2ydx2_2 = d2ydx2(i_2);
    
    % Check the thickness distribution is convex before smoothing
%     d2ydx2_max = max(d2ydx2);
%     if d2ydx2_max > 0; disp(['Max curvature = ' num2str(d2ydx2_max)]); end;
    
    % Calculate quintic polynomial with matrix
    y = @(x) [x^5 x^4 x^3 x^2 x 1];
    dydx = @(x) [5*x^4 4*x^3 3*x^2 2*x 1 0];
    d2ydx2 = @(x) [20*x^3 12*x^2 6*x 2 0 0];
    A = [y(x_1) ; dydx(x_1) ; d2ydx2(x_1) ; y(x_2) ; dydx(x_2) ; d2ydx2(x_2)];
    b = [y_1 ; dydx_1 ; d2ydx2_1 ; y_2 ; dydx_2 ; d2ydx2_2];
    p = A \ b;
    
    % Evaluate and incorporate polynomial into thickness distribution
    thick(i_1+1:i_2-1) = polyval(p,s_cl(i_1+1:i_2-1));
    
end


%% Output results

% Calculate shape space function
S = (thick - s_cl * thick_te) ./ (s_cl.^0.5 .* (1 - s_cl));

% Plot final distributions
if plot_stuff == 1
    
    % Plot shape space function
    if plot_stuff == 1
        subplot(2,1,2); plot(s_cl,S,'k.-'); axis(v);
    end
    
    % Plot thickness distribution
    subplot(2,1,1); 
    if ote ~= 1
        plot(s_cl,0.5 * thick * thick_max,'k.-'); 
    else
        plot(s_cl,thick,'k.-');
    end
end

% Output variables
t.thick = thick; t.s_cl = s_cl; t.S = S;


end

function [thick,S] = calc_thick(s_cl,x,x_f,x_r,s_j,s_split,thick_te)
% Evaluate thickness from stretched polynomials

% Stretch chordwise spacing of camber line points
s_interp = s_cl;
s_interp(s_interp < s_split) = s_interp(s_interp < s_split) - polyval(x,s_interp(s_interp < s_split));

% Evaluate polynomials at stretched spacing
S = zeros(size(s_cl));
S(s_interp < s_j) = polyval(x_f,s_interp(s_interp < s_j));
S(s_interp >= s_j) = polyval(x_r,s_interp(s_interp >= s_j));

% Calculate real thickness
thick = S .* (s_cl.^0.5 .* (1 - s_cl)) + s_cl * thick_te;

end

function F = fun_te(q,s_dim,thick_dim,tchord,plot_stuff)
% Find errors in intersection and angle

% Intialise residuals
F = zeros(2,1);

% Circle centre and radius
xy_cen = [q(1) 0];
rad = s_dim(end) - q(1);

% End of circle location and angle
xy_end = [xy_cen(1) + rad * cos(q(2)) xy_cen(2) + rad * sin(q(2))];
psi_end = q(2) - pi/2;

% Angle on thickness distribution
psi_surf = atan(grad_mg(s_dim,thick_dim));

% Residuals of angle and y distance
F(1) = psi_end - interp1(s_dim,psi_surf,xy_end(1),'pchip');
F(2) = xy_end(2) - interp1(s_dim,thick_dim,xy_end(1),'pchip');

% Plot convergence of trailing edge
if plot_stuff == 1
    
    % Switch to thickness subplot
    subplot(2,1,1);
    
    % Calculate circle coordinates for plotting
    psi = linspace(0,psi_end+pi/2,100).';
    xy_circ = repmat(xy_cen,[100 1]) + rad * [cos(psi) sin(psi)];
    
    % Plot circle coordinates
    plot(xy_circ(:,1)/tchord,xy_circ(:,2),'r-')
    plot(xy_cen(1)/tchord,xy_cen(2),'r.')
end

end
