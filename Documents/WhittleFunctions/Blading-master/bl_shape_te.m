function xy = bl_shape_te(xy,m,plot_stuff)
% BL_SHAPE_TE  Create a biased and thinned trailing edge on a blade section

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Open figure window
if plot_stuff == 1
    figure(); hold on; grid on; box on; axis equal;
end

% Set options for equation solving
if plot_stuff == 1
    options = optimset('Display','iter'); 
else
    options = optimset('Display','off');
end

% Blend resolution
ni_blend = 100;

%% Characterise original trailing edge design

% Parametersise blade geometry
t = bl_parameterise_section(xy); chi_te = t.chi_te; i_circ = sort(t.i_circ);

% Rotate blade to trailing edge metal angle
xy = [xy(:,1) * cosd(-t.chi_te) - xy(:,2) * sind(-t.chi_te) ...
    xy(:,1) * sind(-t.chi_te) + xy(:,2) * cosd(-t.chi_te)];

% Section resolution
ni_old = size(xy,1);

% Plot original blade geometry
if plot_stuff == 1
    plot(xy(:,1),xy(:,2),'k-')
end

% Extract original trailing edge circle
xy_circ = xy(i_circ(1):i_circ(2),:);
ni_circ = size(xy_circ,1);

% Find original trailing edge circle radius and centre
p = CircleFitByTaubin(xy_circ);
r_te = p(3); xy_cen = p(1:2);

% Calculate arc angle of original trailing edge circle
th_te = sum(sum(diff(xy_circ,1,1).^2,2).^0.5) / r_te;

% Plot original trailing edge circle parameters
if plot_stuff == 1
    plot(xy_circ(:,1),xy_circ(:,2),'k.')
    plot(xy_cen(1),xy_cen(2),'k.');
end

% Interpolate blend points
xy_bl_1 = [interp1(t.s_cl_1_raw,xy(t.i_1,1),m.s,'pchip') ...
    interp1(t.s_cl_1_raw,xy(t.i_1,2),m.s,'pchip')];
xy_bl_2 = [interp1(t.s_cl_2_raw,xy(t.i_2,1),m.s,'pchip') ...
    interp1(t.s_cl_2_raw,xy(t.i_2,2),m.s,'pchip')];

% Calculate surface curvature and angle
[dydx_1,d2ydx2_1] = grad_mg(xy(t.i_1,1),xy(t.i_1,2));
[dydx_2,d2ydx2_2] = grad_mg(xy(t.i_2,1),xy(t.i_2,2));
dydx_bl_1 = interp1(t.s_cl_1_raw,dydx_1,m.s,'pchip');
dydx_bl_2 = interp1(t.s_cl_2_raw,dydx_2,m.s,'pchip');
d2ydx2_bl_1 = interp1(t.s_cl_1_raw,d2ydx2_1,m.s,'pchip');
d2ydx2_bl_2 = interp1(t.s_cl_2_raw,d2ydx2_2,m.s,'pchip');

% Plot blend points
if plot_stuff == 1
    plot(xy_bl_1(1),xy_bl_1(2),'ko')
    plot(xy_bl_2(1),xy_bl_2(2),'ko')
end


%% Generate new trailing edge design and insert into section

% Scale circle radius
r_te_new = r_te * m.fr;

% Find coordinates of new trailing edge circle
xy_cen_new = [xy_cen(1) + r_te - r_te_new ...
    xy_cen(2) + m.bs * r_te];

% Determine type of solution to use on first side
if isfield(m,'th') == 0 || isnan(m.th(1)) == 1

    % Solve non-linear equations for blend cubic fitting
    p = [1 1 1 1 pi/2];
    f = @(p) fun_te(p,xy_bl_1,dydx_bl_1,d2ydx2_bl_1,xy_cen_new,r_te_new);
    p_1 = fsolve(f,p,options); th_1 = p_1(5); p_1 = p_1(1:4);
    
else
    
    % Specify trailing edge circle arc directly
    th_1 = m.th(1);
    
    % Calculate boundary conditions at circle
    x_circ_1 = xy_cen_new(1) + r_te_new * cos(th_1);
    y_circ_1 = xy_cen_new(2) + r_te_new * sin(th_1);
    dydx_circ_1 = tan(th_1 - pi/2);
    
    % Blend with a quartic polynomial
    A = [xy_bl_1(1)^4  xy_bl_1(1)^3  xy_bl_1(1)^2  xy_bl_1(1)^1  1  ; ...
        4*xy_bl_1(1)^3  3*xy_bl_1(1)^2  2*xy_bl_1(1)  1  0 ; ...
        12*xy_bl_1(1)^2  6*xy_bl_1(1)  2  0  0  ; ...
        x_circ_1^4  x_circ_1^3  x_circ_1^2  x_circ_1  1  ; ...
        4*x_circ_1^3  3*x_circ_1^2  2*x_circ_1  1  0];
    b = [xy_bl_1(2) ; dydx_bl_1 ; d2ydx2_bl_1 ; y_circ_1 ; dydx_circ_1];
    p_1 = A \ b;
    
end

% Determine type of solution to use on second side
if isfield(m,'th') == 0 || isnan(m.th(2)) == 1

    % Solve non-linear equations for blend cubic fitting
    p = [1 1 1 1 -pi/2];
    f = @(p) fun_te(p,xy_bl_2,dydx_bl_2,d2ydx2_bl_2,xy_cen_new,r_te_new);
    p_2 = fsolve(f,p,options); th_2 = p_2(5); p_2 = p_2(1:4);
    
else
    
    % Specify trailing edge circle arc directly
    th_2 = m.th(2);
    
    % Calculate boundary conditions at circle
    x_circ_2 = xy_cen_new(1) + r_te_new * cos(th_2);
    y_circ_2 = xy_cen_new(2) + r_te_new * sin(th_2);
    dydx_circ_2 = tan(th_2 - pi/2);
    
    % Blend with a quartic polynomial
    A = [xy_bl_2(1)^4  xy_bl_2(1)^3  xy_bl_2(1)^2  xy_bl_2(1)^1  1  ; ...
        4*xy_bl_2(1)^3  3*xy_bl_2(1)^2  2*xy_bl_2(1)  1  0 ; ...
        12*xy_bl_2(1)^2  6*xy_bl_2(1)  2  0  0  ; ...
        x_circ_2^4  x_circ_2^3  x_circ_2^2  x_circ_2  1  ; ...
        4*x_circ_2^3  3*x_circ_2^2  2*x_circ_2  1  0];
    b = [xy_bl_2(2) ; dydx_bl_2 ; d2ydx2_bl_2 ; y_circ_2 ; dydx_circ_2];
    p_2 = A \ b;
    
end

% Generate new trailing edge circle
ni_circ_new = round(abs(th_1 - th_2) * ni_circ / th_te);
th = linspace(th_1,th_2,ni_circ_new)';
xy_circ_new = [xy_cen_new(1) + r_te_new * cos(th) xy_cen_new(2) + r_te_new * sin(th)];

% Plot new trailing edge circle
if plot_stuff == 1
    plot(xy_circ_new(:,1),xy_circ_new(:,2),'b.-')
    plot(xy_cen_new(1),xy_cen_new(2),'b.');
end

% Evaluate blend polynomials
x_1 = linspace(xy_bl_1(1),xy_circ_new(1,1),ni_blend)';
y_1 = polyval(p_1,x_1);
x_2 = linspace(xy_bl_2(1),xy_circ_new(end,1),ni_blend)';
y_2 = polyval(p_2,x_2);

% Plot new blend curves
if plot_stuff == 1
    plot(x_1,y_1,'g.-')
    plot(x_2,y_2,'g.-')
end

% Determine blend indices
i_bl_1 = find(xy(t.i_1,1) < x_1(1),1,'last'); i_bl_1 = t.i_1(i_bl_1);
i_bl_2 = find(xy(t.i_2,1) < x_2(1),1,'last'); i_bl_2 = t.i_2(i_bl_2);

% Assemble into the main section
if i_bl_2 < i_bl_1
    xy = [xy(1:i_bl_2,:) ; x_2 y_2 ; flip(xy_circ_new(2:end-1,:),1) ;...
        flip([x_1 y_1],1) ; xy(i_bl_1:end,:)];
else
    xy = [xy(1:i_bl_1,:) ; x_1 y_1 ; xy_circ_new(2:end-1,:) ; flip([x_2 y_2],1) ; xy(i_bl_2:end,:)];
end

% Smooth new points distribution
i_smooth = find(xy(:,1) > min([x_1(1) x_2(1)]) - 0.1*t.chord);
d = [0 ; cumsum(sum(diff(xy(i_smooth,:),1,1).^2,2).^0.5,1)]; d = d/d(end);
d_smooth = smooth(d,10);

% Ensure correct resolution of new section
ni_new = ni_old + ni_blend;
ni_smooth = length(i_smooth) - size(xy,1) + ni_new;
d_smooth = interp1(linspace(0,1,length(d_smooth))',d_smooth,linspace(0,1,ni_smooth)','pchip');

% Reinterpolate new section
xy_smooth = [interp1(d,xy(i_smooth,1),d_smooth,'pchip') ...
    interp1(d,xy(i_smooth,2),d_smooth,'pchip')];
xy = [xy(1:i_smooth(1)-1,:) ; xy_smooth ; xy(i_smooth(end)+1:end,:)];

% Plot new section
if plot_stuff == 1
    plot(xy(:,1),xy(:,2),'r-');
end

% Rotate back to correct stagger
xy = [xy(:,1) * cosd(chi_te) - xy(:,2) * sind(chi_te) ...
    xy(:,1) * sind(chi_te) + xy(:,2) * cosd(chi_te)];


end

function f = fun_te(p,xy_1,dydx_1,d2ydx2_1,xy_0,r)
% Evaluate residuals on trailing edge blend fitting

% Initialise residuals
f = zeros(5,1);

% Maintain value at blend point
f(1) = p(1) * xy_1(1)^3 + p(2) * xy_1(1)^2 + p(3) * xy_1(1) + p(4) - xy_1(2);

% Maintain tangency at blend point
f(2) = 3 * p(1) * xy_1(1)^2 + 2 * p(2) * xy_1(1) + p(3) - dydx_1;

% Maintain curvature at blend point
f(3) = 6 * p(1) * xy_1(1) + 2 * p(2) - d2ydx2_1;

% Calculate coordinates of new trailing edge circle
y_2 = xy_0(2) + r * sin(p(5));
x_2 = xy_0(1) + r * cos(p(5));

% Maintain value at trailing edge circle
f(4) = p(1) * x_2^3 + p(2) * x_2^2 + p(3) * x_2 + p(4) - y_2;

% Maintain tangency
f(5) = 3 * p(1) * x_2^2 + 2 * p(2) * x_2 + p(3) - tan(p(5) - pi/2);

end
