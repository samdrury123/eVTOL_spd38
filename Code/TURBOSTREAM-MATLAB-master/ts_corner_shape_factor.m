function h_corner = ts_corner_shape_factor(c,xr_hub,xr_cas,plot_stuff)
% TS_CORNER_SHAPE_FACTOR  Calculate corner separation shape factors
%
%   [h_corner] = TS_CORNER_SHAPE_FACTOR(c,xr_hub,xr_cas,plot_stuff)
%
%   c - structured cut of the blade o-mesh
%   xr_hub - matrix of hub line coordinates - [x1 r1 ; x2 r2 ; ...]
%   xr_cas - matrix of casing line coordinates - [x1 r1 ; x2 r2 ; ...]
%   plot_stuff - 0 or 1 for showing working

if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end

%% Prepare data

% Cut blade mesh to suction side only
r_hub = interp1(xr_hub(:,1),xr_hub(:,2),c.x(1,:),'pchip');
r_cas = interp1(xr_cas(:,1),xr_cas(:,2),c.x(1,:),'pchip');
r_nondim = (c.r(1,:) - r_hub) ./ (r_cas - r_hub);
[~,j_mid] = min(abs(r_nondim - 0.5));
t = bl_parameterise_section([c.x(:,j_mid) c.rt(:,j_mid)]);

s_min = 0.15; s_max = 0.95;
[~,i_min] = min(abs(t.s_cl_2_raw - s_min)); i_min = t.i_2(i_min);
[~,i_max] = min(abs(t.s_cl_2_raw - s_max)); i_max = t.i_2(i_max);

C = cell(2,1);
varnames = fieldnames(c);
for v = 1:length(varnames)
    if isfloat(c.(varnames{v})) == 1
        C{1}.(varnames{v}) = c.(varnames{v})(i_min:i_max,1:j_mid);
        C{2}.(varnames{v}) = c.(varnames{v})(i_min:i_max,j_mid:end);
    end
end

for n = 1:length(C)
    c = C{n};
    c = ts_secondary(c);
    varnames = fieldnames(c);
    for v = 1:length(varnames)
        if isnumeric(c.(varnames{v})) == 1
            c.(varnames{v}) = double(c.(varnames{v}));
        end
    end
    C{n} = c;
end

% Flip casing side blade mesh
c = C{2};
for v = 1:length(varnames)
    if isnumeric(c.(varnames{v})) == 1
        c.(varnames{v}) = flipdim(c.(varnames{v}),2);
    end
end
C{2} = c;

%% Hub and casing corner shape factors

% Open a figure window for working
if plot_stuff == 1
    h.corner = figure(); hold on; grid on; box on; axis equal;
    xlabel('Axial Distance'); ylabel('Radius');
end

% Pre-allocate corner data structure
h_corner = zeros(2,1);

% Set absolute distances based upon chord
d_rev = t.chord / 100; d_V = t.chord / 25;

% Loop over both hub and casing corner separations
for n = 1:length(C)
    c = C{n};

    % Check for significant trailing edge spearation at midheight
    if n == 1 && mean(mean(c.Vx(end-9:end,end-4:end))) < 0
        h_corner(1) = -1;
        continue
    elseif n == 2 && mean(mean(c.Vx(end-9:end,1:5))) < 0
        h_corner(2) = -1;
        continue
    end

    % Reshape coordinate and velocity arrays into vectors
    x = reshape(c.x,[],1); r = reshape(c.r,[],1);
    Vx = reshape(c.Vx,[],1); Vr = reshape(c.Vr,[],1); Vt = reshape(c.Vt,[],1);

    % Find points where flow vector is most aligned with spanwise
    % mesh vector
    W = cat(3,diff(c.x,1,2),diff(c.r,1,2),diff(c.rt,1,2));
    W = W ./ repmat(sum(W.^2,3).^0.5,[1 1 3]);
    V = cat(3,c.Vx,c.Vr,c.Vt); 
    V = V ./ repmat(sum(V.^2,3).^0.5,[1 1 3]);
    V_av = 0.5 * (V(:,1:end-1,:) + V(:,2:end,:));

    % Calculate alignment with mesh
    ph = acosd(sum(W.*V_av,3)); 
    
    % Select only points close to wall
    r_hub = interp1(xr_hub(:,1),xr_hub(:,2),c.x(:,1:end-1),'pchip');
    r_cas = interp1(xr_cas(:,1),xr_cas(:,2),c.x(:,1:end-1),'pchip');
    r_nondim = (c.r(:,1:end-1) - r_hub) ./ (r_cas - r_hub);
    if n == 1 
        ph(r_nondim > 0.05) = inf;
    else
        ph(r_nondim < 0.95) = inf;
    end

    % Select only points decelerating axially
    Vx_2 = V(:,2:end,1); Vx_1 = V(:,1:end-1,1);
    ph(Vx_2 > Vx_1) = inf;

    % Find most aligned point
    [ph_min,i_min] = min(ph,[],1); [ph_min,j_min] = min(ph_min); i_min = i_min(j_min);
    q = ph < ph_min + 5;

    % Fit straight line to points of vertical flow
    p_rev = polyfit(c.x(q),c.r(q),1);
    x_min = min(min(c.x(q))); x_max = max(max(c.x(q)));

    % Draw lines above and under flow reversal line
    n_rev = 100;
    x_1 = linspace(x_min,x_max,n_rev).';
    xr_1 = [x_1 polyval(p_rev,x_1) - d_rev];
    xr_2 = [x_1 polyval(p_rev,x_1) + d_rev];

    % Find grid points near above and under lines
    n_grid = length(x);
    d_1 = ( (repmat(x,[1 n_rev]) - repmat(xr_1(:,1).',[n_grid 1])).^2 + ...
        (repmat(r,[1 n_rev]) - repmat(xr_1(:,2).',[n_grid 1])).^2 ).^0.5;
    [~,q_1] = min(d_1,[],1);
    d_2 = ( (repmat(x,[1 n_rev]) - repmat(xr_2(:,1).',[n_grid 1])).^2 + ...
        (repmat(r,[1 n_rev]) - repmat(xr_2(:,2).',[n_grid 1])).^2 ).^0.5;
    [~,q_2] = min(d_2,[],1);            

    % Calculate average unit velocity vectors
    Vh_1 = [Vx(q_1) Vr(q_1) Vt(q_1)]; Vh_1 = Vh_1 ./ repmat(sum(Vh_1.^2,2).^0.5,[1 3]);
    Vh_2 = [Vx(q_2) Vr(q_2) Vt(q_2)]; Vh_2 = Vh_2 ./ repmat(sum(Vh_2.^2,2).^0.5,[1 3]);

    % Calculate maximum angle between vectors
    h_temp = sort(real(acos(sum(Vh_1 .* Vh_2,2))),1,'descend');
    h_corner(n) = mean(h_temp(1:4)); % Mean of worst 4 values
%     h_corner(n) = h_temp(1); % Worst value

    % Plot working for current iteration
    if plot_stuff == 1
        figure(h.corner);
        
        % Contour plot of radial gas angle
        surf(c.x,c.r,zeros(size(c.x)),double(c.Beta))

        % Point most aligned with the spanwise grid
        plot3(c.x(i_min,j_min),c.r(i_min,j_min),0.01,'c.','MarkerSize',20);

        % Line of points where flow reversal occurs
        plot(linspace(x_min,x_max,100),polyval(p_rev,linspace(x_min,x_max,100)),'k-','LineWidth',2)
        plot(x(q_1),r(q_1),'r.'); plot(x(q_2),r(q_2),'g.');

        % Velocity vectors above and below the flow reversal line
        plot([x(q_1) x(q_1)+Vh_1(:,1)*d_V].',[r(q_1) r(q_1)+Vh_1(:,2)*d_V].','r-')
        plot([x(q_2) x(q_2)+Vh_2(:,1)*d_V].',[r(q_2) r(q_2)+Vh_2(:,2)*d_V].','g-')
    end
end
