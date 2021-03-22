function [g,h] = ts_plot_surface(g,col,h,conv_cart,use_mesh,bids)
% TS_PLOT_SURFACE  Plot the mesh surfaces of a given TURBOSTREAM solution
%
%   [g,h] = TS_PLOT_SURFACE(g,col,h,conv_cart,use_mesh,bids)
%
%   g - cell array of TURBOSTREAM solution
%   C - vector of RGB colour intensity
%   h - figure windown handle
%   conv_cart - 0 or 1 to write cartesian coordinates and velocities
%   use_mesh - 0 or 1 to choose between mesh or surface plot
%   bids - vector of block ids to include or leave empty for all
%   g_surface - cell array of surface solutions

% Default to plotting polar coords
if exist('conv_cart','var') == 0
    conv_cart = 0;
end

% Default to using "mesh" function
if exist('use_mesh','var') == 0
    use_mesh = 1;
end

if exist('bids','var') == 0 || isempty(bids) == 1
    bids = ts_get_bids(g);
end

% Default to new figure
if exist('h','var') == 0 || isempty(h)
    h = figure();
elseif ishandle(h) == 0
    if mod(h,1) == 0
        figure(h);
    end
end

% Set figure properties
hold on; grid on; box on; axis equal; view(3);

% Get surface cell array
g = ts_get_surface(g,bids);

% Default to black, or repeat colour vector into colour matrices
if exist('col','var') == 0
    cols = cell(size(g));
    for bid = ts_get_bids(g)
        cols{bid+1} = zeros([size(g{bid+1}.x) 3]);
    end
elseif iscell(col) == 0
    cols = cell(size(g));
    for bid = ts_get_bids(g)
        cols{bid+1} = squeeze(repmat(reshape(col,1,1,1,3),size(g{bid+1}.x)));
    end
else
    cols = col;
end

% Plot all surfaces
for bid = ts_get_bids(g)
    
    % Extract coordinates
    x = squeeze(g{bid+1}.x); r = squeeze(g{bid+1}.r); rt = squeeze(g{bid+1}.rt);
    
    % Determine type of surface to plot
    if use_mesh == 1
        if conv_cart == 0
            h_surf = mesh(x,rt,r,cols{bid+1});
        else
            [y, z, x] = pol2cart(rt./r, r, x);
            h_surf = mesh(x,z,y,cols{bid+1});
        end
        set(h_surf,'LineWidth',0.5);
    elseif use_mesh == 2
        if conv_cart == 0
            h_surf = surf(x,rt,r,cols{bid+1});
        else
            [y,z,x] = pol2cart(rt./r, r, x);
            h_surf = surf(x,z,y,cols{bid+1});
        end
        if size(C,3) == 3
            set(h_surf,'EdgeColor',cols{bid+1}(1,1,:) * 0.5);
        end
    else
        if conv_cart == 0
            h1 = plot3(x,rt,r,'k-');
            h2 = plot3(x.',rt.',r.','k-');
            contourf(repmat(x,[1 1 2]),repmat(r,[1 1 2]),repmat(rt,[1 1 2]),...
                cat(3,zeros(size(x)),ones(size(x))),0.5);
        else
            [y,z,x] = pol2cart(rt./r, r, x);
            h1 = plot3(x,y,z,'k-');
            h2 = plot3(x.',y.',z.','k-');
            fill3(x,z,y,[0.9 0.9 0.9]);
        end
        set(h1,'LineWidth',0.5);
        set(h2,'LineWidth',0.5);
    end
end

end
