function [xy,xy_cam] = bl_construct_section(c,plot_stuff,ote,dev)
% BL_CONSTRUCT_SECTION  Create coordinates of a section from a definition of parameters
%
%   xy = BL_CONSTRUCT_SECTION(c,plot_stuff,ote,dev)
%   
%   c - data structure of section parameters
%   plot_stuff - 0 or 1 for showing working
%   ote - 0 or 1 for open trailing edge to run in MISES
%   dev - optional argument for thickness deviations
%   xy - section coordinate array
%   xy_cam - section camber line coordinates

% Set default to not plot
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Default to closed trailing edge
if exist('ote','var') == 0
    ote = 0;
end

% Set default to no measured deviations
if exist('dev','var') == 0
    dev = [];
end

% Set number of streamwise coordinates
ni = 301;


%% Evaluate camberline & thickness at desired resolution

% Calculate quartic thickness distribution in shape space
thick_param = bl_construct_thickness(c,ni,ote,0);
s_cl = thick_param.s_cl; thick = thick_param.thick * c.thick_max;

% Add on deviations to thickness distribution
if isempty(dev) == 0
    thick = thick + interp1(dev.s,dev.thick,s_cl,'spline') / thick_max;
end

% Evaluate camberline at same chordwise points as thickness distribution
cam_param = bl_construct_camber(c,s_cl); xy_cam = cam_param.xy_cam;


%% Add thickness distribution onto camberline

% Calculate camber line normals
chi = cam_param.chi; n_cam = [sind(chi) -cosd(chi)];

% Add thickness on to camberline to define section top and bottom lines
xy_1 = xy_cam + 0.5 * n_cam .* repmat(thick,[1 2]);
xy_2 = xy_cam - 0.5 * n_cam .* repmat(thick,[1 2]);

% Assemble both sides into one array
if ote == 0
    xy = [xy_1 ; flip(xy_2(1:end-1,:),1)];
else
    xy = [flip(xy_1,1) ; xy_2(2:end,:)];
end


%% Plot the blade section

% Check if plotting is turned on
if plot_stuff == 1
    
    % Open figure window
    figure; hold on; grid on; box on; axis equal;
    
    % Plot camber line and section surface
    plot(xy_cam(:,1),xy_cam(:,2),'b-')
    plot(xy(:,1),xy(:,2),'k.-')
    
end


end