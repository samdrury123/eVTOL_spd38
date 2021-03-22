function [area_av_property,area] = ts_area_average(c,varname,dim)
% Function to calculate the area average of a property from a 2D cut

% dimension = 1 pitchwise averages
% dimension = 2 spanwise averages
% dimension = 3 whole area average
% To operate along a line take a cut 2 cells wide

% Generate property if not already existing
if isfield(c,varname) == 0
    c = ts_secondary(c);
end

% Work in double precision
property_names = fieldnames(c);
for i = 1:length(property_names)
    if strcmp(class(c.(property_names{i})), 'single') == 1
        c.(property_names{i}) = double(c.(property_names{i}));
    end
    
    % Work on 1D lines
    if size(c.(property_names{i}),2) == 1
        if strcmp(property_names{i},'rt') == 0
            c.(property_names{i}) = repmat(c.(property_names{i}),1,2);
        else
            c.rt = [zeros(size(c.rt,1),1) ones(size(c.rt,1),1)];
        end
    end
end

% Assign property of interest
varname = c.(varname);

% Calculate cell centred property
property_av = 0.25 * (varname(1:end-1,1:end-1,:) + varname(2:end,1:end-1,:) + ...
    varname(1:end-1,2:end,:) + varname(2:end,2:end,:));

% Calculate areas
c.t = c.rt ./ c.r;
t_av = 0.25 * (c.t(2:end,2:end) + c.t(1:end-1,1:end-1) + c.t(1:end-1,2:end) + c.t(2:end,1:end-1));
x1 = c.x(2:end,2:end) - c.x(1:end-1,1:end-1);
x2 = c.x(1:end-1,2:end) - c.x(2:end,1:end-1);
r1 = c.r(2:end,2:end) - c.r(1:end-1,1:end-1);
r2 = c.r(1:end-1,2:end) - c.r(2:end,1:end-1);
rt1 = c.r(2:end,2:end) .* (c.t(2:end,2:end) - t_av)...
    - c.r(1:end-1,1:end-1) .* (c.t(1:end-1,1:end-1) - t_av);
rt2 = c.r(1:end-1,2:end) .* (c.t(1:end-1,2:end) - t_av)...
    - c.r(2:end,1:end-1) .* (c.t(2:end,1:end-1) - t_av);

Ax = 0.5 * (r1.*rt2 - r2.*rt1);
Ar = 0.5 * (x2.*rt1 - x1.*rt2);
At = 0.5 * (x1.*r2 - x2.*r1);

% Calculate area of cell faces
mag_dA = (Ax.^2 + Ar.^2 + At.^2).^0.5;

% Choose dimensions to sum over
if dim ~= 3
    dim1 = dim;
    dim2 = dim;
else
    dim1 = 1;
    dim2 = 2;
end

% Integrate area average property and divide by total area
Top = sum(sum(property_av .* mag_dA, dim1),dim2);
Bottom = sum(sum(mag_dA, dim1),dim2);
area_av_property_av = Top ./ Bottom;

area = Bottom;

% Re-assign area averages to nodal coordinates
if dim1 == dim2
    area_av_property = zeros(1, numel(varname)/size(varname,dim));
    area_av_property(2:end-1) = 0.5 * (area_av_property_av(1:end-1) + ...
        area_av_property_av(2:end));

    area_av_property(1) = area_av_property_av(1);
    area_av_property(end) = area_av_property_av(end);
else
    area_av_property = area_av_property_av;
end


end