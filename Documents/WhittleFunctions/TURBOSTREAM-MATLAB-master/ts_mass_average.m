function [mass_av_property,mass_flux,area] = ts_mass_average(c,varname,dimension)
% TS_MASS_AVERAGE  Mass average a property from an unstructured or structured cut
%
%   [mass_av_property,mass_flux,area] = TS_MASS_AVERAGE(c,varname,dimension)
%
%   dimension = 1 :     average in spanwise direction
%   dimension = 2 :     average in pitchwise direction
%   dimension = 3 :     average over whole area

% Generate property if it doesn't already exist
if isfield(c,varname) == 0
    c = ts_secondary(c);
end

% Work in double precision
property_names = fieldnames(c);
for i = 1:length(property_names)
    if isa(c.(property_names{i}),'single') == 1
        c.(property_names{i}) = double(c.(property_names{i}));
    end
end

% Determine whether cut is structured or unstructured
if isfield(c,'faces') == 1
    
    % Unstructured cut
    
    % Get cartesian coordinates and face indices
    y = c.y; z = c.z; 
    faces = c.faces; f1 = faces(:,1); f2 = faces(:,2); f3 = faces(:,3);
    
    % Calculate the area of the face triangles using Heron's formula
    a1 = ((y(f2) - y(f1)).^2 + (z(f2) - z(f1)).^2).^0.5;
    a2 = ((y(f3) - y(f2)).^2 + (z(f3) - z(f2)).^2).^0.5;
    a3 = ((y(f1) - y(f3)).^2 + (z(f1) - z(f3)).^2).^0.5;
    d = sort([a1 a2 a3],2,'descend');
    a1 = d(:,1); a2 = d(:,2); a3 = d(:,3);
    dA = 0.25*sqrt(4*a1.*a1.*a2.*a2 - (a1.*a1 + a2.*a2 - a3.*a3).^2);

    % Correct imaginary areas
    if any(imag(dA)~=0)
        dA = real(dA);
    end

    % Calculate face mass flow rates assuming constant axial cut
    rovx = (c.rovx(f1) + c.rovx(f2) + c.rovx(f3))/3;
    dM = dA .* rovx; 
    
    % Calculate face average property values
    property_av = (c.(varname)(f1) + c.(varname)(f2) + c.(varname)(f3)) / 3;
    
    % Determine type of average to take, total or radial slice
    if dimension == 3
        
        % Calculate total mass, area and average of property
        mass_flux = sum(dM); area = sum(dA);
        mass_av_property = sum(property_av .* dM) / mass_flux;
        
    elseif dimension == 2
        
        % Number of bins to use
        nj = 31;
        
        % Calculate radial bins to divide solution into
        r_av = (c.r(f1) + c.r(f2) + c.r(f3)) / 3;
%         r_bins = linspace(min(r_av),max(r_av),nj);
        r_bins = interp1(linspace(0,1,length(r_av)),sort(r_av),linspace(0,1,nj));
        
        % Preallocate arrays for mass, area and property
        mass_flux = zeros(1,nj-1); area = zeros(1,nj-1); 
        mass_av_property = zeros(1,nj-1);
        
        % Loop over all bins and calculate variables
        for j = 1:length(r_bins)-1
            q = r_av > r_bins(j) & r_av < r_bins(j+1);
            mass_flux(j) = sum(dM(q)); area(j) = sum(dA(q));
            mass_av_property(j) = sum(property_av(q) .* dM(q)) ./ mass_flux(j);
        end
    end
    
else
    
    % Structured cut

    % Work on 1D lines
    if size(c.x,2) == 1
        for i = 1:length(property_names)
            if isa(c.(property_names{i}),'double') == 1
                c.(property_names{i}) = repmat(c.(property_names{i}),1,2);
            end
        end
        c.rt = [zeros(size(c.rt,1),1) c.r(:,1) * 2 * pi];
    end

    % Calculate cell centred averages and relative velocity
    rovx_av = 0.25 * (c.rovx(1:end-1,1:end-1) + c.rovx(2:end,1:end-1) + ...
        c.rovx(1:end-1,2:end) + c.rovx(2:end,2:end));
    rovr_av = 0.25 * (c.rovr(1:end-1,1:end-1) + c.rovr(2:end,1:end-1) + ...
        c.rovr(1:end-1,2:end) + c.rovr(2:end,2:end));
    rpm = c.bv.rpm;
    omega = rpm * 2 * pi() / 60;
    c.rovt = c.rorvt ./ c.r - (c.r .* c.ro * omega);
    rovt_av = 0.25 * (c.rovt(1:end-1,1:end-1) + c.rovt(2:end,1:end-1) + ...
        c.rovt(1:end-1,2:end) + c.rovt(2:end,2:end));
    property_av = 0.25 * (c.(varname)(1:end-1,1:end-1,:) + c.(varname)(2:end,1:end-1,:) + ...
        c.(varname)(1:end-1,2:end,:) + c.(varname)(2:end,2:end,:));

    % Calculate areas by cross product of diagonals
    c.t = c.rt ./ c.r;
    t_av = 0.25 * (c.t(2:end,2:end) + c.t(1:end-1,1:end-1) + ...
        c.t(1:end-1,2:end) + c.t(2:end,1:end-1));
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

    % Find mass flux
    dM = Ax.*rovx_av + Ar.*rovr_av + At.*rovt_av;

    % Calculate area of cell faces
    mag_dA = (Ax.^2 + Ar.^2 + At.^2).^0.5;

    % Choose dimensions to sum over
    if dimension ~= 3
        dim1 = dimension;
        dim2 = dimension;
    else
        dim1 = 1;
        dim2 = 2;
    end

    % Integrate mass average property and divide by total mass flux
    Top = sum(sum(property_av .* dM, dim1),dim2);
    Bottom = sum(sum(dM, dim1),dim2);
    mass_av_property_av = Top ./ Bottom;

    area = sum(sum(mag_dA,dim1),dim2);

    % Re-assign mass averages to nodal coordinates
    if dim1 == dim2
        mass_av_property = zeros(1, numel(c.(varname))/size(c.(varname),dimension));
        mass_av_property(2:end-1) = 0.5 * (mass_av_property_av(1:end-1) + ...
            mass_av_property_av(2:end));

        mass_av_property(1) = mass_av_property_av(1);
        mass_av_property(end) = mass_av_property_av(end);
    else
        mass_av_property = mass_av_property_av;
    end

    mass_flux = Bottom.';

end


end

