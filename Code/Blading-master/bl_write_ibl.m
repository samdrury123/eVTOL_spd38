function [] = bl_write_ibl(filename,b,flip_z,dt,nj)
% BL_WRITE_IBL  Create an ibl file of a blade definition to import into Creo/ProEngineer

% Default no blade flipping
if exist('flip_z','var') == 0
    flip_z = 0;
end

% Default no blade shifting
if exist('dt','var') == 0
    dt = 0;
end

% Default number of sections to write
if exist('nj','var') == 0
    nj = 22;
end

% Initialise file
fid = fopen(filename,'w');
fprintf(fid,'%s\n','Closed Index arclength');

% Get the correct blade definition
write_blade = 0;
if isfield(b,'xrrt_F0') == 1
    xrrt = b.xrrt_F0; write_blade = 1;
elseif isfield(b,'xrrt') == 1
    xrrt = b.xrrt; write_blade = 1;
end

% Check if required to write blade sections
if write_blade == 1

    % Reduce the number of sections
    if nj ~= size(xrrt,2)
        xrrt = xrrt(:,unique(round(linspace(1,size(xrrt,2),nj))),:);
    end
    
    % Shift the blade tangentially
    t = xrrt(:,:,3) ./ xrrt(:,:,2);
    t = t + dt;
    x = xrrt(:,:,1); r = xrrt(:,:,2);
    [y, z] = pol2cart(t, r);

    % Flip if necessary
    if flip_z == 1
        z = -z;
    end

    % Write blade coordinates
    xyz = cat(3,x,y,z); write_surf(fid,xyz);

    % Leading edge lines
    write_curve(fid,squeeze(xyz(1,:,:)));
    write_curve(fid,squeeze(xyz(end,:,:)));

    % Trailing edge line
    ite = round(size(xyz,1)/2);
    write_curve(fid,squeeze(xyz(ite,:,:)));
    
end

% Meridional lines if present
if isfield(b,'xr_hub'); write_curve(fid,b.xr_hub); end
if isfield(b,'xr_cas'); write_curve(fid,b.xr_cas); end

% Disk geometry if present
if isfield(b,'xr_disk') == 1
    
    % Loop over both sides of the disk
    for n = 1:2
        
        % Disk geometry as a curve
        write_curve(fid,b.xr_disk{n});
        
        % End and root geometry as straight lines
        write_lines(fid,b.xr_end{n})
        write_lines(fid,b.xr_root{n})
        
        % Check the disk has a spline definition
        if isfield(b,'xrrt_spl') == 1 && length(b.xrrt_spl) >= n && isempty(b.xrrt_spl{n}) == 0
        
            % Curves for spline surfaces
            xyz = bl_pol2cart(b.xrrt_spl{n});
            write_surf(fid,permute(xyz([1 end],:,:),[2 1 3]))

            % Perpendicular lines for spline surfaces
            write_curve(fid,squeeze(xyz(:,1,:)));
            i_mid = round(size(b.xrrt_spl{n},2)/2);
            write_curve(fid,squeeze(xyz(:,i_mid,:)));
        end        
    end
    
end

% Platform geometry if present
if isfield(b,'xrrt_plat') == 1
    xyz_plat = bl_pol2cart(b.xrrt_plat);
    write_curve(fid,squeeze(xyz_plat(:,1,:)));
    write_curve(fid,squeeze(xyz_plat(:,end,:)));
end

% Splitter geometry if present
if isfield(b,'xr_splitter') == 1
    for j = 1:length(b.xr_splitter)
        write_curve(fid,squeeze(b.xr_splitter{j}));
    end
end

% Close file
fclose(fid);


end


function [] = write_surf(fid,xyz)
% Write a surface to file    

% Loop over sections in the second dimension
for j = 1:size(xyz,2)
    
    % Print section coordinates
    fprintf(fid,'%s\n','begin section');
    fprintf(fid,'%s\n','begin curve');
    fprintf(fid,'%10.11f %10.11f %10.11f\n',squeeze(xyz(:,j,:))' * 1000);
        
end

end

function [] = write_lines(fid,xyz)
% Write segments of straight lines to file

% Increase dimension if meridional coordinates are specified
ni = size(xyz,1);
if size(xyz,2) == 2
    xyz = [xyz zeros(ni,1)];
end

% Print straight line segments to file
fprintf(fid,'begin section\nbegin curve\n%10.11f %10.11f %10.11f\n%10.11f %10.11f %10.11f\n',...
    1000 * xyz(sort([1:ni-1 2:ni]),:)');

end

function [] = write_curve(fid,xyz)
% Write a single curved line to file

% Increase dimension if meridional coordinates are specified
if size(xyz,2) == 2
    xyz = [xyz zeros(size(xyz,1),1)];
end

% Print curve to file
fprintf(fid,'%s\n','begin section');
fprintf(fid,'%s\n','begin curve');
fprintf(fid,'%10.11f %10.11f %10.11f\n',xyz' * 1000);

end