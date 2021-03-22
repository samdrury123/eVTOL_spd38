function g = ts_set_rpms(g,rows,rpm,xr_hub,xr_cas,xr_stationary,xr_rotating)
% TS_SET_RPMS  Set rotating faces for different types of rows and shrouds
%
%   g = TS_SET_RPMS(g,rows,rpm,xr_hub,xr_cas,xr_stationary,xr_rotating)
%
%   g - cell of array for TURBOSTREAM mesh
%   rows - vector corresponding to row type
%   rpm - vector with rpm magnitude for each row
%   xr_hub - array or meridional hub coordinates
%   xr_cas - array or meridional casing coordinates
%   xr_stationary - optional array or stationary shroud faces
%   xr_rotating - optional array or rotating shroud faces
%
%   index to types of blade rows:
%       0 is shrouded stator
%       0.5 is unshrouded stator
%       1 is unshrouded rotor
%       1.5 is shrouded rotor

% List of rpm fields
rpm_ids = {'rpm' , 'rpmi1' , 'rpmi2' , 'rpmj1' , 'rpmj2' , 'rpmk1' , 'rpmk2'};

% Upscale number of points on curves
if exist('xr_stationary','var') ~= 0
    s = dist_2d(xr_rotating,1); s_interp = linspace(0,1,1000)';
    xr_rotating = [interp1(s,xr_rotating(:,1),s_interp,'linear') ...
        interp1(s,xr_rotating(:,2),s_interp,'linear')];
    s = dist_2d(xr_stationary,1);
    xr_stationary = [interp1(s,xr_stationary(:,1),s_interp,'linear') ...
        interp1(s,xr_stationary(:,2),s_interp,'linear')];
end

% Determine which bids fall in which rows
B = ts_get_row_bids(g,xr_hub,xr_cas);

% Default to whole domain if no blades can be found
if isempty(B) == 1
    B{1} = ts_get_bids(g);
end

% Set rpms in rows depending on configuration
for m = 1:length(B)
    if rows(m) == 0
        rpm_flag = [0 0 0 0 0 0 0];
    elseif rows(m) == 1
        rpm_flag = [1 1 1 1 0 1 1];
    elseif rows(m) == 1.5
        rpm_flag = [1 1 1 1 1 1 1];
    elseif rows(m) == 0.5
        rpm_flag = [0 0 0 1 0 0 0];
    end
    for bid = B{m}
        % Set default rpms depending on row geometry
        for n = 1:length(rpm_ids)
            bv_id = rpm_ids{n};
            g{bid+1}.bv.(bv_id) = rpm_flag(n);
        end
        
        % Set cavity well rpms depending on coordinates
        x = sum(sum(sum(g{bid+1}.x))) / numel(g{bid+1}.x);
        r = sum(sum(sum(g{bid+1}.r))) / numel(g{bid+1}.r);
        if r < interp1(xr_hub(:,1),xr_hub(:,2),x) && exist('xr_stationary','var') ~= 0 ...
                && exist('xr_rotating','var') ~= 0            
            xr = cat(4,g{bid+1}.x,g{bid+1}.r);
            
            g{bid+1}.bv.rpmi1 = shroud_face(squeeze(xr(1,:,1,:)),xr_stationary,xr_rotating);
            g{bid+1}.bv.rpmi2 = shroud_face(squeeze(xr(end,:,1,:)),xr_stationary,xr_rotating);
            g{bid+1}.bv.rpmj1 = shroud_face(squeeze(xr(:,1,1,:)),xr_stationary,xr_rotating);
            g{bid+1}.bv.rpmj2 = shroud_face(squeeze(xr(:,end,1,:)),xr_stationary,xr_rotating);

        end
        
        % Ensure fillet blocks remain stationary in cantilevered stator case
        if rows(m) == 0.5
            if r > 0.5 *(interp1(xr_hub(:,1),xr_hub(:,2),x) + interp1(xr_cas(:,1),xr_cas(:,2),x))
                for n = 1:length(rpm_ids)
                    bv_id = rpm_ids{n};
                    g{bid+1}.bv.(bv_id) = 0;
                end
            end
        end
        
        % Ensure fillet blocks remain rotating in cantilevered rotor case
        if rows(m) == 1
            if r < 0.5 *(interp1(xr_hub(:,1),xr_hub(:,2),x) + interp1(xr_cas(:,1),xr_cas(:,2),x)) && ...
                r > interp1(xr_hub(:,1),xr_hub(:,2),x)
                for n = 1:length(rpm_ids)
                    bv_id = rpm_ids{n};
                    g{bid+1}.bv.(bv_id) = 1;
                end
            end
        end
        
        % Set all blocks above casing line stationary
        if r > interp1(xr_cas(:,1),xr_cas(:,2),x)
            for n = 1:length(rpm_ids)
                bv_id = rpm_ids{n};
                g{bid+1}.bv.(bv_id) = 0;
            end
        end
    end
end

% Make sure all periodics are the same on both sides and are equal to the
% block speed
for bid = ts_get_bids(g)
    for pid = ts_get_pids(g,bid)
        if g{bid+1}.patch{pid+1}.attribute.kind == 5
            if g{bid+1}.patch{pid+1}.attribute.ist == 0 && g{bid+1}.patch{pid+1}.attribute.ien == 1
                g{bid+1}.bv.rpmi1 = g{bid+1}.bv.rpm;
            elseif g{bid+1}.patch{pid+1}.attribute.jst == 0 && g{bid+1}.patch{pid+1}.attribute.jen == 1
                g{bid+1}.bv.rpmj1 = g{bid+1}.bv.rpm;
            elseif g{bid+1}.patch{pid+1}.attribute.kst == 0 && g{bid+1}.patch{pid+1}.attribute.ken == 1
                g{bid+1}.bv.rpmk1 = g{bid+1}.bv.rpm;
            elseif g{bid+1}.patch{pid+1}.attribute.ist == g{bid+1}.attribute.ni-1 &&...
                    g{bid+1}.patch{pid+1}.attribute.ien == g{bid+1}.attribute.ni
                g{bid+1}.bv.rpmi2 = g{bid+1}.bv.rpm;
            elseif g{bid+1}.patch{pid+1}.attribute.jst == g{bid+1}.attribute.nj-1 &&...
                    g{bid+1}.patch{pid+1}.attribute.jen == g{bid+1}.attribute.nj
                g{bid+1}.bv.rpmj2 = g{bid+1}.bv.rpm;
            elseif g{bid+1}.patch{pid+1}.attribute.kst == g{bid+1}.attribute.nk-1 &&...
                    g{bid+1}.patch{pid+1}.attribute.ken == g{bid+1}.attribute.nk
                g{bid+1}.bv.rpmk2 = g{bid+1}.bv.rpm;
            end
        end
    end
end

% Connections of shrouds with main passage, ensure rotating on both sides
for m = 1:length(B)
    for bid = B{m}
        x = sum(sum(sum(g{bid+1}.x))) / numel(g{bid+1}.x);
        r = sum(sum(sum(g{bid+1}.r))) / numel(g{bid+1}.r);
        if r < interp1(xr_hub(:,1),xr_hub(:,2),x)
            for pid = ts_get_pids(g,bid)
                if g{bid+1}.patch{pid+1}.attribute.kind == 5
                    nxbid = g{bid+1}.patch{pid+1}.attribute.nxbid;
                    x = sum(sum(sum(g{nxbid+1}.x))) / numel(g{nxbid+1}.x);
                    r = sum(sum(sum(g{nxbid+1}.r))) / numel(g{nxbid+1}.r);
                    if r > interp1(xr_hub(:,1),xr_hub(:,2),x)
                        if rows(m) == 0
                            g{nxbid+1}.bv.rpmj1 = 1;
    %                         g{bid+1}.bv.rpmi1 = 1;
                            g{bid+1}.bv.rpmj2 = 1; % Using knife to knife model
                        elseif rows(m) == 1
                            g{nxbid+1}.bv.rpmj1 = 0;
                            g{bid+1}.bv.rpmj2 = 0;
                        end
                    end
                end
            end
        end
    end
end
    
% Inlet and outlet patches are stationary
for bid = ts_get_bids(g)
    for pid = ts_get_pids(g,bid)
        if g{bid+1}.patch{pid+1}.attribute.kind == 0
            g{bid+1}.bv.rpmi1 = 0;
        end
        if g{bid+1}.patch{pid+1}.attribute.kind == 1
            g{bid+1}.bv.rpmi2 = 0;
        end
    end
end

% Mixing planes are appropriate to row
for bid = ts_get_bids(g)
    for pid = ts_get_pids(g,bid)
        if g{bid+1}.patch{pid+1}.attribute.kind == 2
            if g{bid+1}.patch{pid+1}.attribute.ist == 0
                g{bid+1}.bv.rpmi1 = g{bid+1}.bv.rpm;
            else
                g{bid+1}.bv.rpmi2 = g{bid+1}.bv.rpm;
            end
        end
    end
end

% Set magnitudes of rotation in each row
for n = 1:length(B)
    for bid = B{n}
        for v = 1:length(rpm_ids)
            if g{bid+1}.bv.(rpm_ids{v}) ~= 0
                g{bid+1}.bv.(rpm_ids{v}) = rpm(n);
            end
        end
    end
end


end


function rpm = shroud_face(xr,xr_stationary,xr_rotating)
% Check distance to nearest face and return rotating or stationary

d_rotating = sort(min(sum((repmat(reshape(xr,[1 size(xr,1) 2]),[size(xr_rotating,1) 1 1]) - ...
    repmat(reshape(xr_rotating,[size(xr_rotating,1) 1 2]),[1 size(xr,1) 1]) ).^2,3).^0.5,[],1));

d_stationary = sort(min(sum((repmat(reshape(xr,[1 size(xr,1) 2]),[size(xr_stationary,1) 1 1]) - ...
    repmat(reshape(xr_stationary,[size(xr_stationary,1) 1 2]),[1 size(xr,1) 1]) ).^2,3).^0.5,[],1));

if mean(d_stationary(1:5)) < mean(d_rotating(1:5))
    rpm = 0;
else
    rpm = 1;
end        

end

