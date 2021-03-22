function xrrt = bl_construct_fillet(xrrt,F,xr_hub,xr_cas,plot_stuff)
% BL_CONSTRUCT_FILLET  Add variable radius fillets onto a blade definition
%
%   xrrt = BL_CONSTRUCT_FILLET(xrrt,F,xr_hub,xr_cas,plot_stuff)
%
%   xrrt - 3D float array of blade polar coordinates
%   F - fillet data structure defining geometry
%   xr_hub - 2D float array of hub meridional coordinates
%   xr_cas - 2D float array of casing meridional coordinates
%   plot_stuff - 0 or 1 for showing working
%
%   F is a data structure
%       F.flag is flag for hub or casing fillets present 1x2 vector, 0 or 1
%       F.hub or F.cas for hub or casing definitionas
%       F.hub.s1 is vector of chordwise coordinates on side 1
%       F.hub.s2 is vector of chordwise coordinates on side 2
%       F.hub.rad1 is matching vector of fillet radii on side 1
%       F.hub.rad2 is matching vector of fillet radii on side 2
%       F.tol is float of blade stretch for good meridional intersection

% Default number of sections to use for fillet definition and buffer layer
nj_buff = 6; 
nj_fillet = 35;

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end


%% Prepare blade orientation

% Blade size
[ni,nj,~] = size(xrrt);

% Check nose points down, if not flip in theta direction
[~,i_le] = min(xrrt(:,1,1)); [~,i_te] = max(xrrt(:,1,1));
if xrrt(i_le,1,3) > xrrt(i_te,1,3)
    flip_blade = 1;
    xrrt(:,:,3) = - xrrt(:,:,3);
else
    flip_blade = 0;
end

% Ensure points are specified clockwise by crossing adjacent edges together
a = squeeze(diff(xrrt(:,1,:),1));
b = cross(a(2:end,:),a(1:end-1,:));
if sum(sign(b(:,2))) < 0
    xrrt = flip(xrrt,1);
end


%% Extract regions of blade geometry in the vicinity of the fillets

% Limit hub and casing lines to vicinity of blade
x_min = min(min(xrrt(:,:,1))); x_max = max(max(xrrt(:,:,1)));

% Initialise fillet data for hub and casing side
d = [-1 ; 1]; xr_mer = {xr_hub ; xr_cas}; fname = {'hub' ; 'cas'}; 

% Find midspan meridional line
x = linspace(x_min,x_max,100)';
xr_mid = [x 0.5 * (interp1(xr_hub(:,1),xr_hub(:,2),x) + interp1(xr_cas(:,1),xr_cas(:,2),x))];

% Initialise coordinate arrays
xrrt_end = cell(2,1);

% Loop over hub and casing fillets
for n = find(F.flag == 1)
    
    % Get maximum fillet radius
    rad_max = max(max(F.(fname{n}).rad1),max(F.(fname{n}).rad2));
    
    % Find indices of end parts of blade
    if n == 1
        r_end = max(xr_mer{n}(xr_mer{n}(:,1) > x_min & xr_mer{n}(:,1) < x_max,2)) + rad_max * 3;
        j_end = 1:find(max(xrrt(:,:,2),[],1) >= r_end - F.tol,1,'first');
    elseif n == 2
        r_end = min(xr_mer{n}(xr_mer{n}(:,1) > x_min & xr_mer{n}(:,1) < x_max,2)) - rad_max * 3;
        j_end = find(min(xrrt(:,:,2),[],1) <= r_end + F.tol,1,'last'):nj;
    end

    % Reinterpolate end parts to increased resolution
    xrrt_end{n} = zeros(ni,nj_fillet,3);
    for i = 1:ni
        r_hub_new = linspace(xrrt(i,j_end(1),2),xrrt(i,j_end(end),2),nj_fillet);
        xrrt_end{n}(i,:,1) = interp1(xrrt(i,:,2),xrrt(i,:,1),r_hub_new,'pchip');
        xrrt_end{n}(i,:,2) = interp1(xrrt(i,:,2),xrrt(i,:,2),r_hub_new,'pchip');
        xrrt_end{n}(i,:,3) = interp1(xrrt(i,:,2),xrrt(i,:,3),r_hub_new,'pchip');
    end
end

% Start plotting
if plot_stuff == 1
    
    % Open figure window and set limits
    figure(); hold on; axis equal; view(3); xyz = bl_pol2cart(xrrt);
    axis([min(min(xyz(:,:,1)))-rad_max max(max(xyz(:,:,1)))+rad_max ...
        min(min(xyz(:,:,3)))-rad_max max(max(xyz(:,:,3)))+rad_max...
        min(min(xyz(:,:,2)))-rad_max max(max(xyz(:,:,2)))+rad_max]);
    set(gca,'clippingstyle','rectangle');
%     for n = find(F.flag == 1)
%         mesh(xrrt_end{n}(:,:,1),xrrt_end{n}(:,:,3),xrrt_end{n}(:,:,2))
%     end

    % Plot the meridional surfaces
    for n = find(F.flag == 1)
        rt_max = max(max(xrrt(:,:,3))); rt_min = min(min(xrrt(:,:,3))); drt = rt_max - rt_min;
        xrrt_surf = repmat(reshape(xr_mer{n},[],1,2),[1 30 1]);
        xrrt_surf = cat(3,xrrt_surf,repmat(linspace(rt_min-0.2*drt,rt_max+0.2*drt,30),...
            [size(xr_mer{n},1) 1 1]));
        xyz_surf = bl_pol2cart(xrrt_surf);
        mesh(xyz_surf(:,:,1),xyz_surf(:,:,3),xyz_surf(:,:,2),zeros(size(xyz_surf,1),30,3))
    end
end


%% Calculate fillet geometries

% Preallocate fillet arrays
xrrt_fillet = cell(2,1); xyz_end = cell(2,1);

% Loop over hub and casing fillets
for n = find(F.flag == 1)

    % Calculate blade surface normals
    xyz_end{n} = bl_pol2cart(xrrt_end{n});
    n_xyz = bl_normals(xyz_end{n});

    % Calculate distribution of fillet radii along the blade
    if max(abs(diff(F.(fname{n}).rad1))) ~= 0 || max(abs(diff(F.(fname{n}).rad2))) ~= 0
        rad = zeros(ni,1);
        t = bl_parameterise_section([xrrt_end{n}(:,1,1),xrrt_end{n}(:,1,3)],0,0,0);
        rad(t.i_1) = interp1(F.(fname{n}).s1,F.(fname{n}).rad1,t.s_cl_1_raw,'pchip');
        rad(t.i_2) = interp1(F.(fname{n}).s2,F.(fname{n}).rad2,t.s_cl_2_raw,'pchip');
        rad(1) = rad(2); rad(end) = rad(end-1);
    else
        rad = ones(ni,1) * F.(fname{n}).rad1(1);
    end
    
    % Roll circles along meridional line and spheres along blade to find centre of the fillet circle
    xrrt_cen = zeros(ni,3); xrrt_end_i = zeros(ni,3); xrrt_mer_i = zeros(ni,3);
    for i = 1:ni
        
        % Roll a circle along meridional line and record centre coords
        n_mer = [d(n)*diff(xr_mer{n}(:,2)) -d(n)*diff(xr_mer{n}(:,1))]; 
        n_mer = 0.5 * (n_mer(2:end,:) + n_mer(1:end-1,:));
        n_mer = n_mer ./ repmat(sum(n_mer.^2,2).^0.5,[1 2]);
        xr_mer_temp = xr_mer{n}(2:end-1,:); xr_mer_cen = xr_mer_temp + rad(i) * n_mer;
        
        % Roll a sphere along line of blade sections and get linear intersection
        xyz_bl_cen = squeeze(xyz_end{n}(i,:,:) + rad(i) * n_xyz(i,:,:));
        xrrt_bl_cen = bl_cart2pol(xyz_bl_cen);
        [~, r, ii, jj] = intersections(xrrt_bl_cen(:,1),xrrt_bl_cen(:,2),...
            xr_mer_cen(:,1),xr_mer_cen(:,2));
        [xrrt_cen(i,2), c] = max(r); ii = round(ii(c)); jj = round(jj(c));

        % Record valid nearby indices for interpolation
        ii_min = max(ii-2,1); ii_max = min(ii+2,nj_fillet);
        
        % Increase resolution in the vicinity of the intersection and pinpoint circle centre
        r_1 = linspace(xrrt_bl_cen(ii_min,2),xrrt_bl_cen(ii_max,2),100);
        x_1 = interp1(xrrt_bl_cen(ii_min:ii_max,2),xrrt_bl_cen(ii_min:ii_max,1),r_1,'pchip');
        x_2 = linspace(xr_mer_cen(jj-2,1),xr_mer_cen(jj+2,1),100);
        r_2 = interp1(xr_mer_cen(jj-2:jj+2,1),xr_mer_cen(jj-2:jj+2,2),x_2,'pchip');
        [xrrt_cen(i,1), xrrt_cen(i,2)] = intersections(x_1,r_1,x_2,r_2);
        xrrt_cen(i,3) = interp1(xrrt_bl_cen(ii_min:ii_max,2),xrrt_bl_cen(ii_min:ii_max,3),...
            xrrt_cen(i,2),'pchip');

        % Find coordinates of intersections with the blade
        for m = 1:3
            xrrt_end_i(i,m) = interp1(xrrt_bl_cen(ii_min:ii_max,2),xrrt_end{n}(i,ii_min:ii_max,m),...
                xrrt_cen(i,2),'pchip');
        end

        % Find coordinates of intersection with meridional line
        for m = 1:2
            xrrt_mer_i(i,m) = interp1(xr_mer_cen(jj-2:jj+2,1),xr_mer_temp(jj-2:jj+2,m),...
                xrrt_cen(i,1),'pchip');
        end
        xrrt_mer_i(i,3) = (xrrt_cen(i,3) ./ xrrt_cen(i,2)) * xrrt_mer_i(i,2);
    end
    
    % Calculate normals from centre coordinates to blade intersection coordinates
    n_bl = bl_pol2cart(xrrt_end_i) - bl_pol2cart(xrrt_cen);
    n_bl = n_bl ./ repmat(sum(n_bl.^2,2).^0.5,[1 3]);
    
    % Calculate normals from centre coordinates to meridional intersection coordinates
    n_mer = bl_pol2cart(xrrt_mer_i) - bl_pol2cart(xrrt_cen);
    n_mer = n_mer ./ repmat(sum(n_mer.^2,2).^0.5,[1 3]);    
    
    % Calculate angle extent of fillet
    q = dot(n_bl,n_mer,2); psi_mer = acosd(q); 
    
    % Correct meridional normal by blade normal and angle and renormalise
    n_mer = n_mer - repmat(q,[1 3]) .* n_bl; 
    n_mer = n_mer ./ repmat(sum(n_mer.^2,2).^0.5,[1 3]);

    % Correct meridional curve intersection radius
    xrrt_mer_i(:,2) = interp1(xr_mer{n}(:,1),xr_mer{n}(:,2),xrrt_mer_i(:,1),'pchip');

    % Plot output from intersection search
    if plot_stuff == 1
        xyz = {bl_pol2cart(xrrt_cen) ; bl_pol2cart(xrrt_end_i) ; bl_pol2cart(xrrt_mer_i)};
        cols = [0 0 0 ; 1 0 0 ; 1 0 0];
        for m = 1:3
            plot3(xyz{m}(:,1),xyz{m}(:,3),xyz{m}(:,2),'.-','color',cols(m,:))
        end
    end

    % Calculate spanwise extent and therefore spacing of fillet sections
    r_nondim_blade = (xrrt_end{n}(:,:,2) - ...
        interp1(xr_mer{n}(:,1),xr_mer{n}(:,2),xrrt_end{n}(:,:,1),'pchip')) ./ ...
        (interp1(xr_mid(:,1),xr_mid(:,2),xrrt_end{n}(:,:,1),'pchip') -...
        interp1(xr_mer{n}(:,1),xr_mer{n}(:,2),xrrt_end{n}(:,:,1),'pchip'));
    r_nondim_end = (xrrt_end_i(:,2) - ...
        interp1(xr_mer{n}(:,1),xr_mer{n}(:,2),xrrt_end_i(:,1),'pchip')) ./ ...
        (interp1(xr_mid(:,1),xr_mid(:,2),xrrt_end_i(:,1),'pchip') -...
        interp1(xr_mer{n}(:,1),xr_mer{n}(:,2),xrrt_end_i(:,1),'pchip'));
    [~,i_max] = max(abs(r_nondim_end));
    
    % Define radiussed geometry up to maximum angle
    psi = linspace(psi_mer(i_max),0,nj_fillet);
    xyz_cen = bl_pol2cart(xrrt_cen);
    xrrt_rad = bl_cart2pol(repmat(reshape(xyz_cen(i_max,:),[1 1 3]),[1 nj_fillet 1]) + rad(i_max) ...
        * (repmat(reshape(n_bl(i_max,:),[1 1 3]),[1 nj_fillet 1]) .* repmat(cosd(psi),[1 1 3]) +...
        repmat(reshape(n_mer(i_max,:),[1 1 3]),[1 nj_fillet 1]) .* repmat(sind(psi),[1 1 3])));
    r_nondim_rad = (xrrt_rad(:,:,2) - ...
        interp1(xr_mer{n}(:,1),xr_mer{n}(:,2),xrrt_rad(:,:,1),'pchip')) ./ ...
        (interp1(xr_mid(:,1),xr_mid(:,2),xrrt_rad(:,:,1),'pchip') -...
        interp1(xr_mer{n}(:,1),xr_mer{n}(:,2),xrrt_rad(:,:,1),'pchip'));

    % Plot maximum radius line
    if plot_stuff == 1
        xyz = bl_pol2cart(xrrt_rad);
        plot3(xyz(:,:,1),xyz(:,:,3),xyz(:,:,2),'b.-')
    end

    % Loop over all i locations and calculate fillet geometry
    xrrt_fillet{n} = zeros(ni,nj_fillet,3);
    for i = 1:ni
        
        % Define variation in fillet angle
        psi = linspace(psi_mer(i),0,nj_fillet*5);
        xrrt_temp = bl_cart2pol(repmat(xyz_cen(i,:),[nj_fillet*5 1]) + ...
            rad(i) * (repmat(n_bl(i,:),[nj_fillet*5 1]) .* repmat(cosd(psi'),[1 3]) +...
            repmat(n_mer(i,:),[nj_fillet*5 1]) .* repmat(sind(psi'),[1 3])));
        r_nondim = (xrrt_temp(:,2) - ...
            interp1(xr_mer{n}(:,1),xr_mer{n}(:,2),xrrt_temp(:,1),'pchip'))./ ...
            (interp1(xr_mid(:,1),xr_mid(:,2),xrrt_temp(:,1),'pchip') -...
            interp1(xr_mer{n}(:,1),xr_mer{n}(:,2),xrrt_temp(:,1),'pchip'));

        % Division between radiused fillet and rest of the blade
        j_div = find(max(r_nondim) <= r_nondim_rad,1,'first')-1;
        psi = interp1(r_nondim,psi,r_nondim_rad(1:j_div),'pchip'); psi(1) = psi_mer(i);

        % Fillet geometry
        xrrt_fillet{n}(i,1:j_div,:) = bl_cart2pol(repmat(reshape(xyz_cen(i,:),[1 1 3]),[1 j_div 1]) +...
            rad(i) * (repmat(reshape(n_bl(i,:),[1 1 3]),[1 j_div 1]) .* repmat(cosd(psi),[1 1 3]) +...
            repmat(reshape(n_mer(i,:),[1 1 3]),[1 j_div 1]) .* repmat(sind(psi),[1 1 3])));

        % Rest of geometry above fillet
        xrrt_fillet{n}(i,j_div+1:nj_fillet,1) = interp1(r_nondim_blade(i,:),xrrt_end{n}(i,:,1),...
            r_nondim_rad(j_div+1:nj_fillet),'pchip');
        xrrt_fillet{n}(i,j_div+1:nj_fillet,2) = interp1(r_nondim_blade(i,:),xrrt_end{n}(i,:,2),...
            r_nondim_rad(j_div+1:nj_fillet),'pchip');
        xrrt_fillet{n}(i,j_div+1:nj_fillet,3) = interp1(r_nondim_blade(i,:),xrrt_end{n}(i,:,3),...
            r_nondim_rad(j_div+1:nj_fillet),'pchip');

    end
end


%% Assemble blade and fillet sections into one array with smooth buffer sections

% Loop over hub and casing fillets
for n = find(F.flag == 1)
    
    % Specify maximum index where buffer geometry meets blade
    [~,j] = min(abs(xrrt(i,:,2) - xrrt_fillet{n}(i,end,2))); dr_2 = diff(xrrt(i,j:j+1,2));
    [~,j] = min(abs(xrrt(i,:,2) - (xrrt_fillet{n}(i,end,2) - d(n) * dr_2 * 0.5 * nj_buff)));
    
    % Buffer geometry
    xrrt_buff = zeros(ni,nj_buff,3);
    for i = 1:ni
        
        % Get spacings and coordinates of fillet and blade
        dr_1 = diff(xrrt_fillet{n}(i,end-1:end,2)); r_1 = xrrt_fillet{n}(i,end,2);
        dr_2 = -d(n) * diff(xrrt(i,j:j+1,2)); r_2 = xrrt(i,j,2);

        % Specify polynomial distribution to smooth spacing
        r = poly_bunch(r_1,r_2,dr_1,dr_2,nj_buff);
       
        % Interpolate buffer section
        xrrt_buff(i,:,:) = cat(3,interp1(xrrt(i,:,2),xrrt(i,:,1),r,'pchip'),r,...
            interp1(xrrt(i,:,2),xrrt(i,:,3),r,'pchip'));
    end    

    % Assemble fillet and buffer into blade
    if n == 1
        xrrt = [xrrt_fillet{n} xrrt_buff(:,2:end-1,:) xrrt(:,j:end,:)];
    elseif n == 2
        xrrt = [xrrt(:,1:j,:) xrrt_buff(:,end-1:-1:2,:) flip(xrrt_fillet{n},2)];
    end
end

% Record new blade size
nj = size(xrrt,2);


%% Uncrease surface and smooth spacings on fillet

% Loop over all sections
for j = 1:nj
    
    % Find damaged surface by relative angle
    a = squeeze(xrrt(2:end,j,:) - xrrt(1:end-1,j,:));
    a = a ./ repmat(sum(a.^2,2).^0.5,[1 3]);
    phi = sum(a(1:end-1,:) .* a(2:end,:),2);
        
    % Find centre of damage
    i = find(phi < 0.98);
    if isempty(i) == 0
        i = round(mean(i)); nd = 5;
        i = i-nd:i+nd; %i(i < 1 | i > size(xrrt_blade,1)) = [];
        i(i < 1) = i(i < 1) + size(xrrt,1);
        i(i > size(xrrt,1)) = i(i > size(xrrt,1)) - size(xrrt,1);
    end

    % Check if any damaged parts are found
    if isempty(i) == 0
        
        % Pull out affected coordinates
        xrrt_temp = squeeze(xrrt(i,j,:));

        % Resolved into local coordinates
        s = xrrt_temp(end,:) - xrrt_temp(1,:); s = s ./ sum(s.^2).^0.5;
        xrrt_temp = xrrt_temp - repmat(xrrt_temp(1,:),[length(i) 1]);
        q = sum(xrrt_temp .* repmat(s,[length(i) 1]),2);
        q = q / max(q);
        x = xrrt(i,j,1);
        r = xrrt(i,j,2);
        rt = xrrt(i,j,3);
        
        % Recompute parabolic fit through damaged points
        xrrt(i,j,1) = polyval(polyfit(q,x,2),linspace(0,1,length(i)));
        xrrt(i,j,2) = polyval(polyfit(q,r,2),linspace(0,1,length(i)));
        xrrt(i,j,3) = polyval(polyfit(q,rt,2),linspace(0,1,length(i)));
    end
end

% Sew leading edge
xrrt_temp = 0.5 * sum(xrrt([1 end],:,:),1);
xrrt(1,:,:) = xrrt_temp; xrrt(end,:,:) = xrrt_temp;

% Loop over both fillets and smooth spacing on hub and casing
for n = find(F.flag == 1)
    
    % Determine section indices to work on
    if n == 1
        j_fillet = 1:(nj_fillet+nj_buff);
    elseif n == 2
        j_fillet = nj:-1:nj - (nj_fillet+nj_buff) + 1;
    end

    % Calculate spacing on fillet sections
    s = [zeros(1,length(j_fillet)) ; cumsum(sum(diff(xrrt(:,j_fillet,:),1,1).^2,3).^0.5)]; 
    s = s ./ repmat(s(end,:),[ni 1]);
    
    % Smooth spacings
    s_smooth = zeros(size(s)); 
    for j = 1:size(s_smooth,2); s_smooth(:,j) = smooth(s(:,j),round(ni/20)); end;
    for i = 1:size(s_smooth,1); s_smooth(i,:) = smooth(s_smooth(i,:),5); end;
    
    % Vary quantity of smoothing on each section, reducing where blending into the blade
    frac = linspace(0,1,nj_fillet+nj_buff);
    s_interp = repmat(frac,[ni 1]) .* s + repmat(1 - frac,[ni 1]) .* s_smooth;
    
    % Reinterpolate sections
    for j = 1:length(j_fillet)
        xrrt(:,j_fillet(j),:) = cat(3,interp1(s(:,j),xrrt(:,j_fillet(j),1),s_interp(:,j),'pchip'),...
            interp1(s(:,j),xrrt(:,j_fillet(j),2),s_interp(:,j),'pchip'),...
            interp1(s(:,j),xrrt(:,j_fillet(j),3),s_interp(:,j),'pchip'));
    end
end


%% Tidy up geometry and return new surface

% Correct hub and casing lines stretch blade to lie on hub and casing
for n = find(F.flag == 1)
    if n == 1
        xrrt(:,1,2) = interp1(xr_mer{n}(:,1),xr_mer{n}(:,2),xrrt(:,1,1),'pchip');
    elseif n == 2
        xrrt(:,end,2) = interp1(xr_mer{n}(:,1),xr_mer{n}(:,2),xrrt(:,end,1),'pchip');
    end
end

% Stretch blade very slightly to ensure improved hub and casing intersection
n_sec = size(xrrt,2);
dr = [linspace(-F.tol,0,round(n_sec/2)) linspace(0,F.tol,n_sec-round(n_sec/2))];
xrrt(:,:,2) = xrrt(:,:,2) + repmat(dr,[ni 1 1]);

% Plot blade mesh
if plot_stuff == 1
    xyz = bl_pol2cart(xrrt);
    mesh(xyz(:,:,1),xyz(:,:,3),xyz(:,:,2))
end

% Flip blade back if necessary
if flip_blade == 1
    xrrt(:,:,3) = - xrrt(:,:,3);
end

end