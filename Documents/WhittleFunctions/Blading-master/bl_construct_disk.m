function b = bl_construct_disk(b,plot_stuff)
% BL_CONSTRUCT_DISK  Construct disk geometry from specified parameters

% Default to plot output
if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end

% Find the axial location of the disk centre
if isfield(b,'x_cen') == 0
    x_cen = 0.5 * (b.x_root(1) + b.x_root(2));
else
    x_cen = b.x_cen;
end

% Fillet radius
rad_fil = b.rad_fil;

% Tolerance on evaluating polynomials and splines
r_tol = 0.05e-3; th_tol = 1;

% Repeat end angle if only one is specified
if numel(b.th_end) == 1
    b.th_end = ones(1,2) * b.th_end;
end

% Default splines on both sides of the disc
if isfield(b,'spl_flag') == 0
    b.spl_flag = [1 1];
end


%% Plotting preparation
if plot_stuff == 1

    % Open figure window and add axes labels
    h = figure(); hold on; grid on; box on; axis equal;
    xlabel('Axial Direction') ; ylabel('Radial Direction');

    % Plot blade geometry
    plot(b.xrrt(:,:,1),b.xrrt(:,:,2),'k.-')

    % Plot meridional lines
    plot(b.xr_hub(:,1),b.xr_hub(:,2),'k-');
    plot(b.xr_cas(:,1),b.xr_cas(:,2),'k-');
end


%% End of disk geometry including dovetail for concentrictity and tapered overhang for stiffness

% Check if dovetail is specified
if isfield(b,'d_dove') == 1

    % Define chamfer size
    d_cham = 0.15e-3;

    % Define dovetail fixings at front of the disk
    d = b.d_dove; x = b.xr_hub(1,1); r = b.xr_hub(1,2);
    b.xr_end{1} = [x r ; x r-d ; x+d r-d ; x+d r-2*d ; x r-2*d ; x r-3*d];

    % Define dovetail fixings at rear of the disk
    x = b.xr_hub(end,1); r = b.xr_hub(end,2);
    b.xr_end{2} = [x r ; x r-d ; x+d-d_cham r-d ; x+d r-d-d_cham ; x+d r-2*d+d_cham ; ...
        x+d-d_cham r-2*d ; x r-2*d ; x r-3*d];
    
else
    
    % Define plain end coordinates on front
    x = b.xr_hub(1,1); r = b.xr_hub(1,2); b.xr_end{1} = [x r ; x r-b.d_end];
    
    % Define plain end coordinates on rear
    x = b.xr_hub(end,1); r = b.xr_hub(end,2); b.xr_end{2} = [x r ; x r-b.d_end];
end

% Define overhang to join with centre of disk at front
t = b.t_end; dx = x_cen - t/2 - b.xr_end{1}(end,1);
b.xr_end{1} = [b.xr_end{1} ; b.xr_end{1}(end,1) + dx b.xr_end{1}(end,2) - dx*sind(b.th_end(1))];

% Define overhang to join with centre of disk at rear
dx = b.xr_end{2}(end,1) - x_cen - t/2;
b.xr_end{2} = [b.xr_end{2} ; b.xr_end{2}(end,1) - dx b.xr_end{2}(end,2) - dx*sind(b.th_end(2))];
    

%% Centre of disk geometry including balancing rings

% Define boundary conditions at front of disk
x_root = x_cen - b.t_root/2; r_root = b.r_root;
x_end = b.xr_end{1}(end,1); r_end = b.xr_end{1}(end,2); 

% Calculate parabola at front of disk
r = linspace(r_end,r_root,round((r_end - r_root)/r_tol))';
A = [2*r_end 1 0 ; r_end^2 r_end 1 ; r_root^2 r_root 1]; x = [0 ; x_end ; x_root];
p{1} = A \ x; b.xr_disk{1} = [polyval(p{1},r) r];

% Define boundary conditions at rear of disk
x_root = x_cen + b.t_root/2;
x_end = b.xr_end{2}(end,1); r_end = b.xr_end{2}(end,2); 

% Calculate parabola at rear of disk
r = linspace(r_end,r_root,round((r_end - r_root)/r_tol))';
A = [2*r_end 1 0 ; r_end^2 r_end 1 ; r_root^2 r_root 1]; x = [0 ; x_end ; x_root];
p{2} = A \ x; b.xr_disk{2} = [polyval(p{2},r) r];

% Define balance rings
if isfield(b,'d_bal') == 1
    
    % Loop over both sides of the disk
    dx = [-b.d_bal b.d_bal]; r_bal = b.r_bal; d = b.d_bal; rad = b.rad_bal;
    for n = 1:2

        % Find axial location of face at balance radius
        x_bal = polyval(p{n},r_bal);

        % Draw extent of balance ring
        xr_bal = [x_bal r_bal+d/2 ; x_bal+dx(n) r_bal+d/2 ; x_bal+dx(n) r_bal-d/2 ; x_bal r_bal-d/2];

        % Calculate radiuses on edge of the balance ring
        xr_cen = [xr_bal(2,1) - sign(dx(n)) * rad(1) xr_bal(2,2) - rad(1)];
        xr_fil_1 = fillet_2d(xr_cen,rad(1),xr_bal([1 2],:),xr_bal([2 3],:),0.5);
        xr_cen = [xr_bal(3,1) - sign(dx(n)) * rad(1) xr_bal(3,2) + rad(1)];
        xr_fil_2 = fillet_2d(xr_cen,rad(1),xr_bal([3 4],:),xr_bal([2 3],:),0.5);    
        xr_bal = [xr_bal(1,:) ; xr_fil_1 ; flip(xr_fil_2,1) ; xr_bal(4,:)];

        % Calculate fillet on upper surface
        xr_cen = [x_bal + sign(dx(n)) * rad(2) r_bal + d/2 + rad(2)];
        xr_fil = fillet_2d(xr_cen,rad(2),b.xr_disk{n},xr_bal(1:2,:),0.5);
        xr_bal = [xr_fil ; xr_bal(2:end,:)];

        % Calculate fillet on lower surface
        xr_cen = [x_bal + sign(dx(n)) * rad(2) r_bal - d/2 - rad(2)];
        xr_fil = fillet_2d(xr_cen,rad(2),flip(b.xr_disk{n},1),xr_bal(end:-1:end-1,:),0.5);
        xr_bal = [xr_bal(1:end-1,:) ; flip(xr_fil,1)];

        % Add balance geometry onto parabola
        b.xr_disk{n} = [b.xr_disk{n}(b.xr_disk{n}(:,2) > xr_bal(1,2),:) ; xr_bal ; ...
            b.xr_disk{n}(b.xr_disk{n}(:,2) < xr_bal(end,2),:)];
            plot(xr_bal(:,1),xr_bal(:,2))
    end
end

%% Root of disc geometry including threaded hole and splines

% Check if spline is specified
b.xrrt_spl = cell(1,2);
if isfield(b,'rad_spl') == 1

    % Get spline parameters
    r_root = b.r_root; rad = b.rad_spl; r_hol = b.r_hol; N = b.N_spl; alpha = b.alpha_spl;

    % Calculate the length available for straight sections between the arcs
    l_root = 2 * pi * r_root - 4 * rad * N * sind(alpha);
    l_spl = 2 * pi * r_hol - 4 * rad * N * sind(alpha);

    % Check if the splines fit with current settings
    if l_root < 0; error(['l = ' num2str(l_root) ' Splines do not fit at root']); end;
    if l_spl < 0; error(['l = ' num2str(l_spl) ' Splines do not fit at spline radius']); end;

    % Define radial and tangential coordinates
    nj = round((r_root - r_hol) / r_tol); nk = 2 * alpha / th_tol;
    r = repmat(linspace(r_hol,r_root,nj)',[1 nk]); t = repmat(linspace(0,2*pi/N,nk),[nj 1]);
    rt = r.*t;

    % Initialise axial coordinates
    x = zeros(nj,nk);

    % Calculate extent of chordwise chunks and straight line chunks
    l_arc = 2 * rad * sind(alpha);
    l_lin = pi * r / N - l_arc;

    % Axial displacements of straight line chunks
    dx = l_lin * tand(alpha);

    % Define half of first arc
    rt_cen = 0; q = rt < l_arc/2;
    beta = asind((rt_cen - rt(q)) / rad);
    x(q) = rad * (cosd(beta) - cosd(alpha)) + dx(q) / 2;

    % Define second half of first arc
    rt_cen = repmat(rt(:,end),[1 nk]); q = rt > rt_cen - l_arc/2;
    beta = asind((rt_cen(q) - rt(q)) / rad);
    x(q) = rad * (cosd(beta) - cosd(alpha)) + dx(q) / 2;

    % Define second arc
    rt_cen = repmat(rt(:,end) / 2,[1 nk]); 
    q = rt < rt_cen + l_arc/2 & rt > rt_cen - l_arc/2;
    beta = asind((rt_cen(q) - rt(q)) / rad);
    x(q) = - rad * (cosd(beta) - cosd(alpha)) - dx(q) / 2;

    % Define first straight line segment
    rt_cen = repmat(rt(:,end) / 2,[1 nk]) - l_lin/2 - l_arc/2; 
    q = rt <= rt_cen + l_lin/2 & rt >= rt_cen - l_lin/2;
    x(q) = - (rt(q) - rt_cen(q)) .* dx(q) ./ l_lin(q);

    % Define second straight line segment
    rt_cen = repmat(rt(:,end),[1 nk]) - l_lin/2 - l_arc/2; 
    q = rt <= rt_cen + l_lin/2 & rt >= rt_cen - l_lin/2;
    x(q) = (rt(q) - rt_cen(q)) .* dx(q) ./ l_lin(q);

    % Repeat segments to make up the whole spline face
    dt = 0:2*pi/N:2*pi; dt = repmat(reshape(repmat(dt(1:end-1)',[1 nk-1])',[1 (nk-1)*N]),[nj 1]);
    t = rt(:,1:end-1) ./ r(:,1:end-1); t = repmat(t,[1 N]); t = t + dt;
    x = repmat(x(:,1:end-1),[1 N]); r = repmat(r(:,1:end-1),[1 N]); rt = r.*t;

    % Repeat final row
    x = [x x(:,1)]; r = [r r(:,1)]; rt = [rt 2*pi*r(:,1)];

    % Print spline depth
    if plot_stuff == 1
        disp(['Spline Depth = ' num2str(max(x(:)) - min(x(:)))])
    end

    % Offset segments to inlet and outlet coordinates
    if b.spl_flag(1) == 1; b.xrrt_spl{1} = cat(3,-x + b.x_root(1),r,rt); end;
    if b.spl_flag(2) == 1; b.xrrt_spl{2} = cat(3,x + b.x_root(2),r,rt); end;

    % Plot spline faces
    if plot_stuff == 1
        figure(h); 
        for n = 1:2
            if b.spl_flag(n) == 1
                xyz = bl_pol2cart(b.xrrt_spl{n});
                mesh(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),zeros(size(xyz)));
            end
        end
    end

    % Meridional coordinates of root
    if b.spl_flag(1) == 1
        b.xr_root{1} = [b.xr_disk{1}(end,:) ; min(b.xrrt_spl{1}(end,:,1)) r_root ; ...
            min(b.xrrt_spl{1}(1,:,1)) r_hol];
    else
        b.xr_root{1} = [b.xr_disk{1}(end,:) ; b.x_root(1) b.xr_disk{1}(end,2) ; b.x_root(1) r_hol];
    end
    if b.spl_flag(2) == 1
        b.xr_root{2} = [b.xr_disk{2}(end,:) ; max(b.xrrt_spl{2}(end,:,1)) r_root ; ...
            max(b.xrrt_spl{2}(1,:,1)) r_hol];
    else
        b.xr_root{2} = [b.xr_disk{2}(end,:) ; b.x_root(2) b.xr_disk{2}(end,2) ; b.x_root(2) r_hol];
    end
else
    
    % Plain root coordinates on front
    b.xr_root{1} = [b.xr_disk{1}(end,:) ; b.xr_hub(1,1) b.r_root ; b.xr_hub(1,1) b.r_hol];
    
    % Plain root coordinates on rear
    b.xr_root{2} = [b.xr_disk{2}(end,:) ; b.xr_hub(end,1) b.r_root ; b.xr_hub(end,1) b.r_hol]; 
end


%% Bolt holes through the splines

% Check if bolt holes are specified
if isfield(b,'r_bolt') == 1
    
    % Loop over all bolts
    b.xyz_hol = cell(b.N_bolt,1);
    for n = 1:b.N_bolt
        
        % Polar coordinates
        t = 2 * pi * (n-1) / (b.N_bolt); r = b.r_bolt;
        if b.spl_flag(1) == 1
            x_1 = b.xr_root{2}(2,1); x_2 = b.xr_root{1}(2,1);
        else
            x_1 = b.xr_root{1}(2,1); x_2 = b.xr_root{2}(2,1);
        end
        xrrt = [x_1 r t * r ; x_2 r t * r];
        
        % Cartesian coordinates
        b.xyz_hol{n} = bl_pol2cart(xrrt);
        
        % Plot hole
        if plot_stuff == 1; figure(h); plot3(b.xyz_hol{n}(:,1),b.xyz_hol{n}(:,2),b.xyz_hol{n}(:,3)); end;
        
    end
    
end


%% Add fillets onto joins between sections

% Define sign conventions for both sides of the disk
dx = [1 -1];

% Loop over both sides of the disk
for n = 1:2
    
    % Calculate fillet between end and disk
    xr_cen = [b.xr_end{n}(end,1) - rad_fil * dx(n) b.xr_end{n}(end,2) - rad_fil];
    xr_fil = fillet_2d(xr_cen,rad_fil,b.xr_end{n}(end:-1:end-1,:),flip(b.xr_disk{n},1),0.5);
    b.xr_disk{n} = [xr_fil ; b.xr_disk{n}(b.xr_disk{n}(:,2) < xr_fil(end,2),:)];
    b.xr_end{n}(end,:) = xr_fil(1,:);
    
    % Calculate fillet between disk and root
    xr_cen = [b.xr_root{n}(1,1) - 0.7*rad_fil * dx(n) b.xr_root{n}(1,2) + 0.7*rad_fil];
    xr_fil = fillet_2d(xr_cen,rad_fil,b.xr_disk{n},b.xr_root{n}(1:2,:),0.5);
    b.xr_disk{n} = [b.xr_disk{n}(b.xr_disk{n}(:,2) > xr_fil(1,2),:) ; xr_fil];
    b.xr_root{n}(1,:) = xr_fil(end,:);
    
end

% Join straight lines for dovetail disk or add external radii for simple disk
if isfield(b,'rad_out') == 0

    % Move straight lines onto centre of disk geometry
    for n = 1:2
        b.xr_disk{n} = [b.xr_end{n}(end-1,:) ; b.xr_disk{n} ; b.xr_root{n}(2,:)];
        b.xr_end{n}(end,:) = []; b.xr_root{n}(1,:) = [];
    end

else

    % Add fillets
    for n = 1:2
    
        % Calculate fillet on disk end
        xr_cen = [b.xr_end{n}(2,1) + b.rad_out * dx(n) b.xr_end{n}(2,2) + b.rad_out];
        xr_fil = fillet_2d(xr_cen,b.rad_out,flip(b.xr_end{n}(1:2,:),1),flip(b.xr_end{n}(2:end,:),1),0.5);
        b.xr_end{n} = [b.xr_end{n}(1,:) ; xr_fil ; b.xr_end{n}(end,:)];
        
        % Calculate fillet on root end
        xr_cen = [b.xr_root{n}(2,1) + b.rad_out * dx(n) b.xr_root{n}(2,2) - b.rad_out];
        xr_fil = fillet_2d(xr_cen,b.rad_out,flip(b.xr_root{n}(1:2,:),1),...
            flip(b.xr_root{n}(2:end,:),1),0.5);
        b.xr_root{n} = [b.xr_root{n}(1,:) ; xr_fil ; b.xr_root{n}(end,:)];        
    end
    
    % Join meridional disk coordinates together
    tol = 0.1e-3; varnames = {'xr_end' 'xr_root'};
    for a = 1:2
        for v = 1:length(varnames)
            s = dist_2d(b.(varnames{v}){a}); 
            s_interp = unique([linspace(0,s(2),round(s(2) / tol))' ; s ; ...
                linspace(s(end-1),s(end),round((s(end) - s(end-1)) / tol))']);
            b.(varnames{v}){a} = interp1(s,b.(varnames{v}){a},s_interp);
        end
        b.xr_disk{a} = [b.xr_end{a} ; b.xr_disk{a}(2:end,:) ; b.xr_root{a}(2:end,:)];
        b.xr_end{a} = zeros(0,2); b.xr_root{a} = zeros(0,2); 
    end
end


%% Plot final meridional lines

% Check whether to plot
if plot_stuff == 1
    
    % Change to correct figure window
    figure(h);
    
    % Loop over both sides and plot segments
    for n = 1:2
        plot(b.xr_root{n}(:,1),b.xr_root{n}(:,2),'.-')
        plot(b.xr_disk{n}(:,1),b.xr_disk{n}(:,2),'.-')
        plot(b.xr_end{n}(:,1),b.xr_end{n}(:,2),'.-')
    end
end


end
