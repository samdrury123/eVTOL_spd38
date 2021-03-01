function b = bl_construct_mer(b,d_tol)
% BL_CONSTRUCT_MER  Produce hub and casing lines based upon a set of design parameters

% Specify hub and casing line resolution
if exist('d_tol','var') == 0
    d_tol = 2e-4;
end

% Decide which type of curves to construct
if isstruct(b) == 1
    
    % Construct curves within a single row
    b = construct_cubic(b,d_tol);
    
elseif iscell(b) == 1
    
    % Construct curves between multiple rows
    xr_hub = []; xr_cas = []; B = b; clear b;
    
    % Loop over all blade rows
    for n = 1:length(B)-1

        % Use meridional lines for each blade
        xr_hub = [xr_hub ; B{n}.xr_hub]; xr_cas = [xr_cas ; B{n}.xr_cas];

        % Get boundary conditions of neighbouring blades
        b.A_le = B{n}.A_te; b.r_le = B{n}.r_te; 
        b.x_hub_le = B{n}.xr_hub(end,1); b.x_cas_le = B{n}.xr_cas(end,1); 
        b.drdx_hub_le = B{n}.drdx_hub_te; b.drdx_cas_le = B{n}.drdx_cas_te; 

        b.A_te = B{n+1}.A_le; b.r_te = B{n+1}.r_le; 
        b.x_hub_te = B{n+1}.xr_hub(1,1); b.x_cas_te = B{n+1}.xr_cas(1,1); 
        b.drdx_hub_te = B{n+1}.drdx_hub_le; b.drdx_cas_te = B{n+1}.drdx_cas_le;

        % Generate intermediate cubic
        b = construct_cubic(b,d_tol);

        % Use meridional lines between each blade
        xr_hub = [xr_hub ; b.xr_hub]; xr_cas = [xr_cas ; b.xr_cas];
    end

    % Use final blade row curves
    xr_hub = [xr_hub ; B{end}.xr_hub]; xr_cas = [xr_cas ; B{end}.xr_cas];

    % Add on inlet quadratic curves
    x_1 = B{1}.x_in; drdx_1 = 0;
    
    drdx_hub_2 = (xr_hub(2,2) - xr_hub(1,2)) / (xr_hub(2,1) - xr_hub(1,1));
    xr_hub = [construct_quadratic(x_1,xr_hub(1,1),xr_hub(1,2),drdx_1,drdx_hub_2,d_tol) ; xr_hub];
    
    drdx_cas_2 = (xr_cas(2,2) - xr_cas(1,2)) / (xr_cas(2,1) - xr_cas(1,1));
    xr_cas = [construct_quadratic(x_1,xr_cas(1,1),xr_cas(1,2),drdx_1,drdx_cas_2,d_tol) ; xr_cas];   
    
    % Add on exit quadratic casing
    x_1 = B{end}.x_out; drdx_1 = 0;
    drdx_cas_2 = (xr_cas(end,2) - xr_cas(end-1,2)) / (xr_cas(end,1) - xr_cas(end-1,1));
    xr_cas = [xr_cas ; construct_quadratic(x_1,xr_cas(end,1),xr_cas(end,2),drdx_1,drdx_cas_2,d_tol)];
    [~,i] = unique(xr_cas(:,1)); xr_cas = xr_cas(i,:);
            
    % Add on quadratic hub or quartic nozzle
    if isfield(B{end},'A_noz') == 0
        
        % No nozzle case
        drdx_2 = (xr_hub(end,2) - xr_hub(end-1,2)) / (xr_hub(end,1) - xr_hub(end-1,1));
        xr_hub = [xr_hub ; construct_quadratic(x_1,xr_hub(end,1),xr_hub(end,2),drdx_1,drdx_2,d_tol)];
        
    else
        
        % Nozzle with some hub quadratic for a longer transition
        if isfield(B{end},'x_const') == 1
            drdx_2 = (xr_hub(end,2) - xr_hub(end-1,2)) / (xr_hub(end,1) - xr_hub(end-1,1));
            xr_hub = [xr_hub ; construct_quadratic(B{end}.x_const,xr_hub(end,1),xr_hub(end,2),...
                drdx_1,drdx_2,d_tol)];
        end
        
        % Nozzle boundary conditions
        x_2 = xr_hub(end,1); r_2 = xr_hub(end,2); x_3 = B{end}.x_noz; drdx_3 = 0; 
        drdx_2 = (xr_hub(end,2) - xr_hub(end-1,2)) / (xr_hub(end,1) - xr_hub(end-1,1));
        r_3 = (interp1(xr_cas(:,1),xr_cas(:,2),x_3,'pchip')^2 - B{end}.A_noz / pi)^0.5;
        
        % Generate polynomial coefficients
        A = [4*x_1^3 3*x_1^2 2*x_1 1 0 ; ...
            x_2^4 x_2^3 x_2^2 x_2 1 ; ...
            4*x_2^3 3*x_2^2 2*x_2 1 0 ; ...
            x_3^4 x_3^3 x_3^2 x_3 1 ; ...
            4*x_3^3 3*x_3^2 2*x_3 1 0];
        p = A \ [drdx_1 ; r_2 ; drdx_2 ; r_3 ; drdx_3];
        
        % Polynomial coefficients matching radius of curvature
%         d2rdx2_2 = (xr_hub(end,2) - 2 * xr_hub(end-1,2) + xr_hub(end-2,2)) / ...
%             (xr_hub(end,1) - xr_hub(end-1,1))^2;
%         A = [20*x_2^3 12*x_2^2 6*x_2 2 0 0 ; ...
%             5*x_1^4 4*x_1^3 3*x_1^2 2*x_1 1 0 ; ...
%             x_2^5 x_2^4 x_2^3 x_2^2 x_2 1 ; ...
%             5*x_2^4 4*x_2^3 3*x_2^2 2*x_2 1 0 ; ...
%             x_3^5 x_3^4 x_3^3 x_3^2 x_3 1 ; ...
%             5*x_3^4 4*x_3^3 3*x_3^2 2*x_3 1 0];
%         p = A \ [d2rdx2_2 ; drdx_1 ; r_2 ; drdx_2 ; r_3 ; drdx_3];        
        
        % Generate hub coordinates
        x = sort(linspace(x_1,x_2,round(abs((x_2 - x_1) / d_tol))))';
        xr_hub = [xr_hub ; x polyval(p,x)];
        
    end
    
    % Take unique meridional points on hub
    [~,i] = unique(xr_hub(:,1)); xr_hub = xr_hub(i,:);

    % Return only meridional curves
    clear b; b.xr_hub = xr_hub; b.xr_cas = xr_cas;
end


end

function xr = construct_quadratic(x_1,x_2,r_2,drdx_1,drdx_2,d_tol)

% Generate polynomial coefficients
A = [2*x_1 1 0 ; x_2^2 x_2 1 ; 2*x_2 1 0];
p = A \ [drdx_1 ; r_2 ; drdx_2];
   
% Generate axial coordinates
ni_hub = round(abs((x_2 - x_1) / d_tol));
x = sort(linspace(x_1,x_2,ni_hub))';

% Evaluate polynomial meridional curve
xr = [x polyval(p,x)];
    
end

function b = construct_cubic(b,d_tol)
% Construct a set of cubic meridional lines within a specified row

% Calculate hub and casing radii from mean radii and areas
r_hub_le = b.r_le - b.A_le / (4 * pi * b.r_le);
r_cas_le = b.r_le + b.A_le / (4 * pi * b.r_le);
r_hub_te = b.r_te - b.A_te / (4 * pi * b.r_te);
r_cas_te = b.r_te + b.A_te / (4 * pi * b.r_te);

% Check if axial coordinates need to be calculated from blade definition
if isfield(b,'x_hub_le') == 0
    
    % Get leading and trailing edge blade coordinates
    i_le = 1; i_te = round(size(b.xrrt,1)/2);
    xr_le = squeeze(b.xrrt(i_le,:,1:2)); xr_te = squeeze(b.xrrt(i_te,:,1:2)); 

    % Calculate intersections with hub and casing
%     x_hub_le = interp1(xr_le(:,2),xr_le(:,1),0,'pchip');
%     x_hub_te = interp1(xr_te(:,2),xr_te(:,1),0,'pchip');
%     x_cas_le = interp1(xr_le(:,2),xr_le(:,1),1,'pchip');
%     x_cas_te = interp1(xr_te(:,2),xr_te(:,1),1,'pchip');
    x_hub_le = xr_le(2,1); x_hub_te = xr_te(2,1); x_cas_le = xr_le(end-1,1); x_cas_te = xr_te(end-1,1); 
else
    x_hub_le = b.x_hub_le; x_hub_te = b.x_hub_te; x_cas_le = b.x_cas_le; x_cas_te = b.x_cas_te;
end

% Calculate polynomial coefficients for cubic hub and casing lines
A_hub = [x_hub_le^3 x_hub_le^2 x_hub_le 1 ; ...
    3*x_hub_le^2 2*x_hub_le 1 0 ;...
    x_hub_te^3 x_hub_te^2 x_hub_te 1 ; ...
    3*x_hub_te^2 2*x_hub_te 1 0];
p_hub = A_hub \ [r_hub_le ; b.drdx_hub_le ; r_hub_te ; b.drdx_hub_te];

A_cas = [x_cas_le^3 x_cas_le^2 x_cas_le 1 ; ...
    3*x_cas_le^2 2*x_cas_le 1 0 ;...
    x_cas_te^3 x_cas_te^2 x_cas_te 1 ; ...
    3*x_cas_te^2 2*x_cas_te 1 0];
p_cas = A_cas \ [r_cas_le ; b.drdx_cas_le ; r_cas_te ; b.drdx_cas_te];

% Evaluate hub and casing line at desired resolution
ni_hub = round((x_hub_te - x_hub_le) / d_tol);
ni_cas = round((x_cas_te - x_cas_le) / d_tol);

x_hub = linspace(x_hub_le,x_hub_te,ni_hub)';
x_cas = linspace(x_cas_le,x_cas_te,ni_cas)';

xr_hub = [x_hub polyval(p_hub,x_hub)];
xr_cas = [x_cas polyval(p_cas,x_cas)];

% Return new meridional curves
b.xr_hub = xr_hub; b.xr_cas = xr_cas;

end


