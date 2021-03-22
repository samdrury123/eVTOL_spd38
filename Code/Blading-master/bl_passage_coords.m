function [sr,sx,st] = g_passage_coords(xrrt,xr_hub,xr_cas,xr_le,xr_te,xrrt_ps,xrrt_ss)
% G_PASSAGE_COORDS  Calculate passage relative coordinates for input polar coordinates

% Record size of input coordinates and reshape into a list
nn = size(xrrt); xrrt = reshape(xrrt,[prod(nn(1:end-1)) 3]);

% Calculate non-dimensional radius
r_hub = interp1(xr_hub(:,1),xr_hub(:,2),xrrt(:,1),'linear','extrap');
r_cas = interp1(xr_cas(:,1),xr_cas(:,2),xrrt(:,1),'linear','extrap');
sr = (xrrt(:,2) - r_hub) ./ (r_cas - r_hub);

% Check if arguements for non-dimensional axial coordinate are specified
if exist('xr_le','var') ~= 0

    % Non-dimensional radius of leading and trailing edge lines
    r_hub = interp1(xr_hub(:,1),xr_hub(:,2),xr_le(:,1),'linear','extrap');
    r_cas = interp1(xr_cas(:,1),xr_cas(:,2),xr_le(:,1),'linear','extrap');
    sr_le = (xr_le(:,2) - r_hub) ./ (r_cas - r_hub);
    r_hub = interp1(xr_hub(:,1),xr_hub(:,2),xr_te(:,1),'linear','extrap');
    r_cas = interp1(xr_cas(:,1),xr_cas(:,2),xr_te(:,1),'linear','extrap');
    sr_te = (xr_te(:,2) - r_hub) ./ (r_cas - r_hub);

    % Calculate non-dimensional axial coordinate
    x_le = interp1(sr_le,xr_le(:,1),sr,'linear','extrap');
    x_te = interp1(sr_te,xr_te(:,1),sr,'linear','extrap');
    sx = (xrrt(:,1) - x_le) ./ (x_te - x_le);
    
    % Check if arguments for non-dimensional tangential coordinate are specified
    if exist('xrrt_ps','var') ~= 0

        % Non-dimensional radius of surfaces
        r_hub = interp1(xr_hub(:,1),xr_hub(:,2),xrrt_ps(:,:,1),'linear','extrap');
        r_cas = interp1(xr_cas(:,1),xr_cas(:,2),xrrt_ps(:,:,1),'linear','extrap');
        sr_ps = (xrrt_ps(:,:,2) - r_hub) ./ (r_cas - r_hub);
        r_hub = interp1(xr_hub(:,1),xr_hub(:,2),xrrt_ss(:,:,1),'linear','extrap');
        r_cas = interp1(xr_cas(:,1),xr_cas(:,2),xrrt_ss(:,:,1),'linear','extrap');
        sr_ss = (xrrt_ss(:,:,2) - r_hub) ./ (r_cas - r_hub);

        % Non-dimensional axial coordinate of surfaces
        x_le = interp1(sr_le,xr_le(:,1),sr_ps,'linear','extrap');
        x_te = interp1(sr_te,xr_te(:,1),sr_ps,'linear','extrap');
        sx_ps = (xrrt_ps(:,:,1) - x_le) ./ (x_te - x_le);
        x_le = interp1(sr_le,xr_le(:,1),sr_ss,'linear','extrap');
        x_te = interp1(sr_te,xr_te(:,1),sr_ss,'linear','extrap');
        sx_ss = (xrrt_ss(:,:,1) - x_le) ./ (x_te - x_le);
        
        % Check non-dimensional axial coordinates
%         if any(abs(sx_ps(:) - sx_ss(:)) > 1e-6) == 1
%             error('Difference in surface axial coordinates');
%         end
        
        % Assume surfaces area already plaid and are at the same non-dimensional meridional coords
        i_mid = round(size(sx_ps,1) / 2); j_mid = round(size(sx_ps,2) / 2);
        sx_surf = repmat(sx_ps(:,j_mid),[1 size(sx_ps,2)]); 
        sr_surf = repmat(sr_ps(i_mid,:),[size(sr_ps,1) 1]); 
        
        % Calculate surface tangential coordinates
        f_ps = griddedInterpolant(sx_surf,sr_surf,xrrt_ps(:,:,3),'linear','linear');
        f_ss = griddedInterpolant(sx_surf,sr_surf,xrrt_ss(:,:,3),'linear','linear');
        rt_ps = f_ps(sx,sr); rt_ss = f_ss(sx,sr);
%         rt_ps = interp2(
%         rt_ss = interp2(sx_surf',sr_surf',xrrt_ss(:,:,3)',sx,sr,'linear');
        
        % Calculate non-dimensional tangential coordinate
        st = (xrrt(:,3) - rt_ss) ./ (rt_ps - rt_ss);
           
    else
        st = [];
    end
    
else
    
    sx = []; st = [];

end

% Reshape back to original size
sr = reshape(sr,[nn(1:end-1) 1]); 
if isempty(sx) == 0; sx = reshape(sx,[nn(1:end-1) 1]); end; 
if isempty(st) == 0; st = reshape(st,[nn(1:end-1) 1]); end;


end

