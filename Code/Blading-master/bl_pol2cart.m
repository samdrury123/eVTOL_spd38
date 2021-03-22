function xyz = bl_pol2cart(xrrt,rt_input)
% BL_POL2CART  Convert a set of polar blade coordinates to cartesian
%
%   xyz = BL_POL2CART(xrrt,rt_input)
%
%   xrrt - input array or cell array of polar blade coordinates
%   rt_input - optional input 0 or 1 to determine whether r-theta or theta is input
%   xyz - output array or cell array of cartesian blade coordinates
%
%   Coordinates are spread across the final array dimension

% Default to input rt
if exist('rt_input','var') == 0; rt_input = 1; end;

% Check whether to operate on a cell array or numeric array
if iscell(xrrt) == 1
    
    % Loop over all section cells and call bl_cart2pol again
    xyz = cell(size(xrrt));
    for n = 1:length(xrrt)
        xyz{n} = bl_pol2cart(xrrt{n},rt_input);
    end
    
elseif isnumeric(xrrt) == 1 && size(xrrt,ndims(xrrt)) == 3
    
    % Record size of input array and reshape into a list
    nn = size(xrrt); xrrt = reshape(xrrt,[prod(nn(1:end-1)) 3]);
    
    % Get tangential coordinates
    if rt_input == 1 
        t = xrrt(:,3) ./ xrrt(:,2);
    else
        t = xrrt(:,3);
    end
    
    % Calculate cartesian coordinates
    [y, z] = pol2cart(t, xrrt(:,2));
    xyz = [xrrt(:,1) y z];

    % Reshape back to original size
    xyz = reshape(xyz,nn); 
    
else
    
    % Display error for unrecognised types
    error('Unrecognised input');
    
end


end

