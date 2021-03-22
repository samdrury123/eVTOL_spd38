function xrrt = bl_cart2pol(xyz)
% BL_CART2POL  Convert a set of cartesian blade coordinates to polars
%
%   xrrt = BL_CART2POL(xyz)
%
%   xyz - input array or cell array of cartesian blade coordinates
%   xrrt - output array or cell array of polar blade coordinates
%
%   Coordinates are spread across the final array dimension


% Check whether to operate on a cell array or numeric array
if iscell(xyz) == 1
    
    % Loop over all section cells and call bl_cart2pol again
    xrrt = cell(size(xyz));
    for n = 1:length(xyz)
        xrrt{n} = bl_cart2pol(xyz{n});
    end
    
elseif isnumeric(xyz) == 1 && size(xyz,ndims(xyz)) == 3
    
    % Record size of input array and reshape into a list
    nn = size(xyz); xyz = reshape(xyz,[prod(nn(1:end-1)) 3]);
    
    % Calculate polar coordinates
    [t, r] = cart2pol(xyz(:,2), xyz(:,3));
    xrrt = [xyz(:,1) r r.*t];
    
    % Reshape back to original size
    xrrt = reshape(xrrt,nn); 
    
else
    
    % Display error for unrecognised types
    error('Unrecognised input');
    
end


end

