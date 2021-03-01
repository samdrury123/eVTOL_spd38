function n = unit(n)
% UNIT  Scale a list or matrix of vectors to unit length
%
%   n = UNIT(n);
%
%   n - array of any number of dimensions, last dimension is vector length in each direction

% Divide components of vectors by current magnitude
n = n ./ repmat(sum(n.^2,ndims(n)).^0.5,[ones(1,ndims(n)-1) size(n,ndims(n))]);

end

