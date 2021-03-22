function [varout] = var2matrix(var)
% Function that takes row of elements from 4:22 and converts them into a
% 5x4 matrix suitable for producing contour plots

% Check input array is right size
if numel(var) ~= 19
    error('Incorrect variable size - need 19 elements')
end

% Initialise empty NaN array
mat = NaN.*ones(5,4);
mat(:,4) = var(15:19);
mat(:,3) = var(10:14);
mat(:,2) = var(5:9);
mat(1:4,1) = var(1:4);

varout = mat;
%% chic