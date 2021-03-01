function [xrrt_1,xrrt_2,i_1,i_2] = bl_split(xrrt,i_le,i_te)
% BL_SPLIT  Split a blade or section definition into two halves

% Length of array
ni = size(xrrt,1);

% Check orientation of split indices
if i_te > i_le
    i_1 = [i_le:-1:2 ni:-1:i_te];
    i_2 = i_le:i_te;
else
    i_1 = i_le:-1:i_te;
    i_2 = [i_le:ni 2:i_te];
end

% Extract coordinates
xrrt_1 = xrrt(i_1,:,:);
xrrt_2 = xrrt(i_2,:,:);

end

