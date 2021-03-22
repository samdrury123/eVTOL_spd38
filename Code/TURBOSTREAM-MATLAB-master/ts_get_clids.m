function [clids] = ts_get_clids(varargin)
% Pass a block cell array, or bid and Data, and return clids in it.
if length(varargin) == 1
    b = varargin{1};
elseif length(varargin) == 2
    bid = varargin{1};
    Data = varargin{2};
    b = ts_get_block(bid, Data);
end

if isfield(b, 'clist') == 1
    clids = zeros(1,length(b.clist));
    for n = 1:length(b.clist)
        clids(n) = b.clist{n}.attribute.clid;
    end
else
    clids = [];
end