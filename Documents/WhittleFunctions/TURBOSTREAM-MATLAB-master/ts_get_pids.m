function [pids] = ts_get_pids(varargin)
% Pass a block cell array, or bid and Data, and return pids in it.
if length(varargin) == 1
    b = varargin{1};
elseif length(varargin) == 2
    if iscell(varargin{1})
        Data = varargin{1};
        bid = varargin{2};
    else
        bid = varargin{1};
        Data = varargin{2};
    end
    b = ts_get_block(bid, Data);
end

if isfield(b, 'patch') == 1 && numel(b.patch) > 0
    pids = zeros(1,length(b.patch));
    for n = 1:length(b.patch)
    pids(n) = b.patch{n}.attribute.pid;
    end
else
    pids = [];
end

end