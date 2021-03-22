function g = ts_add_periodic_patch(bid0, bid1, Q0, Q1, dir0, dir1, g, kind)
% Make two periodic patches and add to correct blocks

if exist('kind','var') == 0
    kind = 5;
end

% Get blocks and patch ids
b0 = ts_get_block(bid0,g);
b1 = ts_get_block(bid1,g);

[pids1] = ts_get_pids(b1);

% Create attributes for first patch
attribute0.bid = bid0;
attribute0.pid = -1;    % Assigned when patch is added
attribute0.kind = kind;
attribute0.ist = int32(Q0(1,1)); attribute0.ien = int32(Q0(1,2));
attribute0.jst = int32(Q0(2,1)); attribute0.jen = int32(Q0(2,2));
attribute0.kst = int32(Q0(3,1)); attribute0.ken = int32(Q0(3,2));
attribute0.idir = dir0(1); attribute0.jdir = dir0(2); attribute0.kdir = dir0(3);
attribute0.nxbid = bid1;
if isempty(pids1) == 0 && bid0~=bid1
    attribute0.nxpid = pids1(end) + 1;
elseif bid0 == bid1 && isempty(pids1) == 0
    attribute0.nxpid = pids1(end) + 2;
elseif bid0 == bid1 && isempty(pids1) == 1
    attribute0.nxpid = 1;
else
    attribute0.nxpid = 0;
end
attribute0.nface = 0;
attribute0.nt = 1;

p0.attribute = attribute0;

% Add patch 0 to correct block in Data
[ g ] = ts_add_patch(bid0, p0, g);

% Get pids again in case same block is patched
b0 = ts_get_block(bid0,g);
[pids0] = ts_get_pids(b0);

% Create attributes for second patch
attribute1.bid = bid1;
attribute1.pid = -1;
attribute1.kind = kind;
attribute1.ist = int32(Q1(1,1)); attribute1.ien = int32(Q1(1,2));
attribute1.jst = int32(Q1(2,1)); attribute1.jen = int32(Q1(2,2));
attribute1.kst = int32(Q1(3,1)); attribute1.ken = int32(Q1(3,2));
attribute1.idir = dir1(1); attribute1.jdir = dir1(2); attribute1.kdir = dir1(3);
attribute1.nxbid = bid0;
attribute1.nxpid = pids0(end);
attribute1.nface = 0;
attribute1.nt = 1;

p1.attribute = attribute1;

% Add patch 1 to correct block in Data
[ g ] = ts_add_patch(bid1, p1, g);

end