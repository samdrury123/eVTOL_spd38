function [ g ] = ts_add_fixed_patch(bid, Q, V, g)
% Add a fixed flow patch to a dataset

% Define patch attributes and variables
p.attribute.bid = bid;
p.attribute.pid = -1;
p.attribute.kind = 4;
p.attribute.ist = Q(1,1); p.attribute.ien = Q(1,2);
p.attribute.jst = Q(2,1); p.attribute.jen = Q(2,2);
p.attribute.kst = Q(3,1); p.attribute.ken = Q(3,2);

p.pstagfixed = V(1);
p.tstagfixed = V(2);
p.vxfixed = V(3);
p.vrfixed = V(4);
p.vtfixed = V(5);

% Add patch to dataset
g = ts_add_patch(bid, p, g);

end