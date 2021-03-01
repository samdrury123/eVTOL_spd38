function [bids] = ts_get_bids(Data)
% Pass a data cell array and return bids in it.

bids = zeros(1,length(Data));
for n = 1:length(Data)
    bids(n) = Data{n}.attribute.bid;
end

end