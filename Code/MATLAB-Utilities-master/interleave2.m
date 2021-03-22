function xr = interleave2(xr,ntimes)
% Interleaves more points in two dimensions by ntimes

for j = 1:ntimes
    xr_temp = zeros(size(xr,1)-1,2);
    for n = 1:size(xr,1) - 1
        xr_temp(n,:) = (xr(n,:) + xr(n+1,:))/2;
    end

    xr_new = zeros(size(xr,1)*2-1,2);
    
    for n = 1:size(xr,1)-1
        xr_new(2*n-1:2*n,:) = [xr(n,:) ; xr_temp(n,:)];
    end
    xr_new(end,:) = xr(end,:);
    
    xr = xr_new;
end

end