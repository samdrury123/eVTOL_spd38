function [B] = bernstMatrix(n,x)
% Calculate Bernstein matrix (n
%                             k) x^k*(1-x)^(n-k)
for k=0:n
    B(:,k+1) = nchoosek(n,k)*x.^k.*(1-x).^(n-k);
end
end