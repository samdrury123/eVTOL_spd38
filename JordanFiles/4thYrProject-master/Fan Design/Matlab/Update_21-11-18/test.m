r_h = 20e-3;
r_c = 60e-3;
n = 10;
rtest = zeros(n,n,n);
for i=1:n
    for j = 1:n
    rtest(i,j,:) = linspace(r_h, r_c, n);
    end
end

rtest(1, 1, :)
rtest(2, 1, :)


for i = 1:n
            for j = 1:n
                r_cMESH(i,j,:) = linspace(r_h,r_cMAX, n);
            end
        end