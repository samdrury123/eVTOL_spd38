clear; close;
philist = [0.5:0.1:0.8];
sigmalist = [0.8:0.1:1.2];
i=0;
for pp=1:size(philist,2)
    for ss=1:size(sigmalist,2)
        i=i+1;
        patient(pp,ss).phi = philist(pp);
        patient(pp,ss).sigma = sigmalist(ss);
    end
end

