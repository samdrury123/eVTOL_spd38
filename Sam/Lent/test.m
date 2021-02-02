clear; close;

philist = [0.5:0.1:0.8];
sigmalist = [0.8:0.1:1.2];
Patient([4 5]) = struct();
Patient(2,3).g = 5;
for pp=1:size(philist,2)
    for ss=1:size(sigmalist,2)
        patient(pp,ss).phi = philist(pp);
        patient(pp,ss).sigma = sigmalist(ss);
    end
end
alt = 1000;
atm = Altitude(alt); % At sea level