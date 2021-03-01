load('EXP_DATA_PROP.mat');
load('APC_FIT.mat', 'APC_FIT');
%%
rpmlim = [3000 12000];

testrpm = linspace(rpmlim(1), rpmlim(2), 10000);
testthr = APC_FIT.rpm_thr(testrpm);

for i = 1:length(APC.P)
    apcthr = ones([1, 10000]) .* APC.T(i);
    diff = apcthr - testthr;
    for j = 1:10000
        diffj = diff(j);
        if diffj < 0
            break
        end
    end
    rpm(i) = testrpm(j);
end

APC.rpmmean = rpm;

save('EXP_DATA_PROP.mat', 'APC', 'BASE');

clear rpm diff diffj apcthr testthr testrpm rpmlim APC_FIT i j