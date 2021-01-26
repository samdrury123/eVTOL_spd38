open('EXP_DATA.mat');

figure(101); hold on; plot(sig08.RPM, sig08.mass, ':');plot(sig10.RPM, sig10.mass, '-');plot(sig12.RPM, sig12.mass, '-.');  title('Mass vs RPM'); xlabel('RPM'); ylabel('Mass'); legend('Sigma = 0.8', 'Sigma = 1.0', 'Sigma = 1.2');
figure(102); hold on; plot(sig08.RPM, sig08.FOM, ':');plot(sig10.RPM, sig10.FOM, '-');plot(sig12.RPM, sig12.FOM, '-.');  title('FOM vs RPM'); xlabel('RPM'); ylabel('FOM'); legend('Sigma = 0.8', 'Sigma = 1.0', 'Sigma = 1.2');
figure(103); hold on; plot(sig08.RPM, sig08.P, ':');plot(sig10.RPM, sig10.P, '-');plot(sig12.RPM, sig12.P, '-.');  title('P vs RPM'); xlabel('RPM'); ylabel('P'); legend('Sigma = 0.8', 'Sigma = 1.0', 'Sigma = 1.2');