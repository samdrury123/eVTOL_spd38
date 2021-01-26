function [] = OptimiseMassMESH()
%
%
%% Load geom
global n

load('meshvars.mat','geom');
load('vars.mat', 'powerrated', 'W', 'kv', 'maxcurrent');

% powerrated = maxcurrent*15; % Rated power for const omega set by current limit, less than rated power

phi=geom.phi;
psi=geom.psi;
t=geom.t;
r_h=geom.r_h;
r_cMAX=geom.r_cMAX ;
turboL=geom.turboL ;
inletL=geom.inletL ;
omega=geom.omega ;
rhoPLA=geom.rhoPLA ;
rhoAIR=geom.rhoAIR ;
variables=geom.variables;
plott=geom.plott;
totalmass=geom.totalmass;
numberofprop = geom.numberofprop;
sig=geom.sig;
%% Meshgrid Variables
disp('INITIALISING meshes...');
switch variables
    case 'rcrh'
        [r_cMESH, r_hMESH] = meshgrid(r_crange, r_hrange);
        omegaMESH = omega.*ones(n,n);
        tMESH = t.*ones(n,n);
        phiMESH = phi.*ones(n,n);
        psiMESH = psi.*ones(n,n);
        xplot = r_cMESH;
        yplot = r_hMESH;
        sigmaMESH = sqrt(phiMESH.^2./(2*psiMESH));
    case 'r_ct'
        [r_cMESH, tMESH] = meshgrid(r_crange, trange);
        omegaMESH = omega.*ones(n,n);
        r_hMESH = r_h.*ones(n,n);
        phiMESH = phi.*ones(n,n);
        psiMESH = psi.*ones(n,n);
        xplot = r_cMESH;
        yplot = tMESH;
        sigmaMESH = sqrt(phiMESH.^2./(2*psiMESH));
    case 'r_comega'
        [r_cMESH, omegaMESH] = meshgrid(r_crange, omegarange);
        tMESH = t.*ones(n,n);
        r_hMESH = r_h.*ones(n,n);
        phiMESH = phi.*ones(n,n);
        psiMESH = psi.*ones(n,n);
        xplot = r_cMESH;
        yplot = omegaMESH;
        sigmaMESH = sqrt(phiMESH.^2./(2*psiMESH));
    case 'phipsi'
        phiMESH = zeros(n,n,n);
        psiMESH = zeros(n,n,n);
        r_cMESH = zeros(n,n,n);
        for i = 1:n
            for j = 1:n
                r_cMESH(i,j,:) = linspace(r_h,r_cMAX, n);
                phiMESH(:,i,j) = linspace(0,1, n);
                psiMESH(i,:,j) = linspace(0,0.5, n);
            end
        end
        
        tMESH = zeros(n,n,n) + t;
        r_hMESH = zeros(n,n,n) + r_h;
        omegaMESH = zeros(n,n,n) + omega;
        turboLMESH = zeros(n,n,n) + turboL;
        inletLMESH = zeros(n,n,n) + inletL;
        xplot = phiMESH(:,:,1);
        yplot = psiMESH(:,:,1);
        sigmaMESH = sqrt(phiMESH.^2./(2*psiMESH));
        LMESH = (r_cMESH-r_hMESH)./DF(phiMESH,psiMESH);
        r_icMESH = r_cMESH + (r_cMESH-r_hMESH).*0.25;
        r_ihMESH = r_hMESH - (r_cMESH-r_hMESH).*0.25;
        r_ecMESH = 0.5.*(r_cMESH.*(1+sigmaMESH) + r_hMESH.*(1-sigmaMESH));
        r_ehMESH = 0.5.*(r_cMESH.*(1-sigmaMESH) + r_hMESH.*(1+sigmaMESH));
        
%         [phiMESH, psiMESH] = meshgrid(phirange, psirange);
%         tMESH = t.*ones(n,n);
%         r_hMESH = r_h.*ones(n,n);
%         r_cMESH = r_c.*ones(n,n);
%         omegaMESH = omega.*ones(n,n);
%         turboLMESH = turboL.*ones(n,n);
%         inletLMESH = inletL.*ones(n,n);
%         xplot = phiMESH;
%         yplot = psiMESH;
%         sigmaMESH = sqrt(phiMESH.^2./(2*psiMESH));
%         LMESH = (r_cMESH-r_hMESH)./DFMESH(phiMESH,psiMESH);
%         r_icMESH = r_cMESH + (r_cMESH-r_hMESH).*0.25;
%         r_ihMESH = r_hMESH - (r_cMESH-r_hMESH).*0.25;
%         r_ecMESH = 0.5.*(r_cMESH.*(1+sigmaMESH) + r_hMESH.*(1-sigmaMESH));
%         r_ehMESH = 0.5.*(r_cMESH.*(1-sigmaMESH) + r_hMESH.*(1+sigmaMESH));
    case 'phisigma'
        [phiMESH, sigmaMESH] = meshgrid(phirange, sigmarange);
        tMESH = t.*ones(n,n);
        r_hMESH = r_h.*ones(n,n);
        r_cMESH = r_c.*ones(n,n);
        omegaMESH = omega.*ones(n,n);
        psiMESH = phiMESH.^2 ./ (2*sigmaMESH.^2);
        xplot = phiMESH;
        yplot = sigmaMESH;
        
    case 'phipsi_omega'
        phiMESH = zeros(n,n,n);
        psiMESH = zeros(n,n,n);
        omegaMESH = zeros(n,n,n);
        for i = 1:n
            for j = 1:n
                omegaMESH(i,j,:) = linspace(1000*pi/30,20000*pi/30, n);
                phiMESH(:,i,j) = linspace(0,1, n);
                psiMESH(i,:,j) = linspace(0,0.5, n);
            end
        end
        
        tMESH = zeros(n,n,n) + t;
        r_hMESH = zeros(n,n,n) + r_h;
        r_cMESH = zeros(n,n,n) + 60e-3;
        turboLMESH = zeros(n,n,n) + turboL;
        inletLMESH = zeros(n,n,n) + inletL;
        xplot = phiMESH(:,:,1);
        yplot = psiMESH(:,:,1);
        sigmaMESH = sqrt(phiMESH.^2./(2*psiMESH));
        LMESH = (r_cMESH-r_hMESH)./DF(phiMESH,psiMESH);
        r_icMESH = r_cMESH + (r_cMESH-r_hMESH).*0.25;
        r_ihMESH = r_hMESH - (r_cMESH-r_hMESH).*0.25;
        r_ecMESH = 0.5.*(r_cMESH.*(1+sigmaMESH) + r_hMESH.*(1-sigmaMESH));
        r_ehMESH = 0.5.*(r_cMESH.*(1-sigmaMESH) + r_hMESH.*(1+sigmaMESH));
        r_c = 60e-3;
end

%% MESH CALCULATIONS
disp('BEGINNING mesh calculations...');
tic;
fprintf('- Propulsor mass');
% Mass of propulsor, excluding blades, including motor
casingmassMESH = rhoPLA*pi*turboLMESH.*((r_cMESH + t).^2 - (r_cMESH).^2);
hubmassMESH = rhoPLA*pi*turboLMESH.*(r_hMESH.^2 - (r_hMESH - 10e-3).^2);
unitmassMESH = casingmassMESH + hubmassMESH;
fprintf('    DONE \n');

fprintf('- Diffusor mass');
% Diffuser Casing
ch1MESH = LMESH.*(r_ecMESH + t) ./ (r_ecMESH - r_cMESH);
ch2MESH = LMESH.*(r_ecMESH) ./ (r_ecMESH - r_cMESH);
volcaseMESH = (pi/3)* (ch1MESH.*(r_ecMESH + t).^2 - (ch1MESH-LMESH).*(r_cMESH - t).^2 - ch2MESH.*r_ecMESH.^2 + (ch2MESH-LMESH).*r_cMESH.^2);

% Diffuser Hub 
hh1MESH = LMESH.*(r_hMESH) ./ (r_hMESH - r_ehMESH);
hh2MESH = LMESH.*(r_hMESH - t) ./ (r_hMESH - r_ehMESH);
volhubMESH = (pi/3)* (hh1MESH.*r_hMESH.^2 - (hh1MESH-LMESH).*(r_ehMESH).^2 - hh2MESH.*(r_hMESH - t).^2 + (hh2MESH-LMESH).*(r_ehMESH - t).^2);
fprintf('    DONE \n');

fprintf('- Inlet mass');
% Inlet Casing
ich1MESH = inletLMESH.*(r_icMESH + t) ./ (r_icMESH - r_cMESH);
ich2MESH = inletLMESH.*(r_icMESH) ./ (r_icMESH - r_cMESH);
volicaseMESH = (pi/3)* (ich1MESH.*(r_icMESH + t).^2 - (ich1MESH-inletLMESH).*(r_cMESH - t).^2 - ich2MESH.*r_icMESH.^2 + (ich2MESH-inletLMESH).*r_cMESH.^2);

% Inlet Hub 
ihh1MESH = inletLMESH.*(r_hMESH) ./ (r_hMESH - r_ihMESH);
ihh2MESH = inletLMESH.*(r_hMESH - t) ./ (r_hMESH - r_ihMESH);
volihubMESH = (pi/3)* (ihh1MESH.*r_hMESH.^2 - (ihh1MESH-inletLMESH).*(r_ihMESH).^2 - ihh2MESH.*(r_hMESH - t).^2 + (ihh2MESH-inletLMESH).*(r_ihMESH - t).^2);
% 
% volconeMESH = pi*r_hMESH.^2.*(r_cMESH - r_hMESH)./(3.*DF(phiMESH,psiMESH));
% volshellMESH = 2*pi*tMESH.*sqrt(phiMESH.^2.*(r_cMESH.^2 - r_hMESH.^2)./2.*psiMESH).*sqrt(((r_cMESH - r_hMESH)./DF(phiMESH, psiMESH)).^2 + ((r_cMESH.^2 - r_hMESH.^2).*phiMESH.^2./2.*psiMESH - r_cMESH).^2);
volcaseMESH(imag(volcaseMESH)~=0) = NaN;
volhubMESH(imag(volhubMESH)~=0) = NaN;

exitmassMESH = rhoPLA*(volcaseMESH + volhubMESH);
inletmassMESH = rhoPLA*(volicaseMESH + volihubMESH);
fprintf('    DONE \n');

% Total mass
fprintf('- Total mass');
massMESH = inletmassMESH + exitmassMESH + unitmassMESH;
fprintf('    DONE \n');

% Thrust
fprintf('- Thrust');
%     thrust(i) = phi*sqrt(psi)*rhoAIR*pi*omega^2*((r_crange(i) - r_h)*(r_crange(i) + r_h)^3)/sqrt(2);
thrustMESH = phiMESH.^2.*rhoAIR*pi.*(r_cMESH.^2 - r_hMESH.^2).*(omegaMESH.^2).*((r_cMESH + r_hMESH)./2).^2 ./ sigmaMESH;
fprintf('    DONE \n');

% Thrust Excess
fprintf('- Thrust excess');
deficitMESH = thrustMESH./((massMESH + W.body/numberofprop)*9.81);
fprintf('    DONE \n');
% idx = find(abs(deficitMESH-1)<.003);

% Power
fprintf('- Power');
powerMESH = phiMESH.*psiMESH.*rhoAIR.*pi.*(r_cMESH.^2 - r_hMESH.^2).*(omegaMESH.^3).*((r_cMESH + r_hMESH)./2).^3;
fprintf('    DONE \n');

% Power with Mass
% powerMESH = ((massMESH + totalmass/numberofprop)*9.81).^(3/2) ./ sqrt(4.*sigmaMESH.*rhoAIR.*omegaMESH.*((r_cMESH + r_hMESH)./2));

fprintf('- Figure of merit');
% FoM
meritMESH = (thrustMESH./powerMESH).*sqrt(thrustMESH./(2*rhoAIR*pi.*(r_cMESH.^2 - r_hMESH.^2)));
fprintf('    DONE \n');

switch variables
    case 'phipsi'
        % Critical casing radius
        disp('');
        fprintf('Detrmining critical casing radius for each design point... ');
        for i = 1:n
            for j = 1:n
                for k = 1:n
                    if deficitMESH(i,j,k) > 1
                        r_cCRITMESH(i,j) = r_cMESH(i,j,k);
                        OPthrustMESH(i,j) = thrustMESH(i,j,k);
                        OPpowerMESH(i,j) = powerMESH(i,j,k);
                        OPmeritMESH(i,j) = meritMESH(i,j,k);
                        OPmassMESH(i,j) = massMESH(i,j,k);
                        critidx(i,j) = k;
                        break
                    elseif k == n
                        r_cCRIT(i,j) = r_cMAX;
                    end
                end
            end
        end
        
        fprintf('    DONE \n');

        % Determine if critical sig greater than required sig
        fprintf('Detrmining critical diffusion factor for each design point... ');
        psiCRITMESH = (0.5*xplot.^2)./((r_cCRITMESH + r_h)./(r_cCRITMESH - r_h));
        psitest = abs(yplot./psiCRITMESH - 1);
        
    case 'phipsi_omega'
        % Critical omega
        disp('');
        fprintf('Detrmining critical motor speed for each design point... ');
        for i = 1:n
            for j = 1:n
                for k = 1:n
                    if deficitMESH(i,j,k) > 1
                        omegaCRITMESH(i,j) = omegaMESH(i,j,k);
                        OPthrustMESH(i,j) = thrustMESH(i,j,k);
                        OPpowerMESH(i,j) = powerMESH(i,j,k);
                        OPmeritMESH(i,j) = meritMESH(i,j,k);
                        OPmassMESH(i,j) = massMESH(i,j,k);
                        critidx(i,j) = k;
                        break
                    elseif k == n
                        omegaCRIT(i,j) = 20000*pi/30;
                    end
                end
            end
        end
        fprintf('    DONE \n');

        % Determine if critical sig greater than required sig
        fprintf('Detrmining critical diffusion factor for each design point... ');
        psiCRITMESH = (0.5*xplot.^2)./((r_c + r_h)./(r_c - r_h));
        psitest = abs(yplot./psiCRITMESH - 1);
end
        


for ii = 1:n
    for jj = 1:n
        if ~isnan(psitest(ii,jj))
            if psitest(ii,jj) < 0.01
                psilim(ii) = yplot(ii,jj);
                break
            end
        end
    end
end


OPthrustMESH(OPthrustMESH==0) = NaN;
OPpowerMESH(OPpowerMESH==0) = NaN;
OPmeritMESH(OPmeritMESH==0) = NaN;
OPmassMESH(OPmassMESH==0) = NaN;
psilim(psilim==0) = NaN;
fprintf('    DONE \n');

% Determine power/current limit
fprintf('Detrmining torque limit... ');

switch variables
    case 'phipsi'
        r_cCRITMESH(r_cCRITMESH==0) = NaN;
        nondimpowerlim = (30*maxcurrent)*omegaMESH/(pi*kv*powerrated);
    case 'phipsi_omega'
        omegaCRITMESH(omegaCRITMESH==0) = NaN;
        nondimpowerlim = (30*maxcurrent)*omegaCRITMESH/(pi*kv*powerrated);
end 

for id = 1:n
    for jd = 1:n
        testp = OPpowerMESH(id,jd)/powerrated;
        if testp > nondimpowerlim(id,jd) && ~isnan(prev)
            powerlim(id) = yplot(id,jd);
            break
        elseif testp > nondimpowerlim(id,jd)
            powerlim(id) = NaN;
            break
        elseif jd == n && isnan(testp)
            powerlim(id) = NaN;
            break
        else
            prev = testp;
        end
    end
end
fprintf('    DONE \n');

disp('---END MESH CALC---');

% Plot FoM as a function of only sig, show how it is equal to the above
% version
% meritMESH = sqrt(2.*sigmaMESH);

meshtime = toc;
disp('Time to execute mesh calculations =');disp(meshtime);

system('rm OPMESH*');

%% Save all 
switch variables
    case 'phipsi'
        save('OPMESH','r_cCRITMESH', 'OPthrustMESH', 'OPmassMESH', 'OPpowerMESH', 'OPmeritMESH', 'psilim', 'plott','variables','xplot','yplot');
    case 'phipsi_omega'
        save('OPMESH','omegaCRITMESH', 'OPthrustMESH', 'OPmassMESH', 'OPpowerMESH', 'OPmeritMESH', 'psilim', 'plott','variables','xplot','yplot');
end
        %% Clear memory
clear *MESH
%% Reload data
load('OPMESH.mat');
%% Plot MESH

if plott==1
    switch variables
        case 'rcrh'
            figure(108); [C1,h1] = contourf(xplot, yplot, powerMESH); clabel(C1,h1,0); title('Power r_c vs r_h');xlabel('Casing Radius r_c'); ylabel('Hub Radius r_h'); colorbar();
            figure(109); [C1,h1] = contourf(xplot, yplot, thrustMESH); clabel(C1,h1,0); title('Thrust r_c vs r_h');xlabel('Casing Radius r_c'); ylabel('Hub Radius r_h'); colorbar();
            figure(110); [C1,h1] = contourf(xplot, yplot, massMESH); clabel(C1,h1,0); title('Mass r_c vs r_h');xlabel('Casing Radius r_c'); ylabel('Hub Radius r_h');colorbar();
            figure(111); [C2,h2] = contourf(xplot, yplot, meritMESH); colorbar();
            title('Figure of Merit r_c vs r_h');
            xlabel('Casing Radius r_c'); ylabel('Hub Radius r_h');
        case 'r_comega'
            figure(2); [C1,h1] = contourf(xplot, yplot, powerMESH); clabel(C1,h1,0); title('Power r_c vs r_h');xlabel('Casing Radius r_c'); ylabel('Speed \Omega'); colorbar();
            figure(4); [C1,h1] = contourf(xplot, yplot, thrustMESH); clabel(C1,h1,0); title('Thrust r_c vs r_h');xlabel('Casing Radius r_c'); ylabel('Speed \Omega'); colorbar();
            figure(5); [C1,h1] = contourf(xplot, yplot, massMESH); clabel(C1,h1,0); title('Mass r_c vs r_h');xlabel('Casing Radius r_c'); ylabel('Speed \Omega');colorbar();
            figure(3); [C2,h2] = contourf(xplot, yplot, meritMESH); colorbar();
            title('Figure of Merit r_c vs r_h');
            xlabel('Casing Radius r_c'); ylabel('Hub Radius r_h');
        case 'phipsi'
            figure(112); hold on; [C1,h1] = contourf(xplot, yplot, r_cCRITMESH,'LineColor','none'); clabel(C1,h1,0); title('Casing Radius \phi vs \psi');xlabel('Flow Coefficient (\phi)'); ylabel('Stage Loading (\psi)'); colorbar(); ylim([0.0 0.5]);legend('Casing Radius','Location', 'northwest');
            plot(linspace(0.0, 1, 100), 0.5*linspace(0.0, 1, 100).^2,'k', 'LineWidth', 2); 
            plot(xplot, psilim, '-.k', 'LineWidth', 2); 
            plot(xplot,powerlim,'-r', 'LineWidth', 2); 
            legend('Casing Radius', 'Diffuser = 1', 'Symmetric Diffuser limit', 'Current Limit', 'Location', 'northwest');
            %             plot(linspace(0.0, 1, 100), 0.5*linspace(0.0, 1, 100).^2,'-.r');plot(xplot, psilim, '--r');
            
            figure(108); hold on; [C1,h1] = contourf(xplot, yplot, OPthrustMESH/W.body,'LineColor','none'); clabel(C1,h1,0); title('Thrust \phi vs \psi');xlabel('Flow Coefficient (\phi)'); ylabel('Stage Loading (\psi)'); colorbar();
            plot(linspace(0.0, 1, 100), 0.5*linspace(0.0, 1, 100).^2,'-.r');ylim([0.0 0.5]);plot(xplot, psilim, '--r'); legend('Thrust', 'Diffuser = 1', 'Symmetric Diffuser limit', 'Location', 'northwest');
            
            figure(109); hold on; [Cp,hp] = contourf(xplot, yplot, OPpowerMESH/powerrated,'LineColor','none'); title('Power \phi vs \psi'); xlabel('Flow Coefficient (\phi)'); ylabel('Stage Loading (\psi)'); colorbar(); ylim([0.0 0.5]);
            plot(linspace(0.0, 1, 100), 0.5*linspace(0.0, 1, 100).^2,'k', 'LineWidth', 2); 
            plot(xplot, psilim, '-.k', 'LineWidth', 2); 
            plot(xplot,powerlim,'-r', 'LineWidth', 2);
            legend('Power', 'Diffuser = 1', 'Symmetric Diffuser limit', 'Current Limit', 'Location', 'northwest');
            
            figure(110); hold on; [C1,h1] = contourf(xplot, yplot, OPmassMESH,'LineColor','none'); clabel(C1,h1,0); title('Delta Mass \phi vs \psi');xlabel('Flow Coefficient (\phi)'); ylabel('Stage Loading (\psi)'); colorbar();
            plot(linspace(0.0, 1, 100), 0.5*linspace(0.0, 1, 100).^2,'-.r');ylim([0.0 0.5]);plot(xplot, psilim, '--r');legend('Mass', 'Diffuser = 1', 'Symmetric Diffuser limit', 'Location', 'northwest');
            
            figure(111); hold on; contourf(xplot, yplot, OPmeritMESH, linspace(0,3,9)); colorbar();plot(linspace(0.0, 1, 100), 0.5*linspace(0.0, 1, 100).^2,'-.r'); ylim([0.0 0.5]);
            plot(xplot, psilim, '--r');
            title('Figure of Merit (\phi vs \psi)');legend('M_f', 'Diffuser = 1', 'Symmetric Diffuser limit', 'Location', 'northwest');
            xlabel('Flow Coefficient (\phi)'); ylabel('Stage Loading (\psi)');
%             fignamefom = ['FoM_' num2str(phi) '-' num2str(psi) '_' num2str(sig) '.png'];
%             saveas(gcf,fignamefom);
        case 'phipsi_omega'
            figure(112); hold on; [C1,h1] = contourf(xplot, yplot, omegaCRITMESH*30/pi,'LineColor','none'); clabel(C1,h1); title('Omega \phi vs \psi');xlabel('Flow Coefficient (\phi)'); ylabel('Stage Loading (\psi)'); colorbar('Ticks',[8000 10000 12000 14000 16000 18000 20000]); ylim([0.0 0.5]);legend('Omega','Location', 'northwest');
            plot(linspace(0.0, 1, 100), 0.5*linspace(0.0, 1, 100).^2,'k', 'LineWidth', 2); 
            plot(linspace(0.0, 1, 100), 0.5*linspace(0.0, 1, 100).^2/1.44,'-.k', 'LineWidth', 2); 
%             plot(xplot, psilim, '-.k', 'LineWidth', 2); 
%             plot(xplot,powerlim,'-r', 'LineWidth', 2); 
            legend('Omega', 'Diffuser = 1', 'Diffuser length limit', 'Location', 'northwest');
            %             plot(linspace(0.0, 1, 100), 0.5*linspace(0.0, 1, 100).^2,'-.r');plot(xplot, psilim, '--r');
            
            figure(109); hold on; [Cp,hp] = contourf(xplot, yplot, OPpowerMESH,'LineColor','none'); title('Power \phi vs \psi'); xlabel('Flow Coefficient (\phi)'); ylabel('Stage Loading (\psi)'); colorbar(); ylim([0.0 0.5]);
            plot(linspace(0.0, 1, 100), 0.5*linspace(0.0, 1, 100).^2,'k', 'LineWidth', 2); 
            plot(linspace(0.0, 1, 100), 0.5*linspace(0.0, 1, 100).^2/1.44,'-.k', 'LineWidth', 2); 
%             plot(xplot, psilim, '-.k', 'LineWidth', 2); 
%             plot(xplot,powerlim,'-r', 'LineWidth', 2);
            legend('Power', 'Diffuser = 1', 'Diffuser length limit', 'Location', 'northwest');
            
            figure(111); hold on; contourf(xplot, yplot, OPmeritMESH, linspace(0,3,9)); colorbar();
            plot(linspace(0.0, 1, 100), 0.5*linspace(0.0, 1, 100).^2,'k', 'LineWidth', 2); 
            plot(linspace(0.0, 1, 100), 0.5*linspace(0.0, 1, 100).^2/1.44,'-.k', 'LineWidth', 2); 
            plot(linspace(0.0, 1, 100), 2*linspace(0.0, 1, 100).^2,'-.r', 'LineWidth', 2); 
%             plot(xplot, psilim, '-.k', 'LineWidth', 2); 
%             plot(xplot,powerlim,'-r', 'LineWidth', 2); 
            title('Figure of Merit (\phi vs \psi)');
            legend('M_f', 'Diffuser = 1', 'Diffuser length limit', 'Propeller Max', 'Location', 'northwest');
            xlabel('Flow Coefficient (\phi)'); ylabel('Stage Loading (\psi)'); ylim([0 0.5]);
%             fignamefom = ['FoM_' num2str(phi) '-' num2str(psi) '_' num2str(sig) '.png'];
%             saveas(gcf,fignamefom);
        case 'phisigma'
            figure(2); [C1,h1] = contourf(xplot, yplot, powerMESH); clabel(C1,h1,0); title('Power \phi vs \sig');xlabel('Flow Coefficient (\phi)'); ylabel('Diffusion Factor (\sig)'); colorbar();
            figure(4); [C1,h1] = contourf(xplot, yplot, thrustMESH); clabel(C1,h1,0); title('Thrust \phi vs \sig');xlabel('Flow Coefficient (\phi)'); ylabel('Diffusion Factor (\sig)'); colorbar();
            figure(5); [C1,h1] = contourf(xplot, yplot, massMESH); clabel(C1,h1,0); title('Delta Mass \phi vs \sig');xlabel('Flow Coefficient (\phi)'); ylabel('Diffusion Factor (\sig)'); colorbar();
            figure(3); [C2,h2] = contourf(xplot, yplot, meritMESH); colorbar();
            title('Figure of Merit (\phi vs \sig)');
            xlabel('Flow Coefficient (\phi)'); ylabel('Diffusion Factor (\sig)');

    end 
end

%% Clear all
clear *MESH
end
