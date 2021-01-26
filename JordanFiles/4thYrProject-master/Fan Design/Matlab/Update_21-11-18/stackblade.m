function [R,T,Z] = stackblade(bezier_chi1,bezier_chi2,bezier_c,N_bezier,bladetype,file_name)
% Find B-spline interpolated blade by stacking many 2-D sections
% and convert the curved blade coordinates into Cartesians
global flagplot ac bc cc

r     = [];
chi1d = [];
chi2d = [];
c     = [];
n     = size(bezier_chi1.A1,1);
for i = 1:N_bezier
    r                = [r    ; bezier_chi1.(['A' num2str(i)])(:,2)];
    chi1d            = [chi1d; bezier_chi1.(['A' num2str(i)])(:,1)];
    chi2d            = [chi2d; bezier_chi2.(['A' num2str(i)])(:,1)];
    c                = [c    ; bezier_c.(['A' num2str(i)])(:,1)];
    
    if i<N_bezier
        r(1+(n-1)*i)     = [];
        chi1d(1+(n-1)*i) = [];
        chi2d(1+(n-1)*i) = [];
        c(1+(n-1)*i)     = [];
    end
    
end

save([file_name '_' bladetype '.mat'],'r','chi1d','chi2d', 'c')
N = length(c);

% Loop over radial direction and create a new blade section at every go
count = 0;

for i=1:N
    [xblade,yblade]=bladecentroid(chi1d(i),chi2d(i),c(i),bladetype);
    R(i,:)=r(i)*ones(1,length(xblade));
    
    
%     size(R(i,:))
%     R(i,:)
%     size(yblade)
%     yblade = yblade+0.01;
%     xblade
    
    if strcmpi(bladetype, 'stator')
        yshift = ac*(R(i,:).^2)+bc*R(i,:)+cc;
        yblade = yblade+yshift;
    end
    
    T(i,:)=yblade/r(i);
    Z(i,:)=xblade;
     
    if flagplot == 1
    if (i == 1 | i == 10 | i == 21)
        figure(5);
        count   = count + 1;
        cy      = (count-1) * +0.03;
        cx      = 0.02;
        if strcmp(bladetype,'rotor')
            cx      = -0.05;
        end
        plot(xblade + cx,yblade + cy,'k-','LineWidth',1);
        hold on;
        axis ([-0.06, 0.06, -0.02, 0.1]);
        axis equal;
        xlabel('x');
        ylabel('y');
    end
    end
end
% Take transpose of R,T,Z matrices
R=R';
T=T';
Z=Z';
end