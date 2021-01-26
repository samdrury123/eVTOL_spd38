function [R,T,Z,shiftvary] = stackblade(bezier_chi1,bezier_chi2,bezier_c,N_bezier,bladetype,rc,rh)
% Find B-spline interpolated blade by stacking many 2-D sections
% and convert the curved blade coordinates into Cartesians

global plotflag

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

save([pwd '/Geometry/' bladetype '.mat'],'r','chi1d','chi2d', 'c')
N = length(c);

% Loop over radial direction and create a new blade section at every go
count = 0;

% Stator Shape
theta1 = 25*pi/180;
theta2 = 10*pi/180;
ac = (tan(theta1)+tan(theta2))/(2*rh-2*rc);
bc = tan(theta1)-2*ac*rh;
cc = -ac*rh^2-bc*rh;

for i=1:N
    [xblade,yblade]=bladecentroid(chi1d(i),chi2d(i),c(i),bladetype);
    R(i,:)=r(i)*ones(1,length(xblade));
    
    if strcmpi(bladetype, 'stator')
        yshift = ac*(R(i,:).^2)+bc*R(i,:)+cc;
        yblade = yblade+yshift;
        shiftvary(i) = mean(yshift);
    else
        shiftvary(i) = 0;
    end
    
    T(i,:)=yblade/r(i);
    Z(i,:)=xblade;
     
    if plotflag == 1
    if (i == 1 || i == 10 || i == 21)
        figure(1); subplot(1,2,1); set(gcf, 'Position', [0 320 1500 1000]);
        count   = count + 1;
        cy      = (count-1) * +0.03;
        cx      = 0.02;
        if strcmp(bladetype,'rotor')
            cx      = -.02;
        end
        plot(xblade + cx,yblade + cy,'k-','LineWidth',1);
        hold on;
        xlim ([-0.025, 0.055]);
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