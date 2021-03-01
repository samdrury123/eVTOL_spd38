
% Clear the workspace
close all; clear; clc;
% 
% %Load blade script example
% Bl_Example

close all

% Choose what to make splitter blade for
diffuser_splitter = 0;
inlet_splitter = 0;
diffuser_short_splitter =1;


% If diffuser splitter choose angle of diffuser to design splitter for,
% update file name to save as in script
diffuser_angle = 19.425;
r1 = 87.5;

% r2 based off of equal area, if not choose r2
equal_area = 1;

if equal_area == 0
    r2 = 80;
end 

%% Diffuser splitter

if diffuser_splitter == 1

    
% Calculate r2 of splitter blade for equal area ratio

if equal_area == 1 
    A1a = pi*(120.^2 -r1.^2);
    A1b = pi*(r1.^2 -55.^2);
    num = A1b*(120+130*tan(diffuser_angle*2*pi/360)).^2 + A1a*(55-130*tan(diffuser_angle*2*pi/360)).^2;
    r2 = (num/(A1a+A1b)).^0.5
end

    gtg
% Construct blade geometry
% c = bl_spline_eval(b,0.5); ote = 0;
% % c.s_thick_max = 0.2; 
% c.thick_te = 0.7;
% c.rad_thick_max = 3;
% c.rad_le = 30;
load('thickness_1.mat');m = bl_construct_thickness(c,102,0,1);


% Plot diffuser geometry for equal area 
figure
hold on
line([110,125],[120,120])
line([110 125], [125, 125])
line([125,255],[120, (120+130*tan(diffuser_angle*2*pi/360))])
line([125,255],[125, (125+130*tan(diffuser_angle*2*pi/360))])
line([110,125],[55,55])
line([125,255],[55,(55-130*tan(diffuser_angle*2*pi/360))])

% Plot line for splitter chord
x = linspace(125, 255, 301);
r = linspace(r1, r2, 301);
plot(x,r)

% Add thickness to splitter chord
deltax = x(2:end)-x(1:end-1);
deltar = r(2:end)-r(1:end-1);
p(1,:) = deltax;
p(2,:) = deltar;
clear n
n(1,:)=-deltar;
n(2,:)=deltax;

thick(1) = 0;
thick(2:27) = m.thick(15:40);
thick(27:240)= m.thick(40);
thick(240:299)= m.thick(40:99);
thick(300:300) = 0;

x = linspace(125, 255, 300);
r = linspace(r1, r2, 300);
clear xr
xr(1,:)= x;
xr(2,:)= r;

for i = 1:length(thick)

n(:,i) = n(:,i)/norm(n(:,i));

xr1(:,i) = xr(:,i) + n(:,i)*thick(i);
xr2(:,i) = xr(:,i) - n(:,i)*thick(i);

end

% Plot splitter
plot(xr1(1,:), xr1(2,:));
plot(xr2(1,:), xr2(2,:));
xlim([0 270])
ylim([0 170])

XR = [xr1 flip(xr2,2)];


% Save as .ibl file
statorcenterlinefile = fullfile(pwd,'/15_degree_splitter_3.ibl');

fid = fopen(statorcenterlinefile, 'w');
fprintf(fid,'Closed Index arclength\n');
fprintf(fid,'begin section\n');
fprintf(fid,'begin curve\n');
for ss = 1:size(XR,2)
     fprintf(fid,'%f  %f  %f\n',0,XR(2,ss),XR(1,ss));
end
fclose(fid);


hold off


% For nurbs method
% z1=linspace(-0.001,-0.001,300);
% z2=linspace(0.001, 0.001, 300);
% 
% Z = [xr1(1,:)' xr1(1,:)']
% Y = [xr1(2,:)' xr1(2,:)']
% X = [z1' z2']
% xyz(:,:,1) = X
% xyz(:,:,2) = Y
% xyz(:,:,3) = Z
% srf = nrbloft(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),2);
% igesout(srf, 'split.iges');


end 


%% Diffuser splitter

if diffuser_short_splitter == 1

    
% Calculate r2 of splitter blade for equal area ratio

if equal_area == 1 
    A1a = pi*(120.^2 -r1.^2);
    A1b = pi*(r1.^2 -55.^2);
    num = A1b*(120+65*tan(diffuser_angle*2*pi/360)).^2 + A1a*(55-65*tan(diffuser_angle*2*pi/360)).^2;
    r2 = (num/(A1a+A1b)).^0.5
end

    
% Construct blade geometry
% c = bl_spline_eval(b,0.5); ote = 0;
% % c.s_thick_max = 0.2; 
% c.thick_te = 0.7;
% c.rad_thick_max = 3;
% c.rad_le = 30;
load('thickness_1.mat')
m = bl_construct_thickness(c,102,0,1);


% Plot diffuser geometry for equal area 
figure
hold on
line([110,125],[120,120])
line([110 125], [125, 125])
line([125,190],[120, (120+65*tan(diffuser_angle*2*pi/360))])
line([125,190],[125, (125+65*tan(diffuser_angle*2*pi/360))])
line([110,125],[55,55])
line([125,190],[55,(55-65*tan(diffuser_angle*2*pi/360))])

% Plot line for splitter chord
x = linspace(125, 190, 301);
r = linspace(r1, r2, 301);
plot(x,r)

% Add thickness to splitter chord
deltax = x(2:end)-x(1:end-1);
deltar = r(2:end)-r(1:end-1);
p(1,:) = deltax;
p(2,:) = deltar;
clear n
n(1,:)=-deltar;
n(2,:)=deltax;

thick(1) = 0;
thick(2:27) = m.thick(15:40);
thick(27:240)= m.thick(40);
thick(240:299)= m.thick(40:99);
thick(300:300) = 0;

x = linspace(125, 190, 300);
r = linspace(r1, r2, 300);
clear xr
xr(1,:)= x;
xr(2,:)= r;

for i = 1:length(thick)

n(:,i) = n(:,i)/norm(n(:,i));

xr1(:,i) = xr(:,i) + n(:,i)*thick(i);
xr2(:,i) = xr(:,i) - n(:,i)*thick(i);

end

% Plot splitter
plot(xr1(1,:), xr1(2,:));
plot(xr2(1,:), xr2(2,:));
xlim([0 270])
ylim([0 170])

XR = [xr1 flip(xr2,2)];


% Save as .ibl file
statorcenterlinefile = fullfile(pwd,'/61_short_splitter.ibl');

fid = fopen(statorcenterlinefile, 'w');
fprintf(fid,'Closed Index arclength\n');
fprintf(fid,'begin section\n');
fprintf(fid,'begin curve\n');
for ss = 1:size(XR,2)
     fprintf(fid,'%f  %f  %f\n',0,XR(2,ss),XR(1,ss));
end
fclose(fid);


hold off


% For nurbs method
% z1=linspace(-0.001,-0.001,300);
% z2=linspace(0.001, 0.001, 300);
% 
% Z = [xr1(1,:)' xr1(1,:)']
% Y = [xr1(2,:)' xr1(2,:)']
% X = [z1' z2']
% xyz(:,:,1) = X
% xyz(:,:,2) = Y
% xyz(:,:,3) = Z
% srf = nrbloft(xyz(:,:,1),xyz(:,:,2),xyz(:,:,3),2);
% igesout(srf, 'split.iges');


end 


%% Inlet splitter vane 


if inlet_splitter == 1
    
% Construct blade geometry
% c = bl_spline_eval(b,0.5); ote = 0;
% c.s_thick_max = 0.4;
% c.thick_te = 0.6;
% c.rad_thick_max = 2;
load('thickness_1.mat')
m = bl_construct_thickness(c,102,0,1);
  
% Construct and plot ellipse for inlet casing
figure; hold on
r1 = 55; % x length of ellipse
r2 = 35; %y length of ellipse
theta = 0; %angle you want to rotate the ellipse to relative to x axis
theta = theta * 2*pi/360;
xc = 0;  %offset the ellipse from axis in the x direction
yc = -(120+r2); %offset the ellipse from axis in the y direction

t = linspace(pi/2, pi, 301);
xt = (r1 * cos(t)) + xc;
yt = (r2 * sin(t)) + yc;
% aply rotation by angle theta
cot = cos(theta); sit = sin(theta);
x = xt * cot - yt * sit;
y = xt * sit - yt * cot;
% draw the curve
plot(x, y, '-');


% Construct and plot ellipse for splitter chord
r1 = 45; % x length of ellipse
r2 = 35; %y length of ellipse
theta = 0; %angle you want to rotate the ellipse to relative to x axis
theta = theta * 2*pi/360;
xc = -25;  %offset the ellipse from axis in the x direction
yc = -(100+r2); %offset the ellipse from axis in the y direction

t = linspace(pi/2, 2.7*pi/3, 301);
xt = (r1 * cos(t)) + xc;
yt = (r2 * sin(t)) + yc;
% aply rotation by angle theta
cot = cos(theta); sit = sin(theta);
x = xt * cot - yt * sit;
r = xt * sit - yt * cot;
% draw the curve
plot(x, r, '-');


% Add thickness to splitter chord
deltax = x(2:end)-x(1:end-1);
deltar = r(2:end)-r(1:end-1);
p(1,:) = deltax;
p(2,:) = deltar;
clear n
n(1,:)=-deltar;
n(2,:)=deltax;

thick(1) = 0 ;
thick(2:29) = flip(m.thick(72:99));
thick(29:270)= m.thick(22);
thick(270:299)= m.thick(70:99);
thick(300:300) = 0;

t = linspace(pi/2, 2.7*pi/3, 300);
xt = (r1 * cos(t)) + xc;
yt = (r2 * sin(t)) + yc;
% apply rotation by angle theta
cot = cos(theta); sit = sin(theta);
x = xt * cot - yt * sit;
r = xt * sit - yt * cot;
clear xr
xr(1,:)= x;
xr(2,:)= r;

for i = 1:length(thick)

n(:,i) = n(:,i)/norm(n(:,i));

xr1(:,i) = xr(:,i) + n(:,i)*thick(i)*1.5;
xr2(:,i) = xr(:,i) - n(:,i)*thick(i)*1.5;

end

% Plot splitter
plot(xr1(1,:), xr1(2,:));
plot(xr2(1,:), xr2(2,:));
xlim([-150 0])
ylim([0 170])

XR = [xr1 flip(xr2,2)];

% Save as .ibl file
statorcenterlinefile = fullfile(pwd,'/inlet_45_splitter.ibl');

fid = fopen(statorcenterlinefile, 'w');
fprintf(fid,'Closed Index arclength\n');
fprintf(fid,'begin section\n');
fprintf(fid,'begin curve\n');
for ss = 1:size(XR,2)
     fprintf(fid,'%f  %f  %f\n',0,XR(2,ss),XR(1,ss));
end
fclose(fid);

end


