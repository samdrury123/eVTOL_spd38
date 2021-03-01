function [] = printDiff(r_cCRIT, r_e, r_h, L)
% printDiff 
%   Get geometry of exit.

%% Setup
% Length of turbo machinery section
turboL = 0.01;

% x-index
xturbo = linspace(0,turboL,1000);
xend = turboL + L;
xdiff = linspace(turboL,xend,1000);

% Equation of surfaces
y1 = r_cCRIT + ((r_e - r_cCRIT)/L)*(xdiff-turboL);
y2 = r_h-(r_h/L)*(xdiff-turboL);

%% Plot
% Plot cone and casing
figure; hold on; plot(xdiff, y1, 'b'); plot(xdiff, y2, 'b'); 
plot(xdiff, -y1, 'b'); plot(xdiff, -y2, 'b');

% Plot turbo constant radius section
plot(xturbo, linspace(r_cCRIT,r_cCRIT,1000), 'r'); plot(xturbo, -linspace(r_cCRIT,r_cCRIT,1000), 'r');
plot(xturbo, linspace(r_h,r_h,1000), 'r'); plot(xturbo, -linspace(r_h,r_h,1000), 'r'); pbaspect([1 1 1]);
end