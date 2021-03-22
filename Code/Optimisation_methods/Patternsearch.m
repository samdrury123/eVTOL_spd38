clc; clear; close all;

x1 = linspace(-5.12,5.12,1000); % design space for parameter x1
x2 = linspace(-5.12,5.12,1000); % design space for paramter x2

[X1,X2] = meshgrid(x1,x2); % makes mesh of whole design space 

Z = zeros(size(X1)); % set cost function matrix of whole design space

for i = 1:numel(X1)
    Z(i) = rastriginfcn([X1(i),X2(i)]); % get cost function for cobination
end

fun = @rastriginfcn; % cost function handle

lb = [-100,-100]; % set lower bounds on x1 and x2 search space
ub = [100,100]; % set upper bounds on x1 and x2 search space
x0 = [-5.12,5.12]; % starting point

figure(1)
hold on
contour(X1,X2,Z) % plots function contours

tic
options = optimoptions('patternsearch','Plotfcn',...
    {@psplotbestf},'OutputFcn',{@plotfcn},'Display', 'iter',...
    'UseCompletePoll', true, 'MeshTolerance', 1e-8);
[x,fval,exitflag,output] = patternsearch(fun,x0,...
    [],[],[],[],lb,ub,options);
toc

figure(3)
mesh(X1,X2,Z)