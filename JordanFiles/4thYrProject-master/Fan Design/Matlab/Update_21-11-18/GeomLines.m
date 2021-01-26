function Geom= GeomLines
Geom.make_hub_line      = @make_hub_line;
Geom.make_shroud_line    = @make_shroud_line;


function [totpoints, steps] = make_hub_line(r_hub,nose)
t = linspace (pi/2,pi,25);
r = linspace (r_hub,r_hub,25);
[x, y]=pol2cart(t,r);
lambda = [1.5 0 ; 0 1];
c = lambda*[x;y];
c=rot90(c',2);
c(:,2)=c(:,2)+0.0125-0.00556;
H=[ c(end,1)	0.01250000000
    c(end,1)	0.02481812502
    c(end,1)	0.03731812502
    c(end,1)	0.05775052900
    c(end,1)	0.07775052900];
if nose==1
    steps = [c;H; c(end,1)	0.09];
elseif nose==0
    steps = [r_hub 0.00694-.01; c(end,:);H; c(end,1)	0.09];
else
    disp('whoups')
    return
end
steps = fliplr(steps);
totpoints = length (steps(:,1));

function [totpoints, steps] = make_shroud_line(r_tip,nose)
H=[ r_tip	0.01277064013
    r_tip	0.02101840118
    r_tip	0.03750109976
    r_tip	0.05599234990
    r_tip	0.07775052900];
if nose==1
    steps = [r_tip -0.15 ; H ; r_tip	0.09];
elseif nose==0
    steps = [r_tip 0.00694-.01; r_tip 0.00694 ; H ; r_tip	0.09];
else
    disp('whoups')
    return
end
steps = fliplr(steps);
totpoints = length (steps(:,1));
