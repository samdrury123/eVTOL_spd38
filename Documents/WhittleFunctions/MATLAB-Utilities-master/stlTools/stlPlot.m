function h = stlPlot(v, f, h, col, name, flipdir)
%STLPLOT is an easy way to plot an STL object
%V is the Nx3 array of vertices
%F is the Mx3 array of faces
%col is 1x3 vector rgb colour
%NAME is the name of the object, that will be displayed as a title

if exist('h','var') == 0 || isempty(h) == 1
    h = figure();
else
    figure(h);
end

if exist('col','var') == 0
    col = [0.8 0.8 1];
end

if exist('flipdir','var') == 0
    flipdir = 1;
end

% Flip y and z directions
if flipdir == 1
    v = v(:,[1 3 2]); v(:,2) = - v(:,2);
end

% Patch object
object.vertices = v;
object.faces = f;
if size(f,1) < 20000
    patch(object,'FaceColor',col,'FaceLighting','gouraud','AmbientStrength', 0.4);
else
    patch(object,'FaceColor',col,'EdgeColor','none','FaceLighting','gouraud','AmbientStrength', 0.4);
end

% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');
% view([-135 35]);
grid on;
if exist('name','var') ~= 0
    title(name);
end

% Add axes labels
if flipdir == 1
    xlabel('X / mm'); ylabel('Z / mm'); zlabel('Y / mm');
else
    xlabel('X / mm'); ylabel('Y / mm'); zlabel('Z / mm');
end

end