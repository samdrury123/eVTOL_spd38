function [] = ts_check_normals(g,normal)
% Moves and plots surface in direction of normal by d 
d = 0.0001;

% Plot original surface in black
ts_plot_surface(g, [0 0 0], 1, 0)

% Plot perturbed surface in red
for bid = ts_get_bids(g)
    g{bid+1}.x = g{bid+1}.x + d * normal{bid+1}(:,:,:,1);
    g{bid+1}.r = g{bid+1}.r + d * normal{bid+1}(:,:,:,2);
    g{bid+1}.rt = g{bid+1}.rt + d * normal{bid+1}(:,:,:,3);
end

ts_plot_surface(g, [1 0 0], 1, 0)

end