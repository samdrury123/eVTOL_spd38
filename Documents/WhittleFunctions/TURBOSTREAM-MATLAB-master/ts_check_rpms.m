function h = ts_check_rpms(g,h)
% Plot out all block surfaces to determine whether they are rotating or not

% Open figure window if it is not already specified
if exist('h','var') == 0 || isempty(h) == 1
    h = figure(); hold on; grid on; box on;
    xlabel('Axial Coordinate'); ylabel('Radial Coordinate');
else
    figure(h);
end

% Face names to plot
facenames = {'rpmi1' 'rpmi2' 'rpmj1' 'rpmj2'};

% Loop over all blocks
for bid = ts_get_bids(g)
    
    % Size of current block
    ni = g{bid+1}.attribute.ni; nj = g{bid+1}.attribute.nj;
    
    % Coordinates of each face
    i = {1,ni,1:ni,1:ni}; j = {1:nj,1:nj,1,nj}; k = 1;
    
    % Loop over all faces and plot colour coded lines
    for v = 1:length(facenames)
        
        % Vector of bids
        bids = ones(length(i{v}),length(j{v})) * bid;
        
        % Check whether the face is rotating or not
        if isfield(g{bid+1}.bv,facenames{v}) && g{bid+1}.bv.(facenames{v}) ~= 0
            plot3(g{bid+1}.x(i{v},j{v},k),g{bid+1}.r(i{v},j{v},k),bids,'k-')
        else
            plot3(g{bid+1}.x(i{v},j{v},k),g{bid+1}.r(i{v},j{v},k),bids,'r.')
        end
    end

end

% View from above and make axes equal after bids have been plotted in third dimension
view(2); axis auto; axis equal;

% Add legend
l(1) = plot(nan,nan,'k-'); l(2) = plot(nan,nan,'r.'); 
legend(l,'Rotating','Stationary','location','best')

end

