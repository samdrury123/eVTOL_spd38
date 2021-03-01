function g = ts_load_balance(g,n_proc,plot_stuff)
% TS_LOAD_BALANCE  Use basic cumilative sum to assign a processor id to each block
%
%   g = TS_LOAD_BALANCE(g,n_proc)
%
%   g - is input TURBOSTREAM cell array
%   n_proc - is number of GPUs to be used (not nodes)
%   plot_stuff - optional arugument 0 or 1 to display working

% Default to no plotting
if exist('plot_stuff','var') == 0
    plot_stuff = 0;
end

% Extract the number of cells from the blocks
n_cells = zeros(size(g)); bids = ts_get_bids(g);
for bid = ts_get_bids(g)
    n_cells(bid+1) = numel(g{bid+1}.x);
end

% Target number of cells on each processor
n_total = sum(n_cells);
n_target = round(n_total / n_proc);

% Initial guess from cumilative sum of cells
q = cumsum(n_cells);
i = [-1 round(interp1(q,bids,linspace(n_target,n_total,n_proc)))];

% Apply processor ids to selected blocks
for n = 1:n_proc
    for bid = i(n)+1:i(n+1)
        g{bid+1}.attribute.procid = n - 1;
    end
end

% Sum the number of cells on each processor to check
cell_count = zeros(n_proc,1);
for bid = ts_get_bids(g)
    cell_count(g{bid+1}.attribute.procid+1) = cell_count(g{bid+1}.attribute.procid+1) + n_cells(bid+1);
end

% Display discrepancy from target value
if plot_stuff == 1
    
    % Average discrepancy
    d_avg = mean(abs(n_target - cell_count)) / n_target;
    
    % Maximum discrepancy
    d_max = max(cell_count - n_target) / n_target;
    
    % Display output
    disp(['Load balancing :   Maximum discrepancy = ' num2str(d_max * 100) ...
        ' %   Average discrepancy = ' num2str(d_avg * 100) ' %'])    
end

end
