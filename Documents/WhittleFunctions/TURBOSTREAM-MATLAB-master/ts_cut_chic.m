function [inlet,outlet,pout,rpm] = ts_cut_chic(filename,coords,force_read)
% TS_CUT_CHIC  Take cut planes of various types at inlet and outlet of a machine

% Default to read cut file
if exist('force_read','var') == 0
    force_read = 0;
end

% Check if a cut file already exists
if exist(strrep(filename,'.hdf5','_chic.mat'),'file') == 0 || force_read == 1

    % Read mesh file
    g = ts_read_hdf5(filename);
    
    % Find delta pressure from log file if throttle patch is being used
    [directory,logname,~] = fileparts(filename);
    logname = ['log_' strrep(logname,'_avg','') '.txt'];
    [~,~,dpout] = ts_plot_conv([directory filesep logname],0);

    % Find inlet and exit blocks for the whole domain and take cuts as defaults
    for bid = ts_get_bids(g)
        for pid = ts_get_pids(g,bid)
            if g{bid+1}.patch{pid+1}.attribute.kind == 0
                bid_in = bid; i_in = g{bid+1}.patch{pid+1}.attribute.ist;
            end
            if g{bid+1}.patch{pid+1}.attribute.kind == 1
                bid_out = bid; i_out = g{bid+1}.patch{pid+1}.attribute.ist;
                pout = g{bid+1}.patch{pid+1}.pout + dpout;
            end   
        end
    end

    % Determine which type of cuts to use        
    if exist('coords','var') == 0 || isempty(coords) == 1

        % Take structured cuts close to the inlet and exit
        inlet = ts_secondary(ts_structured_cut(g,bid_in,i_in+1,i_in+1,1,'en',1,'en'));
        outlet = ts_secondary(ts_structured_cut(g,bid_out,i_out+1,i_out+1,1,'en',1,'en'));

    elseif iscell(coords) == 0 && size(coords,2) == 2

        % Extract structured cut information
        bid_in = coords(1,1); i_in = coords(1,2);
        bid_out = coords(2,1); i_out = coords(2,2);

        % Take structured cuts
        inlet = ts_secondary(ts_structured_cut(g,bid_in,i_in+1,i_in+1,1,'en',1,'en'));
        outlet = ts_secondary(ts_structured_cut(g,bid_out,i_out+1,i_out+1,1,'en',1,'en'));

    elseif iscell(coords) == 0 && size(coords,2) == 3

        % Extract and assemble cuts from multiple passage blocks
        for m = 1:size(coords,1)

            % Find blocks in the right region
            C = cell(0,1); x = coords(m,1); r = coords(m,2); i = coords(m,3);
            for bid = ts_get_bids(g)
                if g{bid+1}.x(end,1,1) > x && g{bid+1}.x(1,1,1) < x && g{bid+1}.r(1,end,1) > r
                    C{end+1} = ts_structured_cut(g,bid,i,i,1,'en',1,'en');
                end
            end

            % Join cuts together into one
            c = C{1}; varnames = fieldnames(c);
            for a = 2:length(C)
                for v = 1:length(varnames)
                    if size(c.(varnames{v}),2) > 1
                        c.(varnames{v}) = [c.(varnames{v}) C{a}.(varnames{v})(:,2:end)];
                    end
                end
            end

            % Record cut
            if m == 1; inlet = c; else; outlet = c; end;

        end

    else

        % Extract coordinate information
        bids1 = coords{1}.bids; xr1 = coords{1}.xr;
        bids2 = coords{2}.bids; xr2 = coords{2}.xr;

        % Take structured unstructured cuts
        inlet = ts_secondary(ts_structunstruct_cut(g,xr1,bids1,2,1));
        outlet = ts_secondary(ts_structunstruct_cut(g,xr2,bids2,2,1));

    end
    
    % Calculate machine rpm as maximum value from cut
    varnames = {'rpm' 'rpmi1' 'rpmi2' 'rpmj1' 'rpmj2' 'rpmk1' 'rpmk2'};
    rpm = 0;
    for bid = ts_get_bids(g)
        for v = 1:length(varnames)
            rpm = max(rpm,g{bid+1}.bv.(varnames{v}));
        end
    end

    % Save inlet and outlet cuts
    save(strrep(filename,'.hdf5','_chic.mat'),'inlet','outlet','pout','rpm')
    
else

    % Load the cuts from file
    load(strrep(filename,'.hdf5','_chic.mat'))
    
end
        
end