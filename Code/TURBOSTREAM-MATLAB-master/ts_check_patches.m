function [] = ts_check_patches(g,plot_patches,filename)
% TS_CHECK_PATCHS  Print patch data for a given data cell array

% Default to plot bar charts but not individual patches
if exist('plot_patches','var') == 0
    plot_patches = 2;
end

% Default to print output to screen
if exist('filename','var') == 0
    fid = 1;    
    write_file = 0;
else
    fid = fopen(filename, 'w');
    write_file = 1;
end

% Number of blocks
nB = length(g);

% Preallocate matching arrays
lengths = zeros(nB,1);
corners = zeros(nB,1);
sxs = zeros(nB,1);
srs = zeros(nB,1);
sts = zeros(nB,1);

% Loop over all blocks
for n = 1:nB
    
    % Get current block data
    b0 = g{n};
    nP = b0.attribute.np;
    %nP = size(b0.patch,1);
    fprintf(fid,'\nBlock %i \t np %i \n\n', n-1, nP);
    
    % Loop over all patches
    for i = 1:nP
        
        % Patch attributes and sizes
        p0 = b0.patch{i}.attribute;
        is = [p0.ist p0.ien p0.jst p0.jen p0.kst p0.ken];
        
        % Check patch sizes are correct
        if p0.ist >= p0.ien
            error('\t Patch %i I-size mismatch\n',i-1);
        end
        if p0.jst >= p0.jen
            error('\t Patch %i J-size mismatch\n',i-1);
        end
        if p0.jst >= p0.jen
            error('\t Patch %i K-size mismatch\n',i-1);
        end        
        
        % Check all patch kinds
        if p0.kind == 0 || p0.kind == 1 || p0.kind == 19
            
            % Print inlet or exit patch data
            if p0.kind == 0
                fprintf(fid,'\t Patch %i (Inlet)\n',i-1);
                fprintf(fid,'\t ist = %i, ien = %i; jst = %i, jen = %i; kst = %i, ken = %i\n\n',is);
            elseif p0.kind == 1
                fprintf(fid,'\t Patch %i (Exit)\n',i-1);
                fprintf(fid,'\t ist = %i, ien = %i; jst = %i, jen = %i; kst = %i, ken = %i\n\n',is);
            elseif p0.kind == 19
                fprintf(fid,'\t Patch %i (2D Exit)\n',i-1);
                fprintf(fid,'\t ist = %i, ien = %i; jst = %i, jen = %i; kst = %i, ken = %i\n\n',is);
            end

        elseif p0.kind == 2 
            
            % Print mixing plane data
            b1 = g{p0.nxbid+1};   
            p1 = b1.patch{p0.nxpid+1}.attribute;
            fprintf(fid,'\t Patch %i (Mixing plane)\n',i-1);
            fprintf(fid,'\t ist = %i, ien = %i; jst = %i, jen = %i; kst = %i, ken = %i\n',is);
            fprintf(fid,'\t Next block = %i, next patch = %i\n\n',p1.bid,p1.pid);
            
        elseif p0.kind == 3 
            
            % Print sliding plane data
            fprintf(fid,'\t Patch %i (Sliding plane)\n',i-1);
            fprintf(fid,'\t ist = %i, ien = %i; jst = %i, jen = %i; kst = %i, ken = %i\n',is);
            fprintf(fid,'\t Next block = %i, next patch = %i\n\n',p0.nxbid,p0.nxpid);
            
        elseif p0.kind == 7 
            
            % Print slipwall data
            fprintf(fid,'\t Patch %i (Slipwall)\n',i-1);
            fprintf(fid,'\t ist = %i, ien = %i; jst = %i, jen = %i; kst = %i, ken = %i\n\n',is);
            
        elseif p0.kind == 10 
            
            % Print bleed patch data
            fprintf(fid,'\t Patch %i (Bleed)\n',i-1);
            fprintf(fid,'\t ist = %i, ien = %i; jst = %i, jen = %i; kst = %i, ken = %i\n\n',is);
        
        elseif p0.kind == 15 
            
            % Print non-matching patch data
            b1 = g{p0.nxbid+1};   
            p1 = b1.patch{p0.nxpid+1}.attribute;
            fprintf(fid,'\t Patch %i (No Match)\n',i-1);
            fprintf(fid,'\t ist = %i, ien = %i; jst = %i, jen = %i; kst = %i, ken = %i\n',is);    
            fprintf(fid,'\t Next block = %i, next patch = %i\n\n',p1.bid,p1.pid);
            
        elseif p0.kind == 5 || p0.kind == 16 
            
            % Check periodic or cartesian periodic indices fit on the block
            if p0.ien <= b0.attribute.ni && p0.jen <= b0.attribute.nj && p0.ken <= b0.attribute.nk
                
                % Indices of current patch
                I = p0.ist+1:p0.ien; J = p0.jst+1:p0.jen; K = p0.kst+1:p0.ken;
                
                % Get the current patch coordinates
                p0.x = squeeze(b0.x(I,J,K)); p0.r = squeeze(b0.r(I,J,K));
                p0.rt = squeeze(b0.rt(I,J,K)); p0.t = p0.rt./p0.r;
                p0.ni = p0.ien - p0.ist; p0.nj = p0.jen - p0.jst; p0.nk = p0.ken - p0.kst;
                
                % Get the matching block
                b1 = g{p0.nxbid+1};
                
                % Check the connecting patch exists
                if p0.nxpid+1 <= length(b1.patch)

                    % Get coordinates of the matching patch
                    p1 = b1.patch{p0.nxpid+1}.attribute;
                    p1 = ts_get_match(p0,p1,b1);
                    p1.ni = p1.ien - p1.ist; p1.nj = p1.jen - p1.jst; p1.nk = p1.ken - p1.kst;

                    % Check lengths of the patches match up
                    l1 = sort([p0.ni p0.nj p0.nk]);
                    l2 = sort([p1.ni p1.nj p1.nk]);
                    delta = l2 - l1;
                    
                    % Correct tangential coordinates by a whole number of pitches
                    t0 = mean(p0.t(:)); t1 = mean(p1.t(:)); 
                    dt = 2 * pi / double(g{n}.bv.nblade);
                    dt = round((t1 - t0) / dt) * dt;
                    p1.t = p1.t - dt;

                    % Check differences in coordinates 
                    if numel(p0.x) == numel(p1.x) && size(p1.x,1) == size(p0.x,1)
                        sdx = sum(abs(p0.x(:) - p1.x(:)));
                        sdr = sum(abs(p0.r(:) - p1.r(:)));
                        sdt = sum(abs(p0.t(:) - p1.t(:)));
                        max_dx = max(abs(p0.x(:) - p1.x(:)));
                        max_dr = max(abs(p0.r(:) - p1.r(:)));
                        max_dt = max(abs(p0.t(:) - p1.t(:)));
                    else
                        sdx = 1; sdr = 1; sdt = 1;
                        max_dx = 1; max_dr = 1; max_dt = 1;
                    end

                    % Check 4 corner points of the patch in order to check directions
                    if p0.bid ~= p0.nxbid
                        Q = ts_check_dir(p0,p1);
                    end

                    % Print the kind of patch
                    if p0.kind == 5
                        fprintf(fid,'\t Patch %i (Periodic)\n',i-1);
                    elseif p0.kind == 16
                        fprintf(fid,'\t Patch %i (Periodic Cartesian)\n',i-1);
                    end

                    % Print matching data
                    fprintf(fid,'\t Next block = %i, next patch = %i\n',p0.nxbid,p0.nxpid);
                    fprintf(fid,'\t ist = %i, ien = %i; jst = %i, jen = %i; kst = %i, ken = %i\n',is);
                    fprintf(fid,'\t idir = %i; jdir = %i; kdir = %i\n',p0.idir,p0.jdir,p0.kdir);
                    if p0.bid ~= p0.nxbid
                        fprintf(fid,'\t Corner distances = %8.6f, %8.6f, %8.6f, %8.6f\n',...
                            Q(1),Q(2),Q(3),Q(4));
                    end
                    fprintf(fid,'\t Length difference = %i, %i, %i\n',delta(1), delta(2), delta(3));
                    fprintf(fid,'\t sdx = %12.8f; sdr = %12.8f; sdt = %12.8f\n',sdx,sdr,sdt);
                    fprintf(fid,'\t max_dx = %0.5g; max_dr = %0.5g; max_dt = %0.5g\n\n',...
                        max_dx,max_dr,max_dt);

                    % Save data for plotting
                    sxs(n) = sxs(n) + sdx;
                    srs(n) = srs(n) + sdr;
                    sts(n) = sts(n) + sdt;            

                    % Analyse corner distances
                    lengths(n) = lengths(n) + delta(1) + delta(2) + delta(3);
                    if p0.bid ~= p0.nxbid
                        corners(n) = corners(n) + Q(1) + Q(2) + Q(3) + Q(4);
                    end

                    % Plot meshes of the patch if requested
                    if plot_patches == 1
                        figure(); hold on;
                        mesh(double(p0.x),double(p0.r),double(p0.r).*double(p0.t))
                        mesh(double(p1.x),double(p1.r),double(p1.r).*double(p1.t))
                        axis equal; view(3);
                    end            
                else
                    fprintf(fid,'\t Patch Mismatch\n\n');
                end
            else
                fprintf(fid,'\t Patch %i Block Dimensions Exceeded\n',i-1);
            end
        end
    end
end

% Warn about repeated patches
for bid = ts_get_bids(g);
    pids = ts_get_pids(g{bid+1});
    for pid = pids
        if sum(pid==pids) > 1
            fprintf(fid,'\t Patch %i has been repeated on block %i\n',pid,bid);
        end
    end
end

% Close the output file
if write_file
    fclose(fid);
end

% Plot mismatches as bar charts to see problem patches
if plot_patches == 2
    
    % Open figure window
    figure(); bids = ts_get_bids(g);
    
    % Mismatch in coordinates
    subplot(2,2,1); bar(bids,sxs); ylabel('\Sigmadx');
    subplot(2,2,2); bar(bids,srs); ylabel('\Sigmadr');
    subplot(2,2,3); bar(bids,sts); ylabel('\Sigmadt');

    % Mismatch in indices
    subplot(2,2,4); bar(bids,lengths); ylabel('Lengths');
end

end

            