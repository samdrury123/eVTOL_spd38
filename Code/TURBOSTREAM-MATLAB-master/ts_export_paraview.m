function g = ts_export_paraview(g,filename,bp,bids,conv_cart,repeat_blade)
% TS_EXPORT_PARAVIEW  Write a pair of .hdf5 and .xdmf files for useful visualisation
%
%   g = TS_EXPORT_PARAVIEW(g,filename,bp,bids,conv_cart)
%
%   g - cell array of TURBOSTREAM solution
%   filename - string of .hdf5 filename to write to
%   bp - cell array of block properties to write or predefined set
%   bids - vector of block ids to include or leave empty for all
%   conv_cart - 0 or 1 to write cartesian coordinates and velocities
%   repeat_blade - 1 or 2 to specify blade repeats for single passage or full annulus

% Default to all blocks
if exist('bids','var') == 0 || isempty(bids) == 1
    bids = ts_get_bids(g);
end

% Get secondary coordinates if a flow solution exists
if isfield(g{1},'ro') == 1
    g = ts_secondary(g);
end

% Default to repeat the blades for tracking mode
if exist('repeat_blade','var') == 0 || isempty(repeat_blade) == 1
    repeat_blade = 1;
end


% Cut mesh into surfaces if specified
if strcmp(bp,'Surf') == 1
    g = ts_secondary(ts_get_surface(g,bids)); bp = []; bids = ts_get_bids(g);
end

% Choose defaults if not specified specifically in a cell array
if iscell(bp) == 0
    if strcmp(bp,'Primary') == 1
        bp = {'ro' 'rovx' 'rovr' 'rorvt' 'roe' 'trans_dyn_vis' 'misc'};
    end
    if strcmp(bp,'HighSpeed') == 1
         bp = {'mwall', 'phi', 'M_rel', 'sfunc', 'Alpha_rel', 'P'};
    end
    if strcmp(bp,'LowSpeed') == 1
         bp = {'mwall', 'phi', 'Po_rel', 'Alpha_rel'};
    end
    if strcmp(bp,'Track') == 1
        bp = {'mwall','ro','roe'}; 
%         bp = {'mwall','Po','P','Vn','Vs','Alpha'}; track_mode = 1;
        
        % Put in J and K indices to help recover streamlines later
        for bid = ts_get_bids(g)
            g{bid+1}.I = repmat(double(reshape(0:g{bid+1}.attribute.ni-1,[],1,1)),...
                [1 g{bid+1}.attribute.nj g{bid+1}.attribute.nk]);
            g{bid+1}.J = repmat(double(reshape(0:g{bid+1}.attribute.nj-1,1,[],1)),...
                [g{bid+1}.attribute.ni 1 g{bid+1}.attribute.nk]);
            g{bid+1}.K = repmat(double(reshape(0:g{bid+1}.attribute.nk-1,1,1,[])),...
                [g{bid+1}.attribute.ni g{bid+1}.attribute.nj 1]);
        end
        
        % Repeat blade row at +/- a pitch to enable proper tracking without
        % connectivity data
        if repeat_blade == 1
            g1 = g; g2 = g; 
    %         g3 = g; g4 = g; 
            for bid = ts_get_bids(g)
                nblade = double(g{bid+1}.bv.nblade);
                g1{bid+1}.rt = g1{bid+1}.rt - g1{bid+1}.r * 2 * pi / nblade;
                g2{bid+1}.rt = g2{bid+1}.rt + g2{bid+1}.r * 2 * pi / nblade;
    %             g3{bid+1}.rt = g3{bid+1}.rt - g3{bid+1}.r * 4 * pi / nblade;
    %             g4{bid+1}.rt = g4{bid+1}.rt + g4{bid+1}.r * 4 * pi / nblade;
            end
            g = [g ; g1 ; g2]; 
            bids = [bids bids + length(g2) bids + 2*length(g2)];
    %         g = [g ; g1 ; g2 ; g3 ; g4]; 
    %         bids = [bids bids + length(g2) bids + 2*length(g2) bids +
    %         3*length(g2) bids + 4*length(g2)];
        elseif repeat_blade == 2
            % Full annluus
            g1 = g; g2 = g; 
    %         g3 = g; g4 = g; 
            for bid = ts_get_bids(g)
                g1{bid+1}.rt = g1{bid+1}.rt - g1{bid+1}.r * 2 * pi;
                g2{bid+1}.rt = g2{bid+1}.rt + g2{bid+1}.r * 2 * pi;
    %             g3{bid+1}.rt = g3{bid+1}.rt - g3{bid+1}.r * 4 * pi / nblade;
    %             g4{bid+1}.rt = g4{bid+1}.rt + g4{bid+1}.r * 4 * pi / nblade;
            end
            g = [g ; g1 ; g2]; 
            bids = [bids bids + length(g2) bids + 2*length(g2)];
    %         g = [g ; g1 ; g2 ; g3 ; g4]; 
    %         bids = [bids bids + length(g2) bids + 2*length(g2) bids +
    %         3*length(g2) bids + 4*length(g2)];
        end
    end
end

% Determine cell array for hdf5 file
bp_hdf5_master = {'ro' 'fl' ; 'rovx' 'fl' ; 'rovr' 'fl' ; 'rorvt' 'fl' ; 'roe' 'fl';...
    'mwall' 'in8'; 'tlength' 'fl'; 'trans_dyn_vis' 'fl'; 'xlength' 'fl';  'misc' 'fl';...
    'dwallsq' 'fl' ; 'phi' 'fl' ; 'yplus' 'fl' ; 'egen' 'fl' ; 'tdamp' 'fl'};

bp_hdf5 = cell(length(bp),2);
for n = 1:length(bp)
    if sum(strcmp(bp(n),bp_hdf5_master(:,1))) > 0
        bp_hdf5{n,1} = bp{n};
        bp_hdf5{n,2} = bp_hdf5_master{strcmp(bp(n),bp_hdf5_master(:,1)),2};
    else
        bp_hdf5{n,1} = bp{n};
        bp_hdf5{n,2} = 'fl';
    end
end
bp_hdf5 = [{ 'x' 'fl' ; 'r' 'fl' ; 'rt' 'fl'} ; bp_hdf5];

% Determine cell array for xdmf file
bp_xdmf_master = {'ro' 'Float' 4; 'rovx' 'Float' 4; 'rovr' 'Float' 4; 'rorvt' 'Float' 4;...
    'roe' 'Float' 4; 'mwall' 'Char' 1; 'tlength' 'Float' 4; 'trans_dyn_vis' 'Float' 4;...
    'xlength' 'Float' 4;  'misc' 'Float' 4; 'dwallsq' 'Float' 4;...
    'phi' 'Float' 4; 'yplus' 'Float' 4; 'egen' 'Float' 4; 'tdamp' 'Float' 4};

bp_xdmf = cell(length(bp),3);
for n = 1:length(bp)
    if sum(strcmp(bp(n),bp_xdmf_master(:,1))) > 0
        bp_xdmf{n,1} = bp{n};
        bp_xdmf{n,2} = bp_xdmf_master{strcmp(bp(n),bp_xdmf_master(:,1)),2};
        bp_xdmf{n,3} = bp_xdmf_master{strcmp(bp(n),bp_xdmf_master(:,1)),3};
    else
        bp_xdmf{n,1} = bp{n};
        bp_xdmf{n,2} = 'Float';
        bp_xdmf{n,3} = 4;
    end
end

% Convert to cartesian coordinates and velocities if requested
if exist('conv_cart','var') ~= 0 && conv_cart == 1
    for bid = 1:length(g)
        g{bid}.t = g{bid}.rt ./ g{bid}.r;
        [g{bid}.y, g{bid}.z] = pol2cart(g{bid}.t, g{bid}.r);
        
        if isempty(bp) == 0
            g{bid}.Vr = g{bid}.rovr ./ g{bid}.ro; g{bid}.Vt = g{bid}.rorvt ./ (g{bid}.ro .* g{bid}.r);
            g{bid}.Vy = g{bid}.Vr .* cos(g{bid}.t) - g{bid}.Vt_rel .* sin(g{bid}.t);
            g{bid}.Vz = g{bid}.Vt_rel .* cos(g{bid}.t) + g{bid}.Vr .* sin(g{bid}.t);
        end
        
        g{bid}.rt = g{bid}.z; g{bid}.r = g{bid}.y;
        if isempty(bp) == 0
            g{bid}.Vt = g{bid}.Vz; g{bid}.Vr = g{bid}.Vy;
            g{bid}.Vt_rel = g{bid}.Vt;
        end
    end
end

% Scale velocities to unit size
% if track_mode == 1
%     for bid = 1:length(g)
%         if g{bid}.bv.rpm == 0
%             V_mag = ( g{bid}.Vx.^2 + g{bid}.Vr.^2 + g{bid}.Vt.^2 ).^0.5;
%         else
%             V_mag = ( g{bid}.Vx.^2 + g{bid}.Vr.^2 + g{bid}.Vt_rel.^2 ).^0.5;
%         end
%         g{bid}.Vx = g{bid}.Vx ./ V_mag; g{bid}.Vr = g{bid}.Vr ./ V_mag;
%         g{bid}.Vt = g{bid}.Vt ./ V_mag; g{bid}.Vt_rel = g{bid}.Vt_rel ./ V_mag;
%     end
% end

% Delete pre-existing file
if exist(filename,'file') ~= 0
    delete(filename);
end

% Write pair of .hdf5 and .xdmf files
ts_write_hdf5(g,filename,bp_hdf5,bids);
ts_write_xdmf(g,filename,bp_xdmf,bids);

end