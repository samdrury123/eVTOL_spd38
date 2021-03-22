function g = ts_autogrid_reader(plot_3d_file,bc_file)
% Read Autogrid output and convert into matlab data array for TS use

fid = fopen(plot_3d_file,'r');

% Read block shapes
nb = fscanf(fid,'%i',1);
ijk = reshape(fscanf(fid,'%i',nb*3),[3 nb]).';

% Read and shape mesh geometry
g = cell(nb,1);
for bid = 1:nb
    ni = ijk(bid,1); nj = ijk(bid,2); nk = ijk(bid,3);
    nijk = ijk(bid,1)*ijk(bid,2)*ijk(bid,3);
    yzx = fscanf(fid,'%f',nijk*3);
    y = yzx(1:nijk);
    z = yzx(nijk+1:2*nijk);
    x = yzx(2*nijk+1:end);
    
    g{bid}.x = reshape(x,[ni nj nk]);
    g{bid}.y = reshape(y,[ni nj nk]);
    g{bid}.z = reshape(z,[ni nj nk]);
    
    % Attributes
    g{bid}.attribute.bid = bid-1;
    g{bid}.attribute.ni = ni;
    g{bid}.attribute.nj = nj;
    g{bid}.attribute.nk = nk;
end

fclose(fid);

% Switch i & k directions
for bid = ts_get_bids(g)
    g{bid+1}.x = permute(g{bid+1}.x,[3 2 1]);
    g{bid+1}.y = permute(g{bid+1}.y,[3 2 1]);
    g{bid+1}.z = permute(g{bid+1}.z,[3 2 1]);
    nk = g{bid+1}.attribute.ni; ni = g{bid+1}.attribute.nk;
    g{bid+1}.attribute.ni = ni; g{bid+1}.attribute.nk = nk;
end

% Convert to polar coordinates
for bid = ts_get_bids(g)
    g{bid+1}.r = hypot(g{bid+1}.y,g{bid+1}.z);
    t = atan2(g{bid+1}.y,g{bid+1}.z);
    g{bid+1}.rt = t.*g{bid+1}.r;
    
    g{bid+1} = rmfield(g{bid+1},'y'); g{bid+1} = rmfield(g{bid+1},'z');
end

% Read in patches and apply to cell array
fid = fopen(bc_file,'r');
L = textscan(fid,'%s',inf,'delimiter','\n'); L = L{1};
fclose(fid);

% Collect patch data from file for each block
match = regexp(L,'ROTATION ');
l = find(cellfun('isempty',match) == 0); l(end+1) = length(L);
P = cell(nb,6);
for bid = ts_get_bids(g)
    % Blade numbers
    nblade = regexp(L{l(bid+1)},'\s+','split');
    g{bid+1}.bv.nblade = str2double(nblade{8});
    g{bid+1}.bv.fblade = str2double(nblade{8});
    
    % Collect patch data
    lp = L(l(bid+1)+1:l(bid+2)-3);
    match = regexp(lp,'\d+'); ip = find(cellfun('length',match) == 1);
    np = str2double(lp(ip));
    
    pid = 0;
    for f = 1:6
        P{bid+1,f} = cell(np(f),1);
        for n = 1:np(f)
            temp = regexp(lp{ip(f)+n},'\s+','split');
            
            if sum(strcmp(temp{2},{'PER' 'CON' 'INL' 'OUT' 'NMB' 'ROT' 'PERNM' 'EXT'})) >= 1 
                P{bid+1,f}{n}.type = temp{2};
                P{bid+1,f}{n}.ind = str2double(temp(3:6));
                P{bid+1,f}{n}.pid = pid;
                P{bid+1,f}{n}.bid = bid;

                if sum(strcmp(temp{2},{'PER' 'CON' 'NMB' 'PERNM'})) >= 1 
                    P{bid+1,f}{n}.nxbid = temp{7}-1;
                    P{bid+1,f}{n}.nxp = temp{9}-1;
                end
                
                if strcmp(temp{2},'ROT') == 1
                    M = regexp(temp{1},')','split'); M = regexp(M{1},'(','split');
                    P{bid+1,f}{n}.id = str2double(M{2});
                end
                pid = pid+1;
            end
        end
    end
end

% Find ids, bids & pids of all mixing planes
ibp = [];
for f = 1:length(P)
    for n = 1:length(P{f})
        if strcmp(P{f}{n}.type,'ROT') == 1
            ibp = [ibp ; P{f}{n}.id P{f}{n}.bid P{f}{n}.pid];
        end
    end
end



"PER": ts_tstream_patch_kind.periodic,
"CON": ts_tstream_patch_kind.periodic,
"INL": ts_tstream_patch_kind.inlet,
"OUT": ts_tstream_patch_kind.outlet,
"NMB": ts_tstream_patch_kind.nomatch,
"ROT": ts_tstream_patch_kind.mixing,
"PERNM": ts_tstream_patch_kind.nomatch,
"EXT": ts_tstream_patch_kind.freestream}

% Flip coordinates in i-direction for left handed set
for bid = ts_get_bids(g)
    g{bid+1}.x = flipdim(g{bid+1}.x,1);
    g{bid+1}.r = flipdim(g{bid+1}.r,1);
    g{bid+1}.rt = flipdim(g{bid+1}.rt,1);

    % Change patch starts and ends in i direction
    ni = g{bid+1}.attribute.ni;
    for pid = ts_get_pids(g,bid)
        ist = g{bid+1}.patch{pid+1}.attribute.ist;
        ien = g{bid+1}.patch{pid+1}.attribute.ien;

        if ist ~= 0 || ien ~= ni
            g{bid+1}.patch{pid+1}.attribute.ist = ni - ien;
            g{bid+1}.patch{pid+1}.attribute.ien = ni - ist;
        end
    end
end

end