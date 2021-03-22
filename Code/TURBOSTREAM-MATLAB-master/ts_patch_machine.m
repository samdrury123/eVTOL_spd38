function [g] = ts_patch_machine(g, blade_rows, seals, rotating)
% Patches a machine correctly in the TS data format

% PATCH KINDS
% inlet = 0
% outlet = 1
% mixing = 2
% sliding = 3
% fixed = 4
% periodic = 5
% cooling = 6
% slipwall = 7
% probe = 8
% laminar = 9
% bleed = 10
% cusp = 11
% teforce = 12
% freestream = 13
% bladeforce = 14
% nomatch = 15
% periodic_cartesian = 16
% porous = 17
% isothermal = 18
% outlet2d = 19

% Build kd-tree
[face, face_type] = ts_get_faces(g);
[Tree T] = ts_build_tree(face, g, 0, 1);

% Remove all old patches & clists
nc = 0;
for bid = ts_get_bids(g)
    g{bid+1}.patch = cell(0,1);
%     for clid = ts_get_clids(bid,g)
%         nc = nc + g{bid+1}.clist{clid+1}.attribute.n;
%     end
%     g{bid+1}.clist = cell(0,1);
    if isfield(g{bid+1},'clist')
        g{bid+1} = rmfield(g{bid+1},'clist');
    end
    g{bid+1}.attribute.ncl = 0;
end

%% Inlet
disp('%% Inlet Patch %%')
p.attribute.kind = 0;
p.attribute.bid = blade_rows{1}.inlet;
p.attribute.pid = -1;

p.attribute.ist = 0;
p.attribute.ien = 1;
p.attribute.jst = 0;
p.attribute.jen = g{blade_rows{1}.inlet+1}.attribute.nj;
p.attribute.kst = 0;
p.attribute.ken = g{blade_rows{1}.inlet+1}.attribute.nk;

p.pitch = [];
p.pstag = [];
p.rfin = 0.5;
p.sfinlet = 0;
p.tstag = [];
p.yaw = [];

[ g ] = ts_add_patch(blade_rows{1}.inlet, p, g);
clear p

[g BC] = ts_apply_inlet(g, '/home/jvt24/Documents/TUROBSTREAM/SMURF-TS/BCs/p1b/s2.inflow.mp.dat', 0, -800.0005);

%% Mixing planes
disp('%% Mixing Planes %%')
for row = 2:length(blade_rows)
   
    bid0 = blade_rows{row-1}.exit;
    bid1 = blade_rows{row}.inlet;
    
    Q0 = [g{bid0+1}.attribute.ni-1 g{bid0+1}.attribute.ni ; 0 g{bid0+1}.attribute.nj ;...
        0 g{bid0+1}.attribute.nk];
    dir0 = [0 1 2];
    Q1 = [0 1 ; 0 g{bid1+1}.attribute.nj ; 0 g{bid1+1}.attribute.nk];
    dir1 = [0 1 2];
    g = ts_add_periodic_patch(bid0, bid1, Q0, Q1, dir0, dir1, g, 2);
    
end


%% Exit
disp('%% Exit Patch %%')
p.attribute.kind = 1;
bid = blade_rows{end}.exit;
p.attribute.bid = bid;
p.attribute.pid = -1;

p.attribute.ist = g{bid+1}.attribute.ni - 1;
p.attribute.ien = g{bid+1}.attribute.ni;
p.attribute.jst = 0;
p.attribute.jen = g{bid+1}.attribute.nj;
p.attribute.kst = 0;
p.attribute.ken = g{bid+1}.attribute.nk;

p.ipout = 3;
p.pout = 100000;
p.throttle_type = 0;

[ g ] = ts_add_patch(bid, p, g);
clear p

%% Periodic surfaces
if exist('Tree','var') == 0
    [face, face_type] = ts_get_faces(g);
    [Tree T] = ts_build_tree(face, g, 0, 1);
end

for n = 1:length(blade_rows)
    g = ts_patch_row(g,Tree,T,blade_rows{n}); 
end

if ~isempty(seals)
    for n = 1:length(seals)
        [g] = ts_patch_seal(g,Tree,T,seals{n});
    end
end

%% Rotating Surfaces
if isempty(rotating) ~= 1
    g = ts_set_rpms(g,rotating);
end

%% Set Patch Numbers
for bid = ts_get_bids(g)
    g{bid+1}.attribute.np = length(ts_get_pids(g{bid+1}));
end


























% Connectivity finder unfinished and incorrect

% if find_connectivity == 1
% 
%     % Find data on face_arrays
%     if exist('face','var') == 0
%         disp('%% Finding Faces %%')
%         [face, face_type] = ts_get_faces(g);
%     end
% 
%     % Build kd-tree
%     if exist('Tree','var') == 0
%         disp('%% Building KD-Tree %%')
%         [Tree T] = ts_build_tree(face, g, 0, 1);
%     end
% 
%     % Initialise c-lists
%     disp('%% Initialising C-Lists %%')
%     bids = ts_get_bids(g);
%     for bid = bids
%         g{bid+1}.clist = cell(1,length(bids));
%         for nxbid = bids
%             g{bid+1}.clist{nxbid+1}.attribute.bid = bid;
%             g{bid+1}.clist{nxbid+1}.attribute.clid = nxbid;
%             g{bid+1}.clist{nxbid+1}.attribute.nxbid = nxbid;
%             g{bid+1}.clist{nxbid+1}.attribute.kind = 0;
%             g{bid+1}.clist{nxbid+1}.attribute.n = 0;
% 
%             g{bid+1}.clist{nxbid+1}.ilist = zeros(0,1);
%             g{bid+1}.clist{nxbid+1}.jlist = zeros(0,1);
%             g{bid+1}.clist{nxbid+1}.klist = zeros(0,1);
%             g{bid+1}.clist{nxbid+1}.nxilist = zeros(0,1);
%             g{bid+1}.clist{nxbid+1}.nxjlist = zeros(0,1);
%             g{bid+1}.clist{nxbid+1}.nxklist = zeros(0,1);
%         end
%     end
% 
%     % Find connectivity data and write to clists
%     disp('%% Finding Connectivity %%')
% 
%     n = 1;
%     for bid = ts_get_bids(g)
%         for s = 1:size(face,2)
%             for x = 1:size(face{bid+1,s},1)
%                 for y = 1:size(face{bid+1,s},2)
%                     n = n + 1;
%                 end
%             end
%         end
%     end
%     total_points = n;
%     n = 1;
% 
%     % Loop over all points in faces, check range with kd tree and write to
%     % c-list if coincident
%     for bid = ts_get_bids(g)
%         for s = 1:size(face,2)
%             for x = 1:size(face{bid+1,s},1)
%                 for y = 1:size(face{bid+1,s},2)
% 
%                     XRRT = double(reshape(face{bid+1,s}(x,y,4:6),1,3));
%                     XYZ = reshape(face{bid+1,s}(x,y,7:9),1,3);
%                     IJK = reshape(face{bid+1,s}(x,y,1:3),1,3);
%                     %[index_vals, dist_vals, vector_vals] = kd_rangequery(Tree, double(XRRT), range);
%                     index_vals = kdtree_range(Tree,[XRRT.' - tolerance , XRRT.' + tolerance]);
% 
%                     if ~isempty(index_vals)
%                         for m = 1:length(index_vals)
% 
%                             temp = T(index_vals(m),:);
%                             nxbid = temp(10);
%                             nxIJK = temp(1:3); nxXRRT = temp(4:6); nxXYZ = temp(7:9);
% 
%                             if isequal(IJK,nxIJK) == 0
%                                 if sum(isalmost(XRRT,nxXRRT,tolerance)) == 3
%                                     g{bid+1}.clist{nxbid+1}.ilist = [g{bid+1}.clist{nxbid+1}.ilist ; IJK(1)]; 
%                                     g{bid+1}.clist{nxbid+1}.jlist = [g{bid+1}.clist{nxbid+1}.jlist ; IJK(2)]; 
%                                     g{bid+1}.clist{nxbid+1}.klist = [g{bid+1}.clist{nxbid+1}.klist ; IJK(3)]; 
% 
%                                     g{bid+1}.clist{nxbid+1}.nxilist = [g{bid+1}.clist{nxbid+1}.nxilist ; nxIJK(1)]; 
%                                     g{bid+1}.clist{nxbid+1}.nxjlist = [g{bid+1}.clist{nxbid+1}.nxjlist ; nxIJK(2)]; 
%                                     g{bid+1}.clist{nxbid+1}.nxklist = [g{bid+1}.clist{nxbid+1}.nxklist ; nxIJK(3)]; 
%                                 end
%                             end
%                         end
%                     end
% 
%                     n = n+1;
% 
%                     if mod(n, round(total_points/100)) == 0
%                         disp([num2str(round(n*100/total_points)) ' % done']);
%                     end
%                 end
%             end
%         end
%     end
% 
%     % Unique clists by row
%     for bid = ts_get_bids(g)
%         for c = 1:size(g{bid+1}.clist,2)
%             C = [g{bid+1}.clist{c}.ilist g{bid+1}.clist{c}.jlist g{bid+1}.clist{c}.klist ...
%                 g{bid+1}.clist{c}.nxilist g{bid+1}.clist{c}.nxjlist g{bid+1}.clist{c}.nxklist];
%             C = unique(C,'rows');
%             
%             g{bid+1}.clist{c}.ilist = C(:,1);
%             g{bid+1}.clist{c}.jlist = C(:,2);
%             g{bid+1}.clist{c}.klist = C(:,3);
% 
%             g{bid+1}.clist{c}.nxilist = C(:,4);
%             g{bid+1}.clist{c}.nxjlist = C(:,5);
%             g{bid+1}.clist{c}.nxklist = C(:,6);
%         end
%     end
%     
%     nc_new = 0;
%     % Tidy up empty c-lists, calculate length and assign clids
%     for bid = ts_get_bids(g)
% 
%         for c = size(g{bid+1}.clist,2):-1:1
%             if size(g{bid+1}.clist{c}.ilist,1) == 0
%                 g{bid+1}.clist(c) = [];
%             end
%         end
%         
%         clid = 0;
%         for c = 1:size(g{bid+1}.clist,2)
%             g{bid+1}.clist{c}.attribute.clid = clid;
%             clid = clid+1;
%             g{bid+1}.clist{c}.attribute.n = size(g{bid+1}.clist{c}.ilist,1);
%             nc_new = nc_new + g{bid+1}.clist{c}.attribute.n;
%         end
% 
%     end
% 
% 
% end
% 
% 
%end