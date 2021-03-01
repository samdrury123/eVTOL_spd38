function g = ts_set_dir(g)
% Sets directions on all periodic and mixing plane patches

tol = 1e-5;    

% Loop over blocks and patches and pull out data
for bid = ts_get_bids(g)
    for pid = ts_get_pids(bid,g)
        b0 = g{bid+1};
        p0 = b0.patch{pid+1}.attribute;
        
        if p0.kind == 2 || p0.kind == 5
            nxbid = p0.nxbid;
            nxpid = p0.nxpid;
            b1 = g{nxbid+1};
            p1 = b1.patch{nxpid+1}.attribute;
            
            p0.ni = p0.ien - p0.ist;
            p0.nj = p0.jen - p0.jst;
            p0.nk = p0.ken - p0.kst;
            p0.x = squeeze(b0.x(p0.ist+1:p0.ien,p0.jst+1:p0.jen,p0.kst+1:p0.ken));
            p0.r = squeeze(b0.r(p0.ist+1:p0.ien,p0.jst+1:p0.jen,p0.kst+1:p0.ken));
            p0.t = squeeze(b0.rt(p0.ist+1:p0.ien,p0.jst+1:p0.jen,p0.kst+1:p0.ken))./p0.r;
            
            p1.ni = p1.ien - p1.ist;
            p1.nj = p1.jen - p1.jst;
            p1.nk = p1.ken - p1.kst;
            p1.x = squeeze(b1.x(p1.ist+1:p1.ien,p1.jst+1:p1.jen,p1.kst+1:p1.ken));
            p1.r = squeeze(b1.r(p1.ist+1:p1.ien,p1.jst+1:p1.jen,p1.kst+1:p1.ken));
            p1.t = squeeze(b1.rt(p1.ist+1:p1.ien,p1.jst+1:p1.jen,p1.kst+1:p1.ken))./p1.r;            
            
            % If periodic accross the passage on the same block 
            % set directions equal
            if bid == nxbid
                g{bid+1}.patch{pid+1}.attribute.idir = 0;
                g{bid+1}.patch{pid+1}.attribute.jdir = 1;
                g{bid+1}.patch{pid+1}.attribute.kdir = 2;
                
                g{nxbid+1}.patch{nxpid+1}.attribute.idir = 0;
                g{nxbid+1}.patch{nxpid+1}.attribute.jdir = 1;
                g{nxbid+1}.patch{nxpid+1}.attribute.kdir = 2;
                
            elseif p0.kind == 2
                % Deal with mixing planes by assuming they are in the same
                % directions
                g{bid+1}.patch{pid+1}.attribute.idir = -1;
                g{bid+1}.patch{pid+1}.attribute.jdir = 1;
                g{bid+1}.patch{pid+1}.attribute.kdir = 2;
                
                g{nxbid+1}.patch{nxpid+1}.attribute.idir = -1;
                g{nxbid+1}.patch{nxpid+1}.attribute.jdir = 1;
                g{nxbid+1}.patch{nxpid+1}.attribute.kdir = 2;                
                

            else
       
                % Pull out all corner points from both patches
                P0a = [p0.x(1,1) ; p0.r(1,1) ; p0.t(1,1)];
                P0b = [p0.x(end,1) ; p0.r(end,1) ; p0.t(end,1)];
                P0c = [p0.x(end,end) ; p0.r(end,end) ; p0.t(end,end)];
                P0d = [p0.x(1,end) ; p0.r(1,end) ; p0.t(1,end)];

                P1a = [p1.x(1,1) ; p1.r(1,1) ; p1.t(1,1)];
                P1b = [p1.x(end,1) ; p1.r(end,1) ; p1.t(end,1)];
                P1c = [p1.x(end,end) ; p1.r(end,end) ; p1.t(end,end)];
                P1d = [p1.x(1,end) ; p1.r(1,end) ; p1.t(1,end)];   

                if p0.ni == 1
                    Q0j1 = P0b - P1a;
                    Q0j2 = P0c - P0d;
                    Q0k1 = P0d - P0a;
                    Q0k2 = P0c - P0b;
                    if p1.ni == 1
                        Q1j1 = P1b - P1a;
                        Q1j2 = P1c - P1d;
                        Q1k1 = P1c - P1b;
                        Q1k2 = P1d - P1a;

                        % Q0j for jdir
                        if sum(isalmost(Q0j1,Q1j1,tol)) == 3 || sum(isalmost((Q0j1),Q1j2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,Q1j1,tol)) == 3 || sum(isalmost((Q0j2),Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 1;
                        elseif sum(isalmost(Q0j1,-Q1j1,tol)) == 3 || sum(isalmost(Q0j1,-Q1j2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,-Q1j1,tol)) == 3 || sum(isalmost(Q0j2,-Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 4;
                        elseif sum(isalmost(Q0j1,Q1k1,tol)) == 3 || sum(isalmost(Q0j1,Q1k2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,Q1k1,tol)) == 3 || sum(isalmost(Q0j2,Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 2;
                        elseif sum(isalmost(Q0j1,-Q1k1,tol)) == 3 || sum(isalmost(Q0j1,-Q1k2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,-Q1k1,tol)) == 3 || sum(isalmost(Q0j2,-Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 5;
                        else
                            disp('cant find j-direction')
%                             Q0j
%                             Q1j1
%                             Q1j2
%                             Q1k1
%                             Q1k2
%                             g{bid+1}.patch{pid+1}.attribute
                            disp([bid pid]);
                            g{bid+1}.patch{pid+1}.attribute.jdir = 1;
                        end

                        % Q0k for kdir
                        if sum(isalmost(Q0k1,Q1j1,tol)) == 3 || sum(isalmost(Q0k1,Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 1;
                        elseif sum(isalmost(Q0k1,-Q1j1,tol)) == 3 || sum(isalmost(Q0k1,-Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 4;
                        elseif sum(isalmost(Q0k1,Q1k1,tol)) == 3 || sum(isalmost(Q0k1,Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 2;
                        elseif sum(isalmost(Q0k1,-Q1k1,tol)) == 3 || sum(isalmost(Q0k1,-Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 5;
                        else
                            disp('cant find k-direction')
%                             Q0k
%                             Q1j1
%                             Q1j2
%                             Q1k1
%                             Q1k2
%                             g{bid+1}.patch{pid+1}.attribute
                            disp([bid pid]);
                        end

                    elseif p1.nj == 1
                        Q1i1 = P1b - P1a;
                        Q1i2 = P1c - P1d;
                        Q1k1 = P1c - P1b;
                        Q1k2 = P1d - P1a;

                        % Q0j for jdir
                        if sum(isalmost(Q0j1,Q1i1,tol)) == 3 || sum(isalmost((Q0j1),Q1i2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,Q1i1,tol)) == 3 || sum(isalmost((Q0j2),Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 0;
                        elseif sum(isalmost(Q0j1,-Q1i1,tol)) == 3 || sum(isalmost(Q0j1,-Q1i2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,-Q1i1,tol)) == 3 || sum(isalmost(Q0j2,-Q1i2,tol))
                            g{bid+1}.patch{pid+1}.attribute.jdir = 3;
                        elseif sum(isalmost(Q0j1,Q1k1,tol)) == 3 || sum(isalmost(Q0j1,Q1k2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,Q1k1,tol)) == 3 || sum(isalmost(Q0j2,Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 2; 
                        elseif sum(isalmost(Q0j1,-Q1k1,tol)) == 3 || sum(isalmost(Q0j1,-Q1k2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,-Q1k1,tol)) == 3 || sum(isalmost(Q0j2,-Q1k2,tol)) == 3 
                            g{bid+1}.patch{pid+1}.attribute.jdir = 5;
                        else
                            disp('cant find j-direction')
%                             Q0j
%                             Q1i1
%                             Q1i2
%                             Q1k1
%                             Q1k2
%                             g{bid+1}.patch{pid+1}.attribute
                            disp([bid pid]);
                            g{bid+1}.patch{pid+1}.attribute.jdir = 1;
                        end
                        
                        % Q0k for kdir
                        if sum(isalmost(Q0k1,Q1i1,tol)) == 3 || sum(isalmost(Q0k1,Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 0;
                        elseif sum(isalmost(Q0k1,-Q1i1,tol)) == 3 || sum(isalmost(Q0k1,-Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 3;
                        elseif sum(isalmost(Q0k1,Q1k1,tol)) == 3 || sum(isalmost(Q0k1,Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 2;
                        elseif sum(isalmost(Q0k1,-Q1k1,tol)) == 3 || sum(isalmost(Q0k1,-Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 5;
                        else
                            disp('cant find k-direction')
%                             Q0k
%                             Q1i1
%                             Q1i2
%                             Q1k1
%                             Q1k2
%                             g{bid+1}.patch{pid+1}.attribute
                            disp([bid pid]);
                        end     
                        
                    else
                        Q1i1 = P1b - P1a;
                        Q1i2 = P1c - P1d;
                        Q1j1 = P1c - P1b;
                        Q1j2 = P1d - P1a;  
                        
                        % Q0j for jdir
                        if sum(isalmost(Q0j1,Q1i1,tol)) == 3 || sum(isalmost((Q0j1),Q1i2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,Q1i1,tol)) == 3 || sum(isalmost((Q0j2),Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 0;
                        elseif sum(isalmost(Q0j1,-Q1i1,tol)) == 3 || sum(isalmost(Q0j1,-Q1i2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,-Q1i1,tol)) == 3 || sum(isalmost(Q0j2,-Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 3;
                        elseif sum(isalmost(Q0j1,Q1j1,tol)) == 3 || sum(isalmost(Q0j1,Q1j2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,Q1j1,tol)) == 3 || sum(isalmost(Q0j2,Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 1;
                        elseif sum(isalmost(Q0j1,-Q1j1,tol)) == 3 || sum(isalmost(Q0j1,-Q1j2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,-Q1j1,tol)) == 3 || sum(isalmost(Q0j2,-Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 4;
                        else
                            disp('cant find j-direction')
%                             Q0j
%                             Q1i1
%                             Q1i2
%                             Q1j1
%                             Q1j2
%                             g{bid+1}.patch{pid+1}.attribute
                            disp([bid pid]);
                            g{bid+1}.patch{pid+1}.attribute.jdir = 1;
                        end                        
                        
                        % Q0k for kdir
                        if sum(isalmost(Q0k1,Q1i1,tol)) == 3 || sum(isalmost(Q0k1,Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 0;
                        elseif sum(isalmost(Q0k1,-Q1i1,tol)) == 3 || sum(isalmost(Q0k1,-Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 3;
                        elseif sum(isalmost(Q0k1,Q1j1,tol)) == 3 || sum(isalmost(Q0k1,Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 1;
                        elseif sum(isalmost(Q0k1,-Q1j1,tol)) == 3 || sum(isalmost(Q0k1,-Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 4;
                        else
                            disp('cant find k-direction')
%                             Q0k
%                             Q1i1
%                             Q1i2
%                             Q1j1
%                             Q1j2
%                             g{bid+1}.patch{pid+1}.attribute
                            disp([bid pid]);
                        end     
                    end

                    g{bid+1}.patch{pid+1}.attribute.idir = -1;

                elseif p0.nj == 1
                    Q0i1 = P0c - P0d;
                    Q0k1 = P0d - P0a;
                    if p1.ni == 1
                        Q1j1 = P1b - P1a;
                        Q1j2 = P1c - P1d;
                        Q1k1 = P1c - P1b;
                        Q1k2 = P1d - P1a;

                        % Q0i for idir
                        if sum(isalmost(Q0i1,Q1j1,tol)) == 3 || sum(isalmost(Q0i1,Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 1;
                        elseif sum(isalmost(Q0i1,-Q1j1,tol)) == 3 || sum(isalmost(Q0i1,-Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 4;
                        elseif sum(isalmost(Q0i1,Q1k1,tol)) == 3 || sum(isalmost(Q0i1,Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 2;
                        elseif sum(isalmost(Q0i1,-Q1k1,tol)) == 3 || sum(isalmost(Q0i1,-Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 5;
                        else
                            disp('cant find i-direction')
%                             Q0i
%                             Q1j1
%                             Q1j2
%                             Q1k1
%                             Q1k2
%                             g{bid+1}.patch{pid+1}.attribute
                            disp([bid pid]);
                        end

                        % Q0k for kdir
                        if sum(isalmost(Q0k1,Q1j1,tol)) == 3 || sum(isalmost(Q0k1,Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 1;
                        elseif sum(isalmost(Q0k1,-Q1j1,tol)) == 3 || sum(isalmost(Q0k1,-Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 4;
                        elseif sum(isalmost(Q0k1,Q1k1,tol)) == 3 || sum(isalmost(Q0k1,Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 2;
                        elseif sum(isalmost(Q0k1,-Q1k1,tol)) == 3 || sum(isalmost(Q0k1,-Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 5;
                        else
                            disp('cant find k-direction')

%                             Q0k
%                             Q1j1
%                             Q1j2
%                             Q1k1
%                             Q1k2
%                             g{bid+1}.patch{pid+1}.attribute
                            disp([bid pid]);
                        end
                        
                    elseif p1.nj ==1
                        Q1i1 = P1b - P1a;
                        Q1i2 = P1c - P1d;
                        Q1k1 = P1c - P1b;
                        Q1k2 = P1d - P1a;
                        
                        % Q0i for idir
                        if sum(isalmost(Q0i1,Q1i1,tol)) == 3 || sum(isalmost(Q0i1,Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 0;
                        elseif sum(isalmost(Q0i1,-Q1i1,tol)) == 3 || sum(isalmost(Q0i1,-Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 3;
                        elseif sum(isalmost(Q0i1,Q1k1,tol)) == 3 || sum(isalmost(Q0i1,Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 2;
                        elseif sum(isalmost(Q0i1,-Q1k1,tol)) == 3 || sum(isalmost(Q0i1,-Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 5;
                        else
                            disp('cant find i-direction')
%                             Q0i
%                             Q1i1
%                             Q1i2
%                             Q1k1
%                             Q1k2
%                             g{bid+1}.patch{pid+1}.attribute
                            disp([bid pid]);
                        end       

                        % Q0k for kdir
                        if sum(isalmost(Q0k1,Q1i1,tol)) == 3 || sum(isalmost(Q0k1,Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 0;
                        elseif sum(isalmost(Q0k1,-Q1i1,tol)) == 3 || sum(isalmost(Q0k1,-Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 3;
                        elseif sum(isalmost(Q0k1,Q1k1,tol)) == 3 || sum(isalmost(Q0k1,Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 2;
                        elseif sum(isalmost(Q0k1,-Q1k1,tol)) == 3 || sum(isalmost(Q0k1,-Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 5;
                        else
                            disp('cant find k-direction')
%                             Q0k
%                             Q1i1
%                             Q1i2
%                             Q1k1
%                             Q1k2
%                             g{bid+1}.patch{pid+1}.attribute
                            disp([bid pid]);
                        end                        
                        
                    else
                        Q1i1 = P1b - P1a;
                        Q1i2 = P1c - P1d;
                        Q1j1 = P1c - P1b;
                        Q1j2 = P1d - P1a;                        

                        % Q0i for idir
                        if sum(isalmost(Q0i1,Q1i1,tol)) == 3 || sum(isalmost(Q0i1,Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 0;
                        elseif sum(isalmost(Q0i1,-Q1i1,tol)) == 3 || sum(isalmost(Q0i1,-Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 3;
                        elseif sum(isalmost(Q0i1,Q1j1,tol)) == 3 || sum(isalmost(Q0i1,Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 1;
                        elseif sum(isalmost(Q0i1,-Q1j1,tol)) == 3 || sum(isalmost(Q0i1,-Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 4;
                        else
                            disp('cant find i-direction')
%                             Q0i
%                             Q1i1
%                             Q1i2
%                             Q1j1
%                             Q1j2
%                             g{bid+1}.patch{pid+1}.attribute
                            disp([bid pid]);
                        end    
                        
                        % Q0k for kdir
                        if sum(isalmost(Q0k1,Q1i1,tol)) == 3 || sum(isalmost(Q0k1,Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 0;
                        elseif sum(isalmost(Q0k1,-Q1i1,tol)) == 3 || sum(isalmost(Q0k1,-Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 3;
                        elseif sum(isalmost(Q0k1,Q1j1,tol)) == 3 || sum(isalmost(Q0k1,Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 1;
                        elseif sum(isalmost(Q0k1,-Q1j1,tol)) == 3 || sum(isalmost(Q0k1,-Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.kdir = 4;
                        else
                            disp('cant find k-direction')
%                             Q0k
%                             Q1i1
%                             Q1i2
%                             Q1j1
%                             Q1j2
%                             g{bid+1}.patch{pid+1}.attribute
                            disp([bid pid]);
                        end                                        
                    end

                    g{bid+1}.patch{pid+1}.attribute.jdir = -1;

                else
                    Q0i1 = P0c - P0d;
                    Q0i2 = P0b - P0a;
                    Q0j1 = P0d - P0a;
                    Q0j2 = P0c - P0b;
                    
                    if p1.ni == 1
                        Q1j1 = P1b - P1a;
                        Q1j2 = P1c - P1d;
                        Q1k1 = P1c - P1b;
                        Q1k2 = P1d - P1a;
                        
                        % Q0i for idir
                        if sum(isalmost(Q0i1,Q1j1,tol)) == 3 || sum(isalmost(Q0i1,Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 1;
                        elseif sum(isalmost(Q0i1,-Q1j1,tol)) == 3 || sum(isalmost(Q0i1,-Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 4;
                        elseif sum(isalmost(Q0i1,Q1k1,tol)) == 3 || sum(isalmost(Q0i1,Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 2;
                        elseif sum(isalmost(Q0i1,-Q1k1,tol)) == 3 || sum(isalmost(Q0i1,-Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 5;
                        else
                            disp('cant find i-direction')
%                             Q0i
%                             Q1j1
%                             Q1j2
%                             Q1k1
%                             Q1k2
%                             g{bid+1}.patch{pid+1}.attribute
                            disp([bid pid]);
                        end
                        
                        % Q0j for jdir
                        if sum(isalmost(Q0j1,Q1j1,tol)) == 3 || sum(isalmost((Q0j1),Q1j2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,Q1j1,tol)) == 3 || sum(isalmost((Q0j2),Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 1;
                        elseif sum(isalmost(Q0j1,-Q1j1,tol)) == 3 || sum(isalmost(Q0j1,-Q1j2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,-Q1j1,tol)) == 3 || sum(isalmost(Q0j2,-Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 4;
                        elseif sum(isalmost(Q0j1,Q1k1,tol)) == 3 || sum(isalmost(Q0j1,Q1k2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,Q1k1,tol)) == 3 || sum(isalmost(Q0j2,Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 2;
                        elseif sum(isalmost(Q0j1,-Q1k1,tol)) == 3 || sum(isalmost(Q0j1,-Q1k2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,-Q1k1,tol)) == 3 || sum(isalmost(Q0j2,-Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 5;
                        else
                            disp('cant find j-direction')
%                             Q0j
%                             Q1j1
%                             Q1j2
%                             Q1k1
%                             Q1k2
%                             g{bid+1}.patch{pid+1}.attribute
                            disp([bid pid]);
                            g{bid+1}.patch{pid+1}.attribute.jdir = 1;
                        end     
                        
                    elseif p1.nj == 1
                        Q1i1 = P1b - P1a;
                        Q1i2 = P1c - P1d;
                        Q1k1 = P1c - P1b;
                        Q1k2 = P1d - P1a;   
                        
                        % Q0i for idir
                        if sum(isalmost(Q0i1,Q1i1,tol)) == 3 || sum(isalmost(Q0i1,Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 0;
                        elseif sum(isalmost(Q0i1,-Q1i1,tol)) == 3 || sum(isalmost(Q0i1,-Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 3;
                        elseif sum(isalmost(Q0i1,Q1k1,tol)) == 3 || sum(isalmost(Q0i1,Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 2;
                        elseif sum(isalmost(Q0i1,-Q1k1,tol)) == 3 || sum(isalmost(Q0i1,-Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 5;
                        else
                            disp('cant find i-direction')
%                             Q0i
%                             Q1i1
%                             Q1i2
%                             Q1k1
%                             Q1k2
%                             g{bid+1}.patch{pid+1}.attribute
                            disp([bid pid]);
                        end
                        
                        % Q0j for jdir
                        if sum(isalmost(Q0j1,Q1i1,tol)) == 3 || sum(isalmost((Q0j1),Q1i2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,Q1i1,tol)) == 3 || sum(isalmost((Q0j2),Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 0;
                        elseif sum(isalmost(Q0j1,-Q1i1,tol)) == 3 || sum(isalmost(Q0j1,-Q1i2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,-Q1i1,tol)) == 3 || sum(isalmost(Q0j2,-Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 3;
                        elseif sum(isalmost(Q0j1,Q1k1,tol)) == 3 || sum(isalmost(Q0j1,Q1k2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,Q1k1,tol)) == 3 || sum(isalmost(Q0j2,Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 2;
                        elseif sum(isalmost(Q0j1,-Q1k1,tol)) == 3 || sum(isalmost(Q0j1,-Q1k2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,-Q1k1,tol)) == 3 || sum(isalmost(Q0j2,-Q1k2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 5;
                        else
                            disp('cant find j-direction')
%                             Q0j
%                             Q1i1
%                             Q1i2
%                             Q1k1
%                             Q1k2
%                             g{bid+1}.patch{pid+1}.attribute
                            disp([bid pid]);
                            g{bid+1}.patch{pid+1}.attribute.jdir = 1;
                        end                             

                    else
                        Q1i1 = P1b - P1a;
                        Q1i2 = P1c - P1d;
                        Q1j1 = P1c - P1b;
                        Q1j2 = P1d - P1a;   
                        
                        % Q0i for idir
                        if sum(isalmost(Q0i1,Q1i1,tol)) == 3 || sum(isalmost(Q0i1,Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 0;
                        elseif sum(isalmost(Q0i1,-Q1i1,tol)) == 3 || sum(isalmost(Q0i1,-Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 3;
                        elseif sum(isalmost(Q0i1,Q1j1,tol)) == 3 || sum(isalmost(Q0i1,Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 1;
                        elseif sum(isalmost(Q0i1,-Q1j1,tol)) == 3 || sum(isalmost(Q0i1,-Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.idir = 4;
                        else
                            disp('cant find i-direction')
%                             Q0i
%                             Q1i1
%                             Q1i2
%                             Q1j1
%                             Q1j2
%                             g{bid+1}.patch{pid+1}.attribute
                            disp([bid pid]);
                        end
                        
                        % Q0j for jdir
                        if sum(isalmost(Q0j1,Q1i1,tol)) == 3 || sum(isalmost((Q0j1),Q1i2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,Q1i1,tol)) == 3 || sum(isalmost((Q0j2),Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 0;
                        elseif sum(isalmost(Q0j1,-Q1i1,tol)) == 3 || sum(isalmost(Q0j1,-Q1i2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,-Q1i1,tol)) == 3 || sum(isalmost(Q0j2,-Q1i2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 3;
                        elseif sum(isalmost(Q0j1,Q1j1,tol)) == 3 || sum(isalmost(Q0j1,Q1j2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,Q1j1,tol)) == 3 || sum(isalmost(Q0j2,Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 1;
                        elseif sum(isalmost(Q0j1,-Q1j1,tol)) == 3 || sum(isalmost(Q0j1,-Q1j2,tol)) == 3 ...
                            || sum(isalmost(Q0j2,-Q1j1,tol)) == 3 || sum(isalmost(Q0j2,-Q1j2,tol)) == 3
                            g{bid+1}.patch{pid+1}.attribute.jdir = 4;
                        else
                            disp('cant find j-direction')
%                             Q0j
%                             Q1i1
%                             Q1i2
%                             Q1j1
%                             Q1j2
%                             g{bid+1}.patch{pid+1}.attribute
                            disp([bid pid]);
                            g{bid+1}.patch{pid+1}.attribute.jdir = 1;
                        end
                        
                    end
                    
                    g{bid+1}.patch{pid+1}.attribute.kdir = -1;

                end

            end
        end
    end
end

end