function [] = ts_show_probes(g,write_xdmf,filename)
% Shows the location of probe patches

if exist('write_xdmf','var') == 0
    write_xdmf = 0;
end

ts_plot_surface(g,[0 0 0],[],1,1)
hold on

nth = 3;

for n = 1:length(g)
    for nP = 1:length(g{n}.patch)
        if g{n}.patch{nP}.attribute.kind == 8
            p = g{n}.patch{nP}.attribute;
            
            i = unique(round(linspace(double(p.ist+1),double(p.ien),30)));
            j = unique(round(linspace(double(p.jst+1),double(p.jen),30)));
            k = unique(round(linspace(double(p.kst+1),double(p.ken),30)));
            x = squeeze(g{n}.x(i,j,k));
            r = squeeze(g{n}.r(i,j,k));
            rt = squeeze(g{n}.rt(i,j,k));
            
            [y z] = pol2cart(rt./r,r);
            surf(x,z,y,10*ones(size(x)));
            
        end
    end
end

colormap([0 0 0 ; 1 0 0 ]);
alpha(0.4)
shading flat

if write_xdmf == 1
    g_probes = cell(0); bid = 0;
    for n = 1:length(g)
        for nP = 1:length(g{n}.patch)
            if g{n}.patch{nP}.attribute.kind == 8
                p = g{n}.patch{nP}.attribute;

                i = unique(round(linspace(double(p.ist+1),double(p.ien),30)));
                j = unique(round(linspace(double(p.jst+1),double(p.jen),30)));
                k = unique(round(linspace(double(p.kst+1),double(p.ken),30)));
                i = double(p.ist+1):double(p.ien);
                j = double(p.jst+1):double(p.jen);
                k = double(p.kst+1):double(p.ken);
                x = squeeze(g{n}.x(i,j,k));
                r = squeeze(g{n}.r(i,j,k));
                rt = squeeze(g{n}.rt(i,j,k));

                g_temp.x = x; g_temp.r = r; g_temp.rt = rt;
                g_temp.ro = zeros(size(x)); g_temp.rovx = zeros(size(x));
                g_temp.rovr = zeros(size(x)); g_temp.rorvt = zeros(size(x));
                g_temp.roe = zeros(size(x));
                g_temp.attribute.bid = bid; g_temp.attribute.ni = size(x,1);
                g_temp.attribute.nj = size(x,2); g_temp.attribute.nk = size(x,3);
                g_temp.bv.rpm = 0;
                g_probes = [g_probes ; {g_temp}];
                bid = bid+1;
            end
        end
    end
    
    % Write hdf5 and xdmf files for viewing in paraview
    ts_export_paraview(g_probes,filename,[],[],1)
end

end