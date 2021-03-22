function [chord, bid_os, chi, blades] = ts_find_chord(g, r)

% Find o-mesh blocks
bid_os = ts_find_os(g,r);

chord = zeros(size(bid_os));
chi = zeros(2,size(bid_os,2));
blades = cell(size(bid_os,2),1);

for b = 1:size(bid_os,2)
    bid_o = bid_os(1,b);
    % Find coordinates on nondimensional r = 0.5 line
%     x = zeros(g{bid_o+1}.attribute.ni,1);
%     y = zeros(g{bid_o+1}.attribute.ni,1);
%     z = zeros(g{bid_o+1}.attribute.ni,1);
%     rt = zeros(g{bid_o+1}.attribute.ni,1);
% 
%     if isfield(g{bid_o+1},'y') == 0
%         T = g{bid_o+1}.rt./g{bid_o+1}.r;
%         [Y Z X] = pol2cart(T, g{bid_o+1}.r, g{bid_o+1}.x);
%         g{bid_o+1}.y = Y; g{bid_o+1}.z = Z; g{bid_o+1}.t = T;
%     end
% 
%     for i = 1:g{bid_o+1}.attribute.ni
%         r_nondim = (g{bid_o+1}.r(i,:,1) - min(g{bid_o+1}.r(i,:,1))) / (max(g{bid_o+1}.r(i,:,1)) - min(g{bid_o+1}.r(i,:,1)));
% 
%         x(i,1) = interp1(r_nondim,g{bid_o+1}.x(i,:,1),0.5);
%         y(i,1) = interp1(r_nondim,g{bid_o+1}.y(i,:,1),0.5);
%         z(i,1) = interp1(r_nondim,g{bid_o+1}.z(i,:,1),0.5);
%         rt(i,1) = interp1(r_nondim,g{bid_o+1}.rt(i,:,1),0.5);
%     end

    % Find coordinates on a constant radius line starting at mid height
    J_mid = round(g{bid_o+1}.attribute.nj/2);
    [~, i_min] = min(g{bid_o+1}.x(:,J_mid,1));
    r_mid = (max(g{bid_o+1}.r(i_min,:,1)) + min(g{bid_o+1}.r(i_min,:,1)))/2;
    
    x = zeros(g{bid_o+1}.attribute.ni,1);
    rt = zeros(g{bid_o+1}.attribute.ni,1);
    
    for i = 1:g{bid_o+1}.attribute.ni
        x(i,1) = interp1(g{bid_o+1}.r(i,:,1),g{bid_o+1}.x(i,:,1),r_mid);
        rt(i,1) = interp1(g{bid_o+1}.r(i,:,1),g{bid_o+1}.rt(i,:,1),r_mid);
    end
    if sum(isnan(x)) > 0 || sum(x == 0) > 0
        x = g{bid_o+1}.x(:,J_mid,1);
        rt = g{bid_o+1}.rt(:,J_mid,1);
    end
    % Find maximum distance in mid height blade geometry
%     for i1 = 1:g{bid_o+1}.attribute.ni
%         for i2 = 1:g{bid_o+1}.attribute.ni
%             chord(b) = max(pythag([x(i1,1) y(i1,1) z(i1,1)], [x(i2,1) y(i2,1) z(i2,1)]),chord(b));
%         end
%     end
    % Fit a circle and determine the chord
%     t = bl_parameterise_section([x rt]);
    t = bl_guess_camber([x rt]);
    chord(b) = t.chord;
    chi(1,b) = t.chi_le;
    chi(2,b) = t.chi_te;
    blades{b} = t;

end