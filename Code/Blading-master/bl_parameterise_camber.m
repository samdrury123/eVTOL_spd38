function t = bl_parameterise_camber(xrt_cam)
% BL_PARAMETERISE_CAMBER  Decompose camber line into governing parameters

% Inlet and exit metal angles
t.chi_le = atand(diff(xrt_cam(1:2,2),1,1) ./ diff(xrt_cam(1:2,1),1,1));
t.chi_te = atand(diff(xrt_cam(end-1:end,2),1,1) ./ diff(xrt_cam(end-1:end,1),1,1));

% Non-dimensional distance along the camberline
t.s_cl = [0 ; cumsum(sum(diff(xrt_cam,1,1).^2,2).^0.5,1)]; 
t.s_cl = t.s_cl ./ t.s_cl(end,:);

% Non-dimensional camberline
t.chi = atand(grad_mg(xrt_cam(:,1),xrt_cam(:,2)));
t.cam = (t.chi - t.chi_te) ./ (t.chi_le - t.chi_te);

% Camberline kinky-ness
q = t.s_cl > 0.1 & t.s_cl < 0.9;
p = polyfit(t.s_cl(q),t.cam(q),4);
t.qcam = p(1) - 2*p(2) + p(3);

% Initial camber gradient
q = t.s_cl >= 0.02 & t.s_cl < 0.07;
t.p_le = polyfit(t.s_cl(q),t.cam(q),2);
t.dcam_le = polyval(polyder(t.p_le),0);

% Final camber gradient
q = t.s_cl > 0.93 & t.s_cl <= 0.98;
t.p_te = polyfit(t.s_cl(q),t.cam(q),2);
t.dcam_te = polyval(polyder(t.p_te),1);    


end