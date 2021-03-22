%%
function blade = Make_Blade(chi1,chi2,r,c,bladenum,fignum)
figure(fignum); 
for rr=1:size(chi1,1);
    [x,rt] = bladecentroid(chi1(rr),chi2(rr),c(rr),'rotor');
    y = r(rr) .* sin(rt ./ r(rr));
    z = r(rr) .* cos(rt ./ r(rr));
    if bladenum == 2; x=x+0.18; end
    blade.x(rr,:) = x;
    blade.rt(rr,:) = rt;
    blade.r(rr,:) = r(rr);
    blade.t(rr,:) = rt ./ r(rr);
    blade.y(rr,:) = y;
    blade.z(rr,:) = z;
end

subplot(4,6,[3:4,9:10]);hold on; box on
plot(blade.y',-blade.x','-k') ; axis equal
axis([-.2 .2 -.4 .05])
  
subplot(4,6,[5:6,11:12]);hold on; box on
plot3(blade.y',-blade.x',blade.z','-k') ; axis equal
axis([-.2 .2 -.4 .05 0 0.6]); view(3);
title('3D Plot - Rotate Me')

end
