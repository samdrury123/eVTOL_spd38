function [xblade,yblade] = bladecentroid(chi1d,chi2d,c,bladetype)
% Function that uses the blade function and translates the blade so that
% its area centroid is at the origin
chi1=chi1d*pi/180;
chi2=chi2d*pi/180;
camber=chi1-chi2;
r=c/camber;
% Shape / class function blade definition
xle=linspace(0,0.05,126);
xbody=linspace(0.05,1,96);
x=[xle xbody(2:end)];
C=sqrt(x).*(1-x);
if strcmp(bladetype,'rotor') == 1
    factor = 1.2;                   % Default = 1
    S=factor*(0.01*cos(2*pi*x)-0.11*x + 0.13);% This shape function & zte give
    zte= 0.006*2;                       % blade with tmax=8%c at x=0.19c
    Ntopte=31;                       % T/E at i=244 with 126 LE points+96 body points
    Nbotte=11;
    ii=130;
elseif strcmp(bladetype,'stator') == 1
%     S=7/8*(0.01*cos(2*pi*x)-0.11*x + 0.13); % tmax=7%c at x=0.2c
    factor = 1.5;
    S= factor * 1.2*(0.01*cos(2*pi*x)-0.11*x + 0.13); % tmax=12%c at x=0.2c
%     zte=0.006;                              % T/E at i=245 with 126 LE points+96 body points 1.2%c
    zte=0.006*1.5;                              %1.8%c (0.006 * 1.5)
    Ntopte=25;
    Nbotte=25;
    ii=131;
end
t=S.*C+x*zte; 
% Convert thickness and distance along chord to dimensional values
t=t*c;
x=x*c;
% Find coordinates of top blade surface & camber line
theta=x/r;
xcamber=r*sin(chi1)-r*sin(chi1-theta);
ycamber=r*cos(chi1-theta)-r*cos(chi1);
xbladetop = r*sin(chi1)-(r+t).*sin(chi1-theta);
ybladetop = (r+t).*cos(chi1-theta)-r*cos(chi1);
% Bottom blade surface
xbladebot = r*sin(chi1)-(r-t).*sin(chi1-theta);
ybladebot = (r-t).*cos(chi1-theta)-r*cos(chi1);
% Trailing Edge circle
xdiffend=xbladetop(end)-xbladetop(end-1);
ydiffend=ybladetop(end)-ybladetop(end-1);
xdia=xbladetop(end)-xbladebot(end);
ydia=ybladetop(end)-ybladebot(end);
d=0:0.0001:0.1;
xcentre=xcamber(end)-d*ydia;
ycentre=ycamber(end)+d*xdia;
xrad=xbladetop(end)-xcentre;
yrad=ybladetop(end)-ycentre;
scalarprod=xdiffend*xrad+ydiffend*yrad;
[~,icentre]=min(abs(scalarprod));
xcircTE=xcentre(icentre);
ycircTE=ycentre(icentre);
radiuste=sqrt((xbladetop(end)-xcircTE)^2+(ybladetop(end)-ycircTE)^2);
% Update the camber line to add last point consisting in TE
exte=ydia/sqrt(xdia^2+ydia^2);
eyte=-xdia/sqrt(xdia^2+ydia^2);
xte=xcircTE+radiuste*exte;
yte=ycircTE+radiuste*eyte;
xcamber=[xcamber xte];
ycamber=[ycamber yte];
% Work out the length of camber line and set this equal to chord c
dx=diff(xcamber);
dy=diff(ycamber);
ds=sqrt(dx.^2+dy.^2);
s=sum(ds);
scalefactor=s/c;
% Now add the circular TE
xtopte=linspace(xbladetop(end),xcircTE+radiuste,Ntopte);
xbotte=linspace(xbladebot(end),xcircTE+radiuste,Nbotte);
ytopte=real(ycircTE+sqrt(radiuste^2-(xtopte-xcircTE).^2));
ybotte=real(ycircTE-sqrt(radiuste^2-(xbotte-xcircTE).^2));
% Group coordinates into xblade, yblade
xblade=[xbladetop xtopte(2:end) xbotte(end-1:-1:2) xbladebot(end:-1:1)];
yblade=[ybladetop ytopte(2:end) ybotte(end-1:-1:2) ybladebot(end:-1:1)];
% Translate the blade coordinates to have centroid at (0,0) and then scale
xblade=xblade-xcamber(ii);  % add stacking line lean/sweep f(R)
yblade=yblade-ycamber(ii);
%xblade=xblade/scalefactor;
%yblade=yblade/scalefactor;
% Plot blade surface
if strcmp(bladetype,'stator') == 1
yblade=-yblade;
end

% plot(xblade,yblade)
% hold on
% axis equal

end
