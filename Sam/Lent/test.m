clear; close;

philist = [0.5:0.1:0.8];
sigmalist = [0.8:0.1:1.2];
Patient([4 5]) = struct();
Patient(2,3).g = 5;
for pp=1:size(philist,2)
    for ss=1:size(sigmalist,2)
        patient(pp,ss).phi = philist(pp);
        patient(pp,ss).sigma = sigmalist(ss);
    end
end
alt = 1000;
atm = Altitude(alt); % At sea level

% a([4 5]) = struct();
% p=linspace(1,4,4);
% q=linspace(0.1,0.4,4);
% for i=1:4
%     for j=1:5
%         a(i,j).x.c = i*j;
%         a(i,j).x.d = i*i;
%         a(i,j).y = i;
%     end
% end
% figure; hold on;
% contourf(p,q,a.x); colorbar;

% Nr=31;
% r1 = linspace(20,60,Nr)';
% r2 = linspace(20,60,Nr);
% 
% r5 = linspace(0,1,5); y5 = zeros(5,1);
% r6 = linspace(0,1,6); y6 = zeros(6,1);
% r7 = linspace(0,1,7); y7 = zeros(7,1);
% r8 = linspace(0,1,8); y8 = zeros(8,1);
% figure; hold on; grid off; box on;
% plot(r5,y5, 'ro');
% plot(r6,y6,'ko'); 
% plot(r7,y7,'bo');
% plot(r8,y8,'go');
% rm=(0.5*(0.2^2+0.6^2))^0.5;
% rm_n = (rm-0.2)/(0.6-0.2);
% plot(rm_n,0,'rx');
% legend;

x=0.1:0.1:1;
y=x.^2;
z=2*x;
z2=2.*x;
trial = [0.35,0.45];
Y = interp1(x,y,trial);

NRF.d.a.phi=0.8;
NRF.d.a.sigma=1.1;
NRF.g.rh=0.5;
CRF.d.phi=0.9;
CRF.d.sigma=1.2;
CRF.g.rc=0.7;
a=1;
st=[NRF CRF];

if ~isnan(a)
    for rr=1:2
        b(rr)=1;
    end
end


a=[1 2 3 4 5 6; 7 8 9 10 11 12];
sum(a(2:5,2))




