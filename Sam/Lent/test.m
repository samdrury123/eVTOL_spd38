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


a=[1 2 3 4 5 6; 7 8 9 10 11 12]';
sum(a(2:5,2));

c=[1 2 3 4 5]/10;
endwallfun = @(x,r) x.^2 + c(r);
q = integral(@(x) endwallfun(x,5),0,2);
cas(1:5) = 4;




options = optimset('PlotFcns',@optimplotfval);
fun = @(x)100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
x0 = [-1.2,1];
x = fminsearch(fun,x0,options)














