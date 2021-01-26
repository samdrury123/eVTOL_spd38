function [steps1,steps2,sections,pointspersection,XYZ1,XYZ2] = KON_Make_Blade(chi1,chi2,Radii,c,bladetype,file_name)
% USE bezier spline to convert the 3 section to 21 subsections =(21+1)/2
final_sections = 21;
N_bezier = length(chi1)-1;
[bezier_chi1]       = bezierspline(chi1,Radii,final_sections);
[bezier_chi2]       = bezierspline(chi2,Radii,final_sections);
[bezier_c]             = bezierspline(c,Radii,final_sections);
[R,T,Z] =  stackblade(bezier_chi1,bezier_chi2,bezier_c,N_bezier,bladetype,file_name);
S=size(R);
pointspersection = ceil(S(1)/2);
sections         = S(2);

for i=1:sections
    steps1(:,1,i)=R(1:pointspersection,i);
    steps1(:,2,i)=T(1:pointspersection,i);
    steps1(:,3,i)=Z(1:pointspersection,i);
    steps2(:,1,i)=R(pointspersection:end,i);
    steps2(:,2,i)=T(pointspersection:end,i);
    steps2(:,3,i)=Z(pointspersection:end,i);
    [XYZ1(:,1,i), XYZ1(:,2,i)]=pol2cart(T(1:pointspersection,i),R(1:pointspersection,i));
    XYZ1(:,3,i)=Z(1:pointspersection,i);
    [XYZ2(:,1,i), XYZ2(:,2,i)]=pol2cart(T(pointspersection:end,i),R(pointspersection:end,i));
    XYZ2(:,3,i)=Z(pointspersection:end,i);
end
