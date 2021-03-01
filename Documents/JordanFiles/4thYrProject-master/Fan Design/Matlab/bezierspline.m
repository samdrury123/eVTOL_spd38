function [bezier] = bezierspline(z,R,final_sections)
% Bezier curve script
sections_per_bezie_task = 3; %cant do other

%%
% % % % % r=linspace(0,1,11));
% % % % % bezierbasis=bernsteinMatrix(subsec,r)
% % % % % B are the B-spline control points, C1 and C2 are the control points for
% % % % % the Bezier curves 1 and 2 that make up the B-spline.
% % % % % B=[R z];
% % % % % B2( = -1/4*B(1,:)+3/2*B(2,:)-1/4*B(3,:);
% % % % % % First Bezier branch
% % % % % Cn1=[B(1,:)
% % % % %     2/3*B(1,:)+1/3*B2
% % % % %     1/3*B(1,:)+2/3*B2
% % % % %     B(2,:)];
% % % % % % Second Bezier branch
% % % % % Cn2=[B(2,:)
% % % % %     2/3*B2+1/3*B(3,:)
% % % % %     1/3*B2+2/3*B(3,:)
% % % % %     B(3,:)];
% % % % % % First derivative at joining point
% % % % % % S1 = 3*(Cn1(4,:) - Cn1(3,:))
% % % % % % S2 = 3*(Cn2(2,:) - Cn2(1,:))
% % % % % % Generate Bezier curves from bezierbasis
% % % % % bezier1=bezierbasis*Cn1;
% % % % % bezier2=bezierbasis*Cn2;
% % % % % bezier1=[bezier1(:,2) bezier1(:,1)];
% % % % % bezier2=[bezier2(:,2) bezier2(:,1)];
%%

if length(z)==3
    r=linspace(0,1,ceil(final_sections/(length(z)-1)));
    bezierbasis=bernsteinMatrix(sections_per_bezie_task,r);
    B=[R z];
    B2(1,:) = -1/4*B(1,:)+3/2*B(2,:)-1/4*B(3,:);
    % First Bezier branch
    Cn1=[B(1,:)
        2/3*B(1,:)+1/3*B2
        1/3*B(1,:)+2/3*B2
        B(2,:)];
    %Second Bezier branch
    Cn2=[B(2,:)
        2/3*B2+1/3*B(3,:)
        1/3*B2+2/3*B(3,:)
        B(3,:)];
    bezier.A1=fliplr(bezierbasis*Cn1);
    bezier.A2=fliplr(bezierbasis*Cn2);
    
elseif length(z)==5
    r=linspace(0,1,ceil(final_sections/(length(z)-1)));
    bezierbasis=bernsteinMatrix(sections_per_bezie_task,r);
      B_1=[R(1:3) z(1:3)];
    B2(1,:) = -1/4*B_1(1,:)+3/2*B_1(2,:)-1/4*B_1(3,:);
    
    Cn1=[B_1(1,:)
        2/3*B_1(1,:)+1/3*B2(1,:)
        1/3*B_1(1,:)+2/3*B2(1,:)
        B_1(2,:)];
    
    Cn2=[B_1(2,:)
        2/3*B2(1,:)+1/3*B_1(3,:)
        1/3*B2(1,:)+2/3*B_1(3,:)
        B_1(3,:)];
    
    
    bezier.A1=fliplr(bezierbasis*Cn1);
    bezier.A2=fliplr(bezierbasis*Cn2);
    
    B_2=[R(3:5) z(3:5)];
    B2(2,:) = -1/4*B_2(1,:)+3/2*B_2(2,:)-1/4*B_2(3,:);
    
    Cn1=[B_2(1,:)
        2/3*B_2(1,:)+1/3*B2(2,:)
        1/3*B_2(1,:)+2/3*B2(2,:)
        B_2(2,:)];
    
    Cn2=[B_2(2,:)
        2/3*B2(2,:)+1/3*B_2(3,:)
        1/3*B2(2,:)+2/3*B_2(3,:)
        B_2(3,:)];
    
    bezier.A3=fliplr(bezierbasis*Cn1);
    bezier.A4=fliplr(bezierbasis*Cn2);
    
elseif length(z)==7
    r=linspace(0,1,ceil(final_sections/(length(z)-1)));
    bezierbasis=bernsteinMatrix(sections_per_bezie_task,r);
    B_1=[R(1:3) z(1:3)];
    B2(1,:) = -1/4*B_1(1,:)+3/2*B_1(2,:)-1/4*B_1(3,:);
    
    Cn1=[B_1(1,:)
        2/3*B_1(1,:)+1/3*B2(1,:)
        1/3*B_1(1,:)+2/3*B2(1,:)
        B_1(2,:)];
    
    Cn2=[B_1(2,:)
        2/3*B2(1,:)+1/3*B_1(3,:)
        1/3*B2(1,:)+2/3*B_1(3,:)
        B_1(3,:)];
    
    
    bezier.A1=fliplr(bezierbasis*Cn1);
    bezier.A2=fliplr(bezierbasis*Cn2);
    
    B_2=[R(3:5) z(3:5)];
    B2(2,:) = -1/4*B_2(1,:)+3/2*B_2(2,:)-1/4*B_2(3,:);
    
    Cn1=[B_2(1,:)
        2/3*B_2(1,:)+1/3*B2(2,:)
        1/3*B_2(1,:)+2/3*B2(2,:)
        B_2(2,:)];
    
    Cn2=[B_2(2,:)
        2/3*B2(2,:)+1/3*B_2(3,:)
        1/3*B2(2,:)+2/3*B_2(3,:)
        B_2(3,:)];
    
    bezier.A3=fliplr(bezierbasis*Cn1);
    bezier.A4=fliplr(bezierbasis*Cn2);
    
    B_3=[R(5:7) z(5:7)];
    B2(3,:) = -1/4*B_3(1,:)+3/2*B_3(2,:)-1/4*B_3(3,:);
    
    Cn1=[B_2(1,:)
        2/3*B_2(1,:)+1/3*B2(2,:)
        1/3*B_2(1,:)+2/3*B2(2,:)
        B_2(2,:)];
    
    Cn2=[B_2(2,:)
        2/3*B2(2,:)+1/3*B_2(3,:)
        1/3*B2(2,:)+2/3*B_2(3,:)
        B_2(3,:)];
    
    bezier.A5=fliplr(bezierbasis*Cn1);
    bezier.A6=fliplr(bezierbasis*Cn2);
    
else
    disp('Oups')
    Break
end

% % % %  for i=1:sections_per_bezie_task-1
% % % %      %Bezier branch
% % % %      if mod(i,2)
% % % %          Cn=[B(i,:)
% % % %              2/3*B(i,:)+1/3*B2(ceil(i/2),:)
% % % %              1/3*B(i,:)+2/3*B2(ceil(i/2),:)
% % % %              B(i+1,:)];
% % % %      else
% % % %          Cn=[B(i,:)
% % % %              2/3*B2(ceil(i/2),:)+1/3*B(i+1,:)
% % % %              1/3*B2(ceil(i/2),:)+2/3*B(i+1,:)
% % % %              B(i+1,:)];
% % % %      end
% % % %      % Generate Bezier curves from bezierbasis
% % % %      bezier.(['A' num2str(i)])=fliplr(bezierbasis*Cn);
% % % % end
end