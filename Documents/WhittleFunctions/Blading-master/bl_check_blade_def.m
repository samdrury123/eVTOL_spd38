function [] = bl_check_blade_def(xyz)
% BL_CHECK_BLADE_DEF  Check a blade definition has no overlapping or inverted points
%
%   [] = bl_check_blade_def(xyz)
%   
%   xyz - 3D float array of catersian blade coordinates

% Minimum distance between points
di = sum(squeeze(xyz(2:end,:,:) - xyz(1:end-1,:,:)).^2,3).^0.5;
dj = sum(squeeze(xyz(:,2:end,:) - xyz(:,1:end-1,:)).^2,3).^0.5;

disp(['Min di = ' num2str(min(min(di))) '    Min dj = ' num2str(min(min(dj)))])

% Normal vectors to face
a = xyz(2:end,1:end-1,:) - xyz(1:end-1,2:end,:);
b = xyz(1:end-1,1:end-1,:) - xyz(2:end,2:end,:);

n = cross(a,b);
n = n ./ repmat(sum(n.^2,3).^0.5,[1 1 3]);

Alpha = dot(n(2:end,:,:),n(1:end-1,:,:),3);
Beta = dot(n(:,2:end,:),n(:,1:end-1,:),3);

% Indices
I1 = repmat(reshape(1:size(Alpha,1),[],1),[1 size(Alpha,2)]);
J1 = repmat(reshape(1:size(Alpha,2),1,[]),[size(Alpha,1) 1]);

I2 = repmat(reshape(1:size(Beta,1),[],1),[1 size(Beta,2)]);
J2 = repmat(reshape(1:size(Beta,2),1,[]),[size(Beta,1) 1]);

% Plot change in angle
figure(); 
subplot(2,1,1);
contourf(I1,J1,Alpha);
subplot(2,1,2);
contourf(I2,J2,Beta);

disp(['Min Alpha = ' num2str(min(min(Alpha))) '    Max Alpha = ' num2str(max(max(Alpha)))])
disp(['Min Beta  = ' num2str(min(min(Beta)))  '    Max Beta  = ' num2str(max(max(Beta)))])

end