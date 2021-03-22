function [box_choice,dS] = really_useful_box_selector(xyz,z_weight,fix_orientation,open_front)
% REALLY_USEFUL_BOX_SELECTOR  A really useful function to choose you a really useful box
%
%   [box_choice,dS] = REALLY_USEFUL_BOX_SELECTOR(xyz,z_weight,fix_orientation,open_front)
%
%   xyz - is a vector of your minimum internal box size in mm
%   z_weight - optional argument to define weighting to put on the z dimension
%   fix_orientation - optional argument, either 0 or 1 to allow any orientation of your parts
%   open_front - optional argument, either 0 or 1 to allow or reject open fronted designs
%   dS - representitative value of the left over space

% Default to weight z dimension by half
if exist('z_weight','var') == 0
    z_weight = 0.5;
end

% Default to allow any orientation, if not fix z
if exist('fix_orientation','var') == 0
    fix_orientation = 0;
end

% Default to not use the really useless open front boxes
if exist('open_front','var') == 0
    open_front = 0;
end

% Get the box sizes from sub function
box_sizes = get_sizes;

% Remove open front boxes if necessary
if open_front == 0
    i = regexp(box_sizes(:,1),'open','match');
    box_sizes = box_sizes(cellfun(@isempty,i) == 1,:);
end

% Remove folding boxes
i = regexp(box_sizes(:,1),'folding','match');
box_sizes = box_sizes(cellfun(@isempty,i) == 1,:);

% Parse internal dimension data
[nb,nf] = size(box_sizes);
xyz_temp = regexp(box_sizes(:,3),' x ','split');
xyz_box = reshape(str2double([xyz_temp{:}]),[3 nb])';

% Calculate differnces in all dimensions
if fix_orientation == 0
    dxyz = cat(3,xyz_box,xyz_box(:,[2 3 1]),xyz_box(:,[3 1 2])) - permute(repmat(xyz,[nb 1 3]),[1 3 2]);
else
    dxyz = permute(cat(3,xyz_box,xyz_box(:,[2 1 3])) - repmat(xyz,[nb 1 2]),[1 3 2]);
end

% Look for positive values in the third dimension
q_fit = sum(dxyz >= 0,3) == 3;
no = size(dxyz,2); nq = sum(q_fit(:));
box_fit = reshape(repmat(reshape(box_sizes,[nb 1 nf]),[1 no 1]),[no*nb nf]); 
box_fit = box_fit(reshape(q_fit,[nb*no 1]) == 1,:);

% Characterise the fit by the square of the space in all dimensions
dxyz = reshape(dxyz(repmat(q_fit,[1 1 3])),[nq 3]);
dS = sum(dxyz.^2 .* repmat([1 1 z_weight],[nq 1]),2).^0.5;

% Sort boxes and choose the best
[dS,i] = sort(dS);
dS = dS(1); box_choice = box_fit(i(1),:);

end

function box_sizes = get_sizes

% Define available box size and information
box_sizes = {'0.07 litre Really Useful Box','090 x 065 x 030','061 x 048 x 024','24'
    '0.14 litre Really Useful Box','090 x 065 x 055','064 x 048 x 042','35'
    '0.2 litre Really Useful Box','120 x 085 x 045','090 x 065 x 032','46'
    '0.3 litre Really Useful Box','120 x 085 x 065','090 x 065 x 055','54'
    '0.35 litre Really Useful Box','155 x 100 x 040','120 x 080 x 035','67'
    '0.55 litre Really Useful Box','220 x 100 x 040','185 x 080 x 030','102'
    '0.7 litre Really Useful Box','155 x 100 x 080','120 x 080 x 070','90'
    '0.75 litre Really Useful Box','195 x 135 x 055','150 x 110 x 045','143'
    '0.8 litre Really Useful Box','355 x 100 x 040','320 x 080 x 030','155'
    '0.9 litre Really Useful Box','220 x 100 x 070','185 x 080 x 060','130'
    '1.5 litre Really Useful Box','355 x 100 x 070','320 x 080 x 060','195'
    '1.6 litre Really Useful Box','195 x 135 x 110','150 x 110 x 095','200'
    '1.7 litre folding Really Useful Box','180 x 135 x 090','160 x 115 x 087','105'
    '1.75 litre Really Useful Box','245 x 180 x 078','205 x 150 x 060','306'
    '2.1 litre Really Useful Box','240 x 130 x 125','190 x 100 x 104','241'
    '2.5 litre Really Useful Box','340 x 200 x 70','290 x 160 x 50','365'
    '3 litre Really Useful Box','245 x 180 x 160','190 x 143 x 120','500'
    '4 litre Really Useful Box','395 x 255 x 088','348 x 220 x 068','700'
    '5 litre Really Useful Box','340 x 200 x 125','285 x 160 x 105','540'
    '5 litre XL Really Useful Box','340 x 200 x 150','285 x 160 x 135','570'
    '6 litre Really Useful Box','465 x 270 x 085','400 x 225 x 070','890'
    '6.5 litre Really Useful Box','430 x 180 x 160','380 x 143 x 120','668'
    '7 litre Really Useful Box','400 x 350 x 085','340 x 307 x 070','800'
    '8 litre open front Really Useful Box','340 x 200 x 175','290 x 165 x 153','590'
    '9 litre Really Useful Box','395 x 255 x 155','335 x 210 x 140','900'
    '"9 litre XL Really Useful Box (7"" small)"','395 x 255 x 205','335 x 210 x 180','1000'
    '10 litre Really Useful Box','520 x 340 x 085','460 x 298 x 070','1200'
    '11 litre Really Useful Box','456 x 356 x 120','375 x 310 x 091','1040'
    '11 litre XL Really Useful Box','455 x 356 x 148','375 x 310 x 125','1215'
    '12 litre Really Useful Box','465 x 270 x 150','415 x 235 x 140','998'
    '14 litre open front Really Useful Box','395 x 255 x 210','336 x 214 x 190','950'
    '18 litre Really Useful Box','480 x 390 x 200','395 x 335 x 170','1300'
    '"18 litre XL Really Useful Box (7"" medium)"','480 x 390 x 230','395 x 335 x 200','1500'
    '19 litre Really Useful Box','395 x 255 x 290','315 x 205 x 270','1200'
    '19 litre XL Really Useful Box (LP small)','395 x 255 x 330','315 x 205 x 310','1300'
    '20 litre Really Useful Box','710 x 440 x 120','610 x 370 x 090','2150'
    '21 litre Really Useful Box','450 x 350 x 200','370 x 300 x 175','1330'
    '21 litre XL Really Useful Box','450 x 350 x 230','370 x 300 x 210','1520'
    '22 litre Really Useful Box','820 x 255 x 155','750 x 210 x 140','1410'
    '24 litre Really Useful Box','465 x 270 x 290','380 x 230 x 270','1450'
    '24.5 litre Really Useful Box','600 x 400 x 155','530 x 350 x 130','1576'
    '25 litre Really Useful Box','395 x 255 x 360','330 x 220 x 345','1362'
    '32 litre folding Really Useful Box','470 x 342 x 234','443 x 315 x 225','945'
    '33 litre Really Useful Box','710 x 440 x 165','605 x 370 x 145','2233'
    '35 litre Really Useful Box','480 x 390 x 310','370 x 310 x 280','2100'
    '35 litre XL Really Useful Box (LP medium)','480 x 390 x 345','370 x 310 x 310','2300'
    '42 litre Really Useful Box','520 x 440 x 310','405 x 365 x 280','2700'
    '45 litre folding Really Useful Box','390 x 570 x 285','360 x 520 x 280','1490'
    '48 litre Really Useful Box','600 x 400 x 315','530 x 325 x 289','2575'
    '48 litre XL Really Useful Box','600 x 400 x 350','530 x 325 x 330','2900'
    '50 litre Really Useful Box','710 x 440 x 230','605 x 370 x 200','2700'
    '64 litre Really Useful Box','710 x 440 x 310','605 x 370 x 280','3100'
    '64 litre open front Really Useful Box','710 x 440 x 310','605 x 370 x 280','3100'
    '70 litre Really Useful Box','810 x 620 x 225','710 x 545 x 190','4055'
    '77 litre Really Useful Box','1201 x 270 x 355','1135 x 232 x 340','3935'
    '84 litre Really Useful Box','710 x 440 x 380','605 x 370 x 355','3300'
    '145 litre Really Useful Box','810 x 620 x 430','710 x 545 x 400','4164'};

end
