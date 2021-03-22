function [color,name] = alphabet(num_or_str,pie_chart)

% alphabet  a colormap with the alphabet of colors
%   
%   alphabet provides a colormap with 26 clearly distinguishable
%   colors corresponding to the letters of the alphabet, see
%   http://en.wikipedia.org/wiki/Help:Distinguishable_colors.
%
%   alphabet complies with standard syntax for colormaps. In addition,
%   [color,name] = alphabet(num_or_str,pie_chart) works as follows:
%
%   num_or_str is either the wanted number of colors or a string 
%   whose characters are the initials of the colors to be selected. If
%   pie_chart is nonzero, a pie chart with colors and names is drawn. 
%   The default values are num_or_str = 26 and pie_chart = 0.
%   color is a matrix of rgb colors, and name is a cell array of names.
%
% See also: colormap, alphabetdemo

% Ulrich Reif
% March 23, 2015

if nargin==0
  num_or_str = 26;
end

if ischar(num_or_str)
  I = lower(num_or_str) - 'a' + 1;
else
  I = mod(0:num_or_str-1,26) + 1;
end

C = [...
240,163,255;0,117,220;153,63,0;76,0,92;25,25,25;0,92,49;43,206,72;...
255,204,153;128,128,128;148,255,181;143,124,0;157,204,0;194,0,136;...
0,51,128;255,164,5;255,168,187;66,102,0;255,0,16;94,241,242;0,153,143;...
224,255,102;116,10,255;153,0,0;255,255,128;255,255,0;255,80,5]/255;
N = {...
'Amethyst','Blue','Caramel','Damson','Ebony','Forest','Green',...
'Honeydew','Iron','Jade','Khaki','Lime','Mallow',...
'Navy','Orpiment','Pink','Quagmire','Red','Sky','Turquoise',...
'Uranium','Violet','Wine','Xanthin','Yellow','Zinnia'};

color = C(I,:);
name  = N(I);

if nargin==2 && any(pie_chart(:))
  figure
  pie(0*I+1,name(end:-1:1))
  colormap(color(end:-1:1,:))
end
