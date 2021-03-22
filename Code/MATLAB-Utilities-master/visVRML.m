function handles = visVRML(varargin)
%VISVRML

%defaults
bkcolor = [1, 1, 1];
signs = [1, -1, -1];

%color strings
colnames = 'brgcmyk';
coldef = [     0     0     1
     1     0     0
     0     1     0
     0     1     1
     1     0     1
     1     1     0
     0.5   0.5   0.5];
 
%manages input arguments
if ~isstruct(varargin{1})
    errstr = ['1st input argument must be a VRML structure.'];
    error(errstr);
    
end
iplot = 0;
for ii = 1:nargin
    switch class(varargin{ii});
        case {'struct'}
            iplot = iplot + 1;
            wrl(iplot) = varargin{ii};
            colspec(iplot,:) = [0.5, 0.5, 0.5]; 
            FaceAlpha(iplot) = 1.0;
            
        case {'char'}
            iicol = strfind(colnames, varargin{ii});
            colspec(iplot,:) = coldef(iicol,:); 
            
        case {'double'}
            FaceAlpha(iplot,1) = varargin{ii};
            
        otherwise
            errstr = ['Input type ',class(varargin{ii}),' not supported.'];
            error(errstr)
    end
  
end
nplot = iplot; 

%manages the single plot case
if nplot == 1
    colspec = repmat(colspec,2,1);
    colvals = 0.5;
else
    colvals = [0:nplot-1]'./(nplot-1);
end

%sets up figure & axes
handles.fig = gcf;
clf(handles.fig);
set(handles.fig,'color',bkcolor);

%sets up axes
handles.ax = axes;
handles.light = light;
set(handles.ax,'color',bkcolor);
colormap(colspec);

%plots surfaces
hold on
for ii = 1:nplot
    %generates a colour vector
    npts = length( wrl(ii).pts );
    C = colvals(ii).*ones( npts,1);

    %plots surface
    ConnectionData = wrl(ii).knx;
    X = wrl(ii).pts(:,1) .* signs(1);
    Y = wrl(ii).pts(:,2) .* signs(2);
    Z = wrl(ii).pts(:,3) .* signs(3);
    handles.surfs(ii) = trisurf( ConnectionData, X, Y, Z, C);
    
    %surface specific formatting
    caxis([0 1]);
    set(handles.surfs(ii),'FaceAlpha',FaceAlpha(ii));
    
end
hold off

%formatting
axis equal
shading flat

%axis labels
labstr = {'X','Y','Z'};
for ii = 1:3
    if signs(ii) < 1 
        labstr{ii} = ['-',labstr{ii}];
    end
end
handles.xlab = xlabel(labstr{1});
handles.ylab = ylabel(labstr{2});
handles.zlab = zlabel(labstr{3});

%ends function VISVRML
end