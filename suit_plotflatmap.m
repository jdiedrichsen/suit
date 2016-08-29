function varargout=suit_plotflatmap(data,varargin)
% suit_plotflatmap('overlay',file,thresholds);
% Visualised cerebellar cortical acitivty on a flatmap in a matlab window
% INPUT:
%    data                Vector(s) to be plotted. If overlay is NaN,
%                        underlay is shown. Needs to be a 28935x1 vector 
% VARARGIN:
%  'underlay',file        Specific a metric or surface_shape file that
%                         dictates the shape of the grey underlay
%  'underscale',[min max] Color scale to determine the value to color mapping 
%                         for the underlay
%  'undermap',map         Color map for the underlay,a N x 3 Matrix specifying RGB values [0..1]
%  'type':                Data type to be shown, either
%                           func:  Functional data
%                           label: discrete labels
%  'cmap',C             Color map for overlap,a N x 3 Matrix specifying RGB values [0..1]
%  'cscale',[min max]   Color scale: determines the mapping of values to
%                        color map
%  'border',borderfile  Specifies a borderfile to plot. Use [] for no borders
%  'bordersize',pt      Border point size in pt (default 8) 
%  'threshold',val      Shows only values above a certain threshold
%  'xlims',[xmin xmax]  Limits of the x-coordinate
%  'ylims',[ymin ymax]  Limits of the y-coordinate
%
% (c) j.diedrichsen@ucl.ac.uk, 2014
% EXAMPLES:
% 1. Plot a functional volume at a certain treshold + and color scale
%   Data = suit_map2surf('name.nii','space','SUIT');
%   suit_plotflatmap(Data,'threshold',0.01,'cscale',[0.01 0.5]);
% 2. Plot flatmap overlayed with the lobules
%   Data = suit_map2surf('Cerebellum-lobules.paint','stats',@(x)nanmedian(x,2));
%   suit_plotflatmap(Data,'interpol',0,'cmap','Cerebellum-lobules.color.txt');
% 3. Plot flatmap with the 17 resting state networks by Buckner
%   Data = suit_map2surf('Buckner_17Networks.paint');
%   suit_plotflatmap('overlay',Data,'cmap','Buckner_17Networks.color.txt');
% -----------------------------------------------------------------

% -----------------------------------------------------------------
%   Set Default values and deal with user options
% -----------------------------------------------------------------
global defaults;
if (isempty(defaults))
    spm('Defaults','fmri');
    global defaults;
end; 
flat_dir   = [];                    %
coord      = 'FLAT.coord.gii';      % Coordinate file for flat map
topo       = 'CUT.topo.gii';        % Topology file for flat map
underlay   = 'SUIT.shape.gii';      % File determining colring of underlay
underscale = [-1 0.5];              % Color scale [min max] for the overlay
undermap   = gray;                    % Color map for underlay
type       = 'func';                 %Either funtional activation or label
threshold  = [];                    % Threshold for functional overlay
cscale     = [];                    % Color scale [min max] for the overlay
cmap       =  colormap;             % Use current colormap by default
border     = 'fissures_flat.mat';   % File containing fissure information
bordersize = 8;                     % Size of the border points 
xlims      =  [-100 100];           % X-limits for Figure
ylims      =  [-100 100];           % Y-limits for Figure
vararginoptions(varargin,{'coord','topo','underlay','underscale','undermap',...
    'type','threshold','cscale','cmap','border','bordersize','xlims','ylims','flat_dir'});

SCCSid   = '3.1';
SPMid    = spm('FnBanner',mfilename,SCCSid);

% -----------------------------------------------------------------
%  Determine directory for the flatmap files
% -----------------------------------------------------------------
spmVer=spm('Ver');
if ~strcmp(spmVer,'SPM12')
    error('plot flatmap requires SPM12 (for gifti support)');
end;

if (isempty(flat_dir))
    flat_dir=fullfile(defaults.tbx.dir{1},'suit','flatmap');
end;

% -----------------------------------------------------------------
%   Load the surface and determine X,Y coordinates for all tiles
% -----------------------------------------------------------------
C=gifti(fullfile(flat_dir,coord));
T=gifti(fullfile(flat_dir,topo));
P=size(C.vertices,1);
X(1,:)=C.vertices(T.faces(:,1),1);
X(2,:)=C.vertices(T.faces(:,2),1);
X(3,:)=C.vertices(T.faces(:,3),1);
Y(1,:)=C.vertices(T.faces(:,1),2);
Y(2,:)=C.vertices(T.faces(:,2),2);
Y(3,:)=C.vertices(T.faces(:,3),2);

% Find all tiles that have a single vertex (or more) in the image
k=find(any(X>xlims(1) & X<xlims(2),1) & any(Y>ylims(1) & Y<ylims(2),1));
X=X(:,k);
Y=Y(:,k);

% -----------------------------------------------------------------
%  Determine the underlay and assign color
% -----------------------------------------------------------------
UN=gifti(fullfile(flat_dir,underlay));

% Determine the shading of the faces by the vertices and scale the color assignment
d=[UN.cdata(T.faces(k(:),1),1) UN.cdata(T.faces(k(:),2),1) UN.cdata(T.faces(k(:),3),1)]';
M=size(undermap,1);
dindx=ceil((d-underscale(1))/(underscale(2)-underscale(1))*M);
dindx(dindx<1)=1;
dindx(dindx>M)=M;

% Now assign the RGB color to each face
for i=1:3
    for j=1:size(dindx,1)
        COL(j,:,i) = undermap(dindx(j,:),i)';
    end;
end;

% -----------------------------------------------------------------
%  determine the overlay and assign color
% -----------------------------------------------------------------

% If input data is empty, make vector of NaNs; 
if (isempty(data))
    data=nan(P,1); 
end; 

% Check input data
if (size(data,1))~=P 
    error(sprintf('Input data must be a numVert (%d) x 1 vector',P)); 
end; 

% Determine colormap 
if (ischar(cmap))
    CM=load(cmap);
    cmap=[CM(:,2) CM(:,3) CM(:,4)]/255;
end;
colormax = size(cmap,1);

% Determine the color of the overlay
switch (type)
    case 'label'
        % if plotting labels, use the numerical vlaues themselves
        dindx=[data(T.faces(k(:),1),1) data(T.faces(k(:),2),1) data(T.faces(k(:),3),1)]';
    case 'func'
        % Otherwise take the mean value
        d=[data(T.faces(k(:),1),1) data(T.faces(k(:),2),1) data(T.faces(k(:),3),1)]';
        if (isempty(cscale) || any(isnan(cscale)))
            cscale=[min(d(:)) max(d(:))];
        end;
        dindx=ceil((d-cscale(1))/(cscale(2)-cscale(1))*colormax);
        dindx(dindx<=0)=1;
        dindx(dindx>colormax)=colormax;
    otherwise
        error('unknown data type');
end;


% Assign the color for each tile from the color map
for i=1:3
    for j=1:size(dindx,1)
        % Apply map thresholding
        if (isempty(threshold) || isnan(threshold))
            indx = find(dindx(j,:)>0 & dindx(j,:)<=colormax);
        else
            indx = find(d(j,:)>threshold);
        end;
        COL(j,indx,i) = cmap(dindx(j,indx),i)';
    end;
end;

% -----------------------------------------------------------------
%   Draw the actual patches
% -----------------------------------------------------------------
drawmode=get(gca,'NextPlot'); 
if (~strcmp(drawmode,'add'))
    cla; 
end; 
p=patch(X,Y,COL);
set(gca,'XLim',xlims,'YLim',ylims);
set(p,'LineStyle','none');
set(p,'EdgeAlpha',0);
axis equal;

% -----------------------------------------------------------------
%   Plot the border
% -----------------------------------------------------------------
if (~isempty(border))
    hold on;
    load(fullfile(flat_dir,border));
    for i=1:length(Border)
        xB=Border(i).data(:,1);
        yB=Border(i).data(:,2);
        indx=find(xB>xlims(1) & xB<xlims(2) & yB>ylims(1) & yB<ylims(2));
        p=plot(xB(indx),yB(indx),'k.');
        set(p,'MarkerSize',bordersize);
    end;
    if (~strcmp(drawmode,'add'))
        hold off;
    end; 
end;
