function D=suit_ROI_summarize(images,varargin);
% function varargout=suit_ROI_summarize(images,varargin);
% Uses cerebellar ROIs from an atlas image to summarize cerebellar data 
%__________________________________________________________________________
% INPUT:
%   images: character or cell array of images to summarize
%__________________________________________________________________________
% OPTIONS:
%   'atlas',atlas_image: Atlas image (defaults_suit.summarize.atlas) 
%   'regionname', names:   Cell array of region names 
%   'outfilename',name : Name under which the structure will be saved as txt-file
%   'stats',{'nanmean','max',...}: Name of statistical function to be
%                         performed on data within lobuli. A field of the
%                         same name will be added to the structure and
%                         text-file.
%__________________________________________________________________________
% OUTPUT:
%   Structure with fields
%       D.image:    Image number
%       D.region:   Region number
%       D.regionname: Region name (if given)
%       D.size:     Size of Region in mm3
%       D.nanmean:  Mean excluding nans
%       D.max:      Maximal value
%       D.....:     Addtional fields for any additional statistics
%__________________________________________________________________________
% Example:
%   
% -------------------------------------------------------------------------
% Copyright (C) 2010 
% Joern Diedrichsen (joern.diedrichsen@googlemail.com)
% Diedrichsen, J. (2006). A spatially unbiased atlas template of the human
% cerebellum. Neuroimage.

% v.2.3: First released version
% v.2.4: Compatibility update for SPM8 
%       + bugfix for missing spmj_affine_transform 
% v.2.5: made compatible with SPM8-jobmanager (26/08/2010) 
% v.3.3: Introduced suit_ROI_summarize
% v.3.4: 

global defaults_suit;
if isempty(defaults_suit)
    suit_defaults;
end


atlas=defaults_suit.summarize.atlas{1};

stats={'nanmean'};
regionname={}; 
outfilename=[];

SCCSid   = '3.7';
SPMid    = spm('FnBanner',mfilename,SCCSid);


vararginoptions(varargin,{'outfilename','atlas','stats','regionname'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin<1 || isempty(images))
    images=spm_select(inf,'image','Select images to do statistics on');
end;

if (~isstruct(atlas))
    if (~exist(atlas,'file'))
        error(sprintf('Atlas file: %s not found. \nYou may have to download github/DiedrichsenLab/cerebellar_atlases, \n or set the location of the atlas directory in suit_defaults.m',atlas));
    end
    Vatlas=spm_vol(atlas);
else
    Vatlas=atlas;
end;
if (ischar(images))
    V=spm_vol(images);
elseif (iscell(images))
    V=spm_vol(char(images));
elseif (isstruct(images))
    V=images;
else 
    error('Images must be string array, cell arrary, or structure(s)');  
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get (x,y,z) - get coordinates for all locations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=spm_read_vols(Vatlas);
numregions = max(X(:)); 
for r=1:numregions
    Xindx=find(X==r);
    [x{r},y{r},z{r}]=ind2sub(size(X),Xindx);
    [x{r},y{r},z{r}]=spmj_affine_transform(x{r},y{r},z{r},Vatlas.mat);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now loop over images and get data from all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D=[];
for i=1:length(V)
    for r=1:numregions
        [a,b,c]=spmj_affine_transform(x{r},y{r},z{r},inv(V(i).mat));
        Data=spm_sample_vol(V(i),a,b,c,1);
        T.image=i;
        T.region=r;
        if (~isempty(regionname) && r<=length(regioname)); 
            T.regionname={regionname{r}};
        else 
            T.regionname={sprintf('region %d',r)}; 
        end; 
        T.size=length(x{r})*abs(det(Vatlas.mat));
        if T.size<1
            error('no voxels found, please check that you have the correct atlas file');
        else
            for s=1:length(stats)
                T.(stats{s})=eval([stats{s} '(Data)';]);
            end;
            D=addstruct(D,T);
        end;
    end;
end;

if (nargout==0) 
    if (isempty(outfilename))
        outfilename=uiputfile({'*.txt', 'Tab delimted text file'},'Save Table as');
    end;
    dsave(outfilename,D);
end; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dsave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dsave(filename,D)
%  function dsave(filename,D)
%       writes out a Data structure as a tab-delimted worksheet
%       filename if empty, function will promt user for filename
%       D: data-struct
%           each field is a column
% v1.0: 2002 Berkeley
% v1.2: Speed up code (9/21/05)
%  Joern Diedrichsen (jdiedric@jhu.edu)
if nargin==1
    D=filename;
    filename='';
end;
if (~isstruct(D))
    error('Cannot save matrix or empty structure');
end;
if isempty(filename)
    [F,P]=uiputfile('*.*','Save Cell Array as');
    filename = [P,F];
end
fid=fopen(filename,'wt');
if (fid==-1)
    error(sprintf('Error opening file %s\n',filename));
end;
names=fieldnames(D);
linefeed=sprintf('\n');
tab=sprintf('\t');

% put in variable names
numvar=length(names);
for v=1:numvar
    fprintf(fid,'%s',names{v});
    if(v<numvar)
        fprintf(fid,'\t');
    end;
    var_length(v)=size(getfield(D,names{v}),1);
end;
fprintf(fid,'%s',linefeed);

% Check if all variables have the same length
if (any(var_length~=var_length(1)))
    error('dsave: all variables must have the same length');
end;
% now get out all the variables
TEXT=[];
for v=1:numvar
    var=getfield(D,names{v});
    if (iscell(var))
        TEXT=[TEXT char(var)];
    elseif isnumeric(var) || islogical(var)
        TEXT=[TEXT num2str(var)];
    elseif ischar(var)
        TEXT=[TEXT var];
    else
        error('variable types have to be cell, numeric or character');
    end;
    if (v<numvar)
        TEXT=[TEXT ones(var_length(1),1).*tab];
    end;
end;
for l=1:size(TEXT,1)
    fprintf(fid,'%s\n',TEXT(l,:));
end;
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addstruct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D=addstruct(D,A,type,force)
% function D=addstruct(D,A,type,force)
% Adds the fields of structure A to the fields of structure D
% adds the field as a row or column depending on type
% last concatinates 2-d along columns, N-d along the last dimension
% type = 'row' / 'column' / 'last'
%   row: (DEFAULT) add as rows
%   column: add as columns
%   last: add as the last dimension
% 'force': Forces all fields to be added in same row/column
%          This means, if a field is not exisitent or shorter, it will
%          be padded with NaNs
%          If a field is no exisitent in the added structure, it will also
%          be padded with NaNs;
% Joern Diedrichsen
% v.1.1 09/18/05: added support for cell arrays
% v.1.2 12/14/05: taken out reference to header
% v.1.3 06/07/06: force option that includes NaN for
%               so far non-existing fields
% -------------------------------------------------------
names=fieldnames(A);
Dnames=fieldnames(struct(D));
if (nargin <3 | strcmp(type,'row'))
    dim=1;
elseif (strcmp(type,'column'))
    dim=2;
elseif (strcmp(type,'last'))
    dim=length(size(eval(['A.' names{1} ';'])));
else
    error('unknown option (row/column/last)');
end;

if (nargin >3 & strcmp(force,'force'))
    force=1;
    if (~isempty(Dnames))
        Dlength=size(getfield(D,Dnames{end}),dim);
    else
        Dlength=0;
    end;
else
    force=0;
end;

for (i=1:length(names))
    if (~isfield(D,names{i}) & ~force)
        eval(['D.' names{i} ' = A.' names{i} ';']);
    elseif (~isfield(D,names{i}) & force)
        if (strcmp(type,'row'))
            eval(['Alength= size(A.' names{i} ',2);']);
            eval(['D.' names{i} '= ones(Dlength,Alength).*NaN;']);
        elseif (strcmp(type,'column'))
            eval(['Alength= size(A.' names{i} ',1);']);
            eval(['D.' names{i} '= ones(Alength,Dlength).*NaN;']);
        else
            error ('force does only work with column / row');
        end;
        eval(['D.' names{i} '= cat(dim,D.' names{i} ' ,A.' names{i} ');']);
    else
        eval(['D.' names{i} '= cat(dim,D.' names{i} ' ,A.' names{i} ');']);
    end;
end;

if (force)
    for i=1:length(Dnames)
        if (isempty(strmatch(Dnames{i},names)))
            if (strcmp(type,'row'))
                eval(['cases= size(A.' names{1} ',1);']);
                eval(['width= size(D.' Dnames{i} ',2);']);
                eval(['D.' Dnames{i} '= cat(dim,D.' Dnames{i} ' ,ones(cases,width).*NaN);']);
            elseif (strcmp(type,'column'))
                eval(['cases= size(A.' names{1} ',2);']);
                eval(['leng= size(D.' Dnames{i} ',1);']);
                eval(['D.' Dnames{i} '= cat(dim,D.' Dnames{i} ' ,ones(leng,cases).*NaN);']);
            else
                error ('force does only work with column / row');
            end;
        end;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spmj_affine_transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y1,y2,y3] = spmj_affine_transform(x1,x2,x3,M)
% function [y1,y2,y3] = affine_transform(x1,x2,x3,M)
% -----------------------------------------------------------------
% Affine Transform for input stuff in any format (N-dim strcutures)
% -----------------------------------------------------------------
y1 = M(1,1)*x1 + M(1,2)*x2 + M(1,3)*x3 + M(1,4);
y2 = M(2,1)*x1 + M(2,2)*x2 + M(2,3)*x3 + M(2,4);
y3 = M(3,1)*x1 + M(3,2)*x2 + M(3,3)*x3 + M(3,4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nanmean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = nanmean(x)
    i=find(~isnan(x)); 
    y=sum(x(i))./length(i); 
    

