function suit_normalize(name,varargin);
% function suit_normalize(name,optionname1,value1,...);
%   Normalizes the cerebellum of an individual to the SUIT-template
% defaults can be overwritten by passing the var_name and value
% Uses the cosine basis function approach introduced by J. Ashburner(1998)
% and makes heavy use of the normalisation fuctions of SPM
% Needs SPM 5-12 to be started
%
% Without values it will prompt for the
% 1. cropped individual image (name)
% 2. The mask that defines the cerebellum & brainstem (cerebellum_isolate)
%__________________________________________________________________________
% INPUT:
%       c_<name>.nii: cropped cerebellum from suit_isolate
%       c_<name>_pcereb_corr.nii: Hand corrected probability 
%__________________________________________________________________________
% OUTPUT:
%       Normalised image is written as wsuit_<name>.nii
%       And the normalization parameters are written as <name>_snc.mat
%__________________________________________________________________________
% OPTIONS:
%   'template',filename:        Specify a template other than SUIT
%   'template_weight',filename: Specify a weight image for the template
%   'source_weight',filename:   Specify a weight image for the source
%   'lesion_mask', filename:    Specifies a lesion mask for the source (1-weight)
%   'lesion_rim',2:             Specifies the width of the safety rim around the lesion
%   'reg',value:                Specifies the regularisation factor (default 1 - higher value, less deformation)
%   'param_postfix',string:     String that is attached to parameters _sn
%   'prefix',string:            String before resliced image
%   'mask',filename:            Name of corrected segmentation
%   'smooth_mask',value:       Gaussian kernel for smoothing of mask
%                               0 for no smoothing
%   'SUITs':                    sets template, template_weight and postfix
%                               for fitting to SUIT* template
%   'regularisation',reg
%   'estimate',struct:          flags for estimation
%   'write','struct:            flags for writing
%   'outfilename',name :        Name of written file<s>
%__________________________________________________________________________
% EXAMPLE:
%       suit_normalize('c_<name>.nii','mask', ...
%            'c_<name>_pcereb_corr.nii');
% ------------------------------------------------------------------------
% Copyright (C) 2010 
% Joern Diedrichsen (j.diedrichsen@ucl.ac.uk)
% Diedrichsen, J. (2006). A spatially unbiased atlas template of the human
% cerebellum. Neuroimage.

% v.1.2
% v.1.2.1: If images do not have enough overlap by themselves, use a
%          prealignment
% v.2.0: one-file nii support
% v.2.1: :Lesion mask support for masking of source image
%         Bounding box and Template are located 5mm higher in z-direction!
%         regularisation factor made flexible as an option
% v.2.2: Uses fullfile throughout for platform compatibility
% v.2.3: Masking does no longer assume same orientation of masking images
%        Lesion masking for lesion overlap analysis improved
% v.2.4: Compatibility with SPM8
% v.2.5: Make compatible with the job manager of SPM8
% v.2.5: Compatibility with SPM job manager (26/08/2010) 
% v.2.6: Compatibility with SPM 12b
% v.2.7: Further compartibility with SPM 12b and improved documentation 
%           George Prichard (dgmprichard@gmail.com)
% v.3.0: Compartibility with SPM 12
% ------------------------------------------------------------------------


global defaults;
template_weight=[];
global defaults_suit; 
if (~isstruct(defaults_suit))
    suit_defaults;
end;



names=fieldnames(defaults_suit.normalise); 
for i=1:length(names)
    eval([names{i} '= defaults_suit.normalise.' names{i} ';']); 
end;

% Single parameters for estimate
reg=estimate.reg;
nits=estimate.nits;

vararginoptions(varargin,{'mask','smooth_mask','template','template_weight','source_weight',...
    'lesion_mask','lesion_rim','reg','prefix','param','outfilename','param_postfix','regularisation','estimate',...
    'write','reg','nits'},{'SUITs'});

estimate.reg=reg;
estimate.nits=nits;

SCCSid   = '3.0';
SPMid    = spm('FnBanner',mfilename,SCCSid);


% ----------------------------------------------------------------------
% Get the image and mask if necessary
% ----------------------------------------------------------------------
spmVer=spm('Ver');
if (nargin<1 | isempty(name))
    name=spm_select(1,'image','Get Source Image');
end;

[source_dir,midname,postname,source_num]=spm_fileparts(name);
name=[midname postname];

if isempty(source_dir)
    source_dir=pwd;
end;

if (~exist('mask','var'))
    mask=spm_select([0 1],'image','Mask? (press Done for none)');
end;
[mask_dir,mask,postmask,mask_num]=spm_fileparts(mask);
mask=[mask postmask];
if isempty(mask_dir)
    mask_dir=pwd;
end;
Vsource=spm_vol(fullfile(source_dir,name));

% ----------------------------------------------------------------------
% Check images and transfer from cell to char
% ----------------------------------------------------------------------
if (iscell(template));template=char(template);end;
if (iscell(template_weight));template_weight=char(template_weight);end;
if (iscell(source_weight));source_weight=char(source_weight);end;
if (iscell(lesion_mask));lesion_mask=char(lesion_mask);end;
if (iscell(mask));mask=char(mask);end;


% ----------------------------------------------------------------------
% Mask the Source image
% ----------------------------------------------------------------------
if ~isempty(mask)
    mname=['m' midname];
    if (smooth_mask>0)
        spm_smooth(spm_vol(fullfile(mask_dir,mask)),fullfile(mask_dir,['s' mask]),[1 1 1]*smooth_mask);
        mask=['s' mask];
    end;
    VV(1)=Vsource;
    VV(2)=spm_vol(fullfile(mask_dir,mask));
    X=spm_read_vols(VV(2));
    scale=max(X(:));
    name=[mname postname];
    Vsource.fname=fullfile(source_dir,name);
    Vsource=spm_imcalc(VV,Vsource,sprintf('i1.*i2./%2.5f',scale),{0,0,1});
end;

% ----------------------------------------------------------------------
% Generate or load source_weight
% ----------------------------------------------------------------------
if (~isempty(source_weight))
    VWF=spm_vol(source_weight);
elseif (~isempty(lesion_mask))
    if (lesion_rim>0)
        lesion_mask=spmj_extender(lesion_mask,lesion_rim);
    end;
    VV=spm_vol(lesion_mask);
    % smooth lesion mask before
    % subtracting
    VWF = struct(	'fname',	'./test4.nii',...
        'dim',		Vsource.dim(1:3),...
        'dt',		[2 spm_platform('bigend')],...
        'mat',		Vsource.mat,...
        'descrip',	'spm - algebra');
    VWF=spm_imcalc(VV,VWF,'1-(i1>0)',{0,0,0});
else
    VWF=[];
end;

% ----------------------------------------------------------------------
% Do the Normalization: find parameters
% ----------------------------------------------------------------------
if (isempty(param))
    param=fullfile(source_dir,[name(1:end-4) param_postfix '.mat']);
    if(~isempty(template_weight))
        VWG=spm_vol(template_weight);
    else
        VWG=[];
    end;
    Vtemplate=spm_vol(template);
    % Check whether you have enough overlap between images
    cogsource=find_COG(Vsource);
    cogtemplate=find_COG(Vtemplate);
    d=cogtemplate-cogsource;
    if sqrt(d'*d)>30
        fprintf(['Images have little overlap, improving by prealignmen. \n'...
            'It is recommented to set (0,0,0) to the AC in the original image.\n']);
        Vsource.mat=spm_matrix(d')*Vsource.mat;
        estimate.graphics=0;        % Prevent graphics display
    end;
    spm_normalise(Vtemplate,Vsource,param,VWG,VWF,estimate);
    if sqrt(d'*d)>30
        load(param);
        VF.mat=spm_matrix(-d')*VF.mat;
        save(param,'Affine','Tr','VF','VG','flags');
        spm_normalise_disp(param,VF);
    end;
end;

% ----------------------------------------------------------------------
% Reslice the Source image into the Template
% ----------------------------------------------------------------------
V=spm_write_sn(fullfile(source_dir,name),param,write);
if (isempty(outfilename))
    outfilename=fullfile(source_dir,[prefix name]);
end;
V.fname=outfilename;
spm_write_vol(V,V.dat);

% ----------------------------------------------------------------------
% Find center of gravity for the image
% ----------------------------------------------------------------------
function coord=find_COG(V);
X=spm_read_vols(V);
x=squeeze(mean(mean(X>0,2),3));
y=squeeze(mean(mean(X>0,1),3))';
z=squeeze(mean(mean(X>0,1),2));
x=x.*[1:V.dim(1)]'./sum(x);
y=y.*[1:V.dim(2)]'./sum(y);
z=z.*[1:V.dim(3)]'./sum(z);
coord=[sum(x);sum(y);sum(z);1];
coord=V.mat*coord;
coord=coord(1:3);
