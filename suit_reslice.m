function suit_reslice(source,param,varargin); 
% function suit_reslice(source,param,varname1,value1,...); 
%   Reslice a set of volumes (source) to the SUIT template 
%   With optinal masking by segmentation before 
%__________________________________________________________________________
% INPUT: 
%   source: image(s) to be resliced: character/cell array 
%   param: mat file with the deformations from normalization
%   (mc_<name>_snc.mat).
%__________________________________________________________________________
% OUTPUT:
%       is written as wc_<name> (or as <Outfilename> if given)
% If the image is also masked, the masked image is m<name>
% And the masked normalized image is written as wcm<name>
%__________________________________________________________________________
% OPTIONS:
%       mask: image to mask the image with (segementation), '' for none
%       smooth_mask: smoothing kernel for mask (default: 2mm) 
%       bb: Bounding box (default: [[-70 -100 -80];[70 -6 6]];)
%       vox: voxel size (default: [2 2 2])
%       prefix: prefix for reslices image (default: 'wc')
%       outfilename: If given it uses these as outfilenames: character
%       array / cell array 
%       interp: Interpolation: 0 nearest neighbourgh, 1: trilinear, etc.
%__________________________________________________________________________
% EXAMPLE:
%       suit_reslice('<source>','mc_<name>_snc.mat','mask', ... 
%           'c_<name>_pcereb_corr');
% -------------------------------------------------------------------------
% Copyright (C) 2010 
% Joern Diedrichsen (j.diedrichsen@ucl.ac.uk) 12/01/10
% Diedrichsen, J. (2006). A spatially unbiased atlas template of the human
% cerebellum. Neuroimage.

% v.2.0
% v.2.1 Compatibility issue with directory names on Mac-platform 
%       new Bounding box for corrected suit-template (5mm shift in
%       z-direction) 
% v.2.2 Compatibility with Linux resolved (Thanks to John Schlerf)
% v.2.3 write options as optinal commands 
% v.2.4 compatible with SPM8
% v.2.5 Job manager support in SPM8 (24/08/10)
% v.2.6 Support for 4D-niftis included 
% v.2.7 Support for SPM12b, removed support for SPM2
% v.3.0 Support for SPM12
% -------------------------------------------------------------------------

global defaults; 
global defaults_suit; 
if (~isstruct(defaults_suit))
    suit_defaults;
end;

interp=1; % Bug fix to avoid problems on some systems in which matlab appears to want to use it as a function
names=fieldnames(defaults_suit.reslice); 
for i=1:length(names)
    eval([names{i} '= defaults_suit.reslice.' names{i} ';']); 
end;
outfilename=[];

vararginoptions(varargin,{'mask','smooth_mask',...
    'prefix','bb','vox','outfilename','interp','preserve'});

flags.write.wrap = [0 0 0];
flags.write.vox=vox;
flags.write.bb=bb;
flags.write.interp=interp;
flags.write.preserve=preserve; 



SCCSid   = '3.0';
SPMid    = spm('FnBanner',mfilename,SCCSid);


% ----------------------------------------------------------------------
% Get the image and mask if necessary 
% ----------------------------------------------------------------------
spmVer=spm('Ver');
if (nargin<1 | isempty(source)) 
	source=spm_select([1 inf],'image','Get Source Image(s)');
end;


if (nargin<2 | isempty(param)) 
    param=spm_select(1,'mat','Deformation map (_snc.mat)');
end;

if (~exist('mask','var'))
    mask=spm_select([0 1],'image','Mask (press Done for none)');
end; 

%---------------------------------------------------------------------
% Make sure file names are character arrays instead of cells 
% ---------------------------------------------------------------------
if(iscell(source))
    source=char(source); 
end; 
if (iscell(param))
    param=char(param);
end; 
if (iscell(outfilename))
    outfilename=char(outfilename);
end; 
% ----------------------------------------------------------------------
% Mask the Source images
% ----------------------------------------------------------------------
if ~isempty(mask)
    if (iscell(mask))
        mask=mask{1}; 
    end; 
    [mask_dir,mask,postmask]=fileparts(mask);
    mask=[mask postmask];
    if isempty(mask_dir)
        mask_dir=pwd;
    end;

    if (smooth_mask>0) 
       spm_smooth(spm_vol(fullfile(mask_dir, mask)),fullfile(mask_dir,['s' mask]),[1 1 1]*smooth_mask);
       mask=['s' mask];
    end;
    VM=spm_vol(fullfile(mask_dir,mask)); 
    X=spm_read_vols(VM);
    scale=max(X(:));
    for i=1:size(source,1)
        [source_dir,fname,postfix,postnum]=spm_fileparts(source(i,:));
        if isempty(source_dir)
            source_dir=pwd;
        end;
        mname{i}=fullfile(source_dir,['m' fname postfix postnum]);       
        V(1)=spm_vol(source(i,:));
        V(2)=VM;
        Vo=V(1);
        Vo.fname=fullfile(source_dir,['m' fname postfix]);
        spm_imcalc(V,Vo,sprintf('i1.*i2./%f',scale),{0,0,4,[]});
    end;
    source=char(mname);
end;

% ----------------------------------------------------------------------
% Reslice the Source image into the Template 
% ----------------------------------------------------------------------

for i=1:size(source,1)
    [source_dir,filename,postname,num]=spm_fileparts(source(i,:));
    if (isempty(source_dir))
        source_dir=pwd;
    end;
    V=spm_write_sn(source(i,:),param,flags.write);
    if (isempty(outfilename))
        V.fname=fullfile(source_dir, [prefix filename postname]);
    else 
        V.fname=outfilename(i,:);
    end;
    spm_write_vol(V,V.dat);
end;
