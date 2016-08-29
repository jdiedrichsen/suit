function varargout=suit_isolate(Source,varargin);
% function varargout=suit_isolate(Source,varargin);
% Isolation algorithm for the cerebellum and brainstem for subsequent
% normalization to a specialized template
% The Algorithm relies on
% a) A Bayesian prior in MNI space of where the cerebellum is
% b) A segmentation into issue types
% c) A seperation of the cortical and cerebellar white matter that can then
%   be used to classify the gray matter in terms of distance from white
%   matter
%
% For the computation, the alogrithm goes through the following steps:
% 1. Finding a affine alignment to whole brain template
% 2. Segment source image into white and gray matter, using John Ashburners
%       algorithm
% 3. Crop the image, so that it includes only the cerebellum, this also
%   brings cropped images in LPI orientation
% 4. seperate white matter into cortical and cerebellar white matter
% 5. Integrate the information sources
% 6. Repeat steps 1-5, but replace the whole brain template with a new
%       cerebellar specific template. The segmentation acquired in the last
%       step is used to find the alignment.
%__________________________________________________________________________
% INPUT:
%       Whole brain anatomical scan
%__________________________________________________________________________
% OUTPUT:
%       The cropped anatomical containing the cerebellum
%   c_<Source>_pcereb
%       The final output of the algorithm is a probability map that indicates the
%       posterior probability of a voxel belonging to the cerebellum
%   c_<Source>_pcereb_corr
%       The same map, thresholded at p=0.5 for hand correction
%   <Source>_seg1 - <Source>_seg3
%       Probabilities of each voxel to be Gray matter, White matter or CSF
%__________________________________________________________________________
% OPTIONS:
%       'crop',{'none','reslice','crop'}: How to crop images
%                       'crop' preserves the original image but requires LPI
%                       orientation
%       'bb',2x3 array:    Bounding box for cropping (defined in MNI)
%       'Ms2a',mat:     Tranformation matrix for affine alignement to atlas
%       'prior',filename: Prior probability image (p(c))
%       'keeptempfiles',0/1: Keep temporary files? default: no
%       'new_segmentation',0/1: use new segmentation routine in SPM5/8?
%__________________________________________________________________________
% EXAMPLE:
%       suit_isolate('<name>.nii'); % <name> = name of whole-brain
%       structural.
% ------------------------------------------------------------------------
% Copyright (C) 2010 
% Joern Diedrichsen (j.diedrichsen@ucl.ac.uk)
% Use of the software should be cited as:
% Diedrichsen, J. (2006). A spatially unbiased atlas template of the human
% cerebellum. Neuroimage.

% Version history 
%        2/01/06
% 1.2    16/06/06 compatibility for both SPM2 / SPM5
%        6/07/06 compatibility bug with matlab 6.5, bigger bounding box
%        19/07/06 additionally to cropping the script can now resample
%                 into bounding box
% 1.2.1  31/08/06 Use for second affine registration the already found
%                 affreg to the MNI template
% 1.2.2  9/5/07   Add cerebral range and cerebellar range to vararinoptions
% 2.0    21/7/07  One-file Nifti (*.nii) support:
%                   If the Source file is a .img pair, then all will be
%                   saved as *.hdr/*.img pair
%                   If the Source file is a .nii, then all will be
%                   saved as *.nii files
%                 -Change in white-matter classifaction algorithm to a
%                   volume-based clustering
%                 - Full parameter options to be handed into varargin
%                 - Acceleration of Smoothing and whitematter range
%                   calculation
%                 - Support for any image orientation, not only LPI
%                 - reslicing now is done while trying to preserved
%                   voxel-sizes
% 2.1     20/3/08 No changes to suit_isolate, but updated with new toolbox version
% 2.2     15/4/08 - in SPM5 now uses the segmentation + normalization rather
%                   than the standard segmentation
%                 - uses fullfile to be platform independent  
% 2.4     12/01/10 - compatibility with SPM8 
% 2.5     24/08/10 - Job manager support in SPM8
% 2.5.2   27/12/10 - bug fixes with cropping
% 2.6     7/04/13  - update for SPM12b
% 2.7:    13/09/13 - Further compartibility with SPM 12b and improved 
%                     documentation. George Prichard (dgmprichard@gmail.com)
% 3.0:    27/02/15 - Compatibility with SPM12 
%                    Cerebellar flatmap added 
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% Set the defaults.
% ------------------------------------------------------------------------
% All defaults can be changed by calling the functions with
% suit_isolate(Source,var_name,value,.....);
crop='crop';
global defaults;
if (~isstruct(defaults))
    error('Start SPM to use SUIT_isolate');
end;

global defaults_suit; 
if (~isstruct(defaults_suit))
    suit_defaults;
end;

names=fieldnames(defaults_suit.isolate); 
for i=1:length(names)
    eval([names{i} '= defaults_suit.isolate.' names{i} ';']); 
end;

% Conditional probabilities of being inside or outside the range
pTnc=0.3709;                % Probability of white or gray matter, given it's not the cerebellum s
pRc=0.0959;pnRc=1-pRc;      % Probability of being inside cort. range, given cerebellum
pRnc=0.9431;pnRnc=1-pRnc;   % Probability of being inside cort. range, given ~cerebellum
pAc=0.5174;pnAc=1-pAc;      % Probability of being inside cereb. range, given cerebellum
pAnc=0.0139;pnAnc=1-pAnc;   % Probability of being inside cereb. range, given ~cerebellum

% Flags: Force to redo steps that have been accomplished?
force_align=0;
force_range=0;
force_segment=0;
new_segmentation=1; 

% Get the defaults from varargin
vararginoptions(varargin,{'Ms2a','bb','affreg','prior1','prior2','gray_prior','white_prior','csf_prior','iterations',...
    'template_dir','prior_dir','template1','template2','template_weight2',...
    'crop','keeptempfiles','banner','cerebral_range','cerebellar_range','affregsmooth2','white_cutoff','range_FWHM','new_segmentation'},...
    {'force_range','force_segment','force_align'});

% -----------------------------------------------------------------
% check spm-version
spmVer=spm('Ver');
SCCSid   = '3.0';
SPMid    = spm('FnBanner',mfilename,SCCSid);

% -----------------------------------------------------------------
% Get Source image and Atlas Image
% -----------------------------------------------------------------
if (nargin<1 | isempty(Source))
    Source=spm_select(1,'image','Get Source Image');
end;

[source_dir,Source,ext,num]=spm_fileparts(Source);

    
if (isempty(source_dir))
    source_dir=pwd;
end;

Sourcec=fullfile(source_dir,Source);
V_atlas=spm_vol(fullfile(template_dir,template1));
V_source=spm_vol([Sourcec ext]);

% determine voxelsize
[u,s,v]=svd(V_source.mat(1:3,1:3));
[dummy,indx]=max(abs(u)');
voxelsize=diag(s);
voxelsize=voxelsize(indx)';

% if (~strcmp(crop,'none') & (any(diag(V_source.mat)<0) |  V_source.mat(1,2)~=0 | V_source.mat(1,3)~=0 | V_source.mat(2,3)~=0))
%     fprintf('Source not in LPI orientation (oblique or mirrored)\nWill reslice instead of crop\n');
%     crop='reslice';
% end;

% -----------------------------------------------------------------
% Get alignment with the whole-brain template
% -----------------------------------------------------------------

if (~exist('Ms2a'))
    fprintf('Get alignment...');
    Mname=[Sourcec '_2atlas.mat'];
    if (exist(Mname,'file') & ~force_align)
        load(Mname);
    else
        Ms2a=spmj_get_affine_mapping(V_source,V_atlas,affreg);
    end;
    Ma2s=inv(Ms2a*V_source.mat)*V_atlas.mat;
    save(Mname,'Ms2a');
    fprintf('done\n');
end;
Ma2s=inv(Ms2a*V_source.mat)*V_atlas.mat;

% -----------------------------------------------------------------
% segment source volume in gray and white matter
% In SPM5 uses the normalization and segmentation algorithm
% -----------------------------------------------------------------
seg1_image=fullfile(source_dir,[Source,'_seg1',ext]); 
seg2_image=fullfile(source_dir,[Source,'_seg2',ext]); 
if ~exist(seg1_image) | ~exist(seg2_image) | force_segment
    fprintf('Segmentation...');
    if new_segmentation==0

        flags.estimate.priors = str2mat(...
            fullfile(prior_dir,gray_prior),...
            fullfile(prior_dir,white_prior),...
            fullfile(prior_dir,csf_prior));
        VO=spm_segment(V_source,Ms2a,flags);
        VO(1).fname=seg1_image;
        VO(2).fname=seg2_image;
        VR= spm_create_vol(VO(1));
        VR=spm_write_vol(VR,double(VO(1).dat).*VO(1).pinfo(1));
        VR= spm_create_vol(VO(2));
        spm_write_vol(VR,double(VO(2).dat).*VO(2).pinfo(1));
    elseif new_segmentation==1
        switch spmVer
            case {'SPM12b','SPM12'}
                opts=defaults.old.preproc;
            case {'SPM5','SPM8'}
                opts=defaults.preproc;
        end
        outp.GM=[0 0 1];
        outp.WM=[0 0 1];
        outp.CSF=[0 0 0];
        outp.biascor=1;
        outp.cleanup=0;
        
        if (isfield(opts,'tpm') && iscell(opts.tpm))
            opts.tpm=strvcat(opts.tpm{:});
        end; 
        res        = spm_preproc(V_source,opts);
        [sn,isn]   = spm_prep2sn(res);
        savefields(fullfile(source_dir,[Source '_seg_sn.mat']),sn);
        savefields(fullfile(source_dir,[Source '_seg_inv_sn.mat']),isn);
        spm_preproc_write(sn,outp);
        movefile(fullfile(source_dir,['c1' Source ext]),seg1_image);
        movefile(fullfile(source_dir,['c2' Source ext]),seg2_image);
        if (strcmp(ext,'.img'))
            movefile(fullfile(source_dir,['c1' Source '.hdr']),fullfile(source_dir,[Source,'_seg1.hdr']));
            movefile(fullfile(source_dir,['c2' Source '.hdr']),fullfile(source_dir,[Source,'_seg2.hdr']));
        end;
    end;
    fprintf('done\n');
end;

% -----------------------------------------------------------------
% crop the volumes to a standard bounding box
% 'crop': This just crops the image, preserving the orientation and
%   voxelsize of the original image. IMPORTANT: the image has to have LPI
%   orientation do do this
% 'reslice': To deal with oblique or mirrored orientations of the
%   original image, the cropping is performed by resampling into an image
%   that is in the standard cerebellar bounding box and would have a 1x1x1 mm
%   in MNI space (but not in inidividual space)
% -----------------------------------------------------------------
if (strcmp(crop,'crop'))
    fprintf('Cropping...');
    BOUND=[bb(1,1) bb(1,2) repmat(mean(bb(1,:)),1,4);...
        mean(bb(2,:)) mean(bb(2,:)) bb(2,1) bb(2,2) mean(bb(2,:)) mean(bb(2,:));...
        repmat(mean(bb(3,:)),1,4) bb(3,1) bb(3,2) ;...
        ones(1,6)];
    BOUND=round(Ma2s*inv(V_atlas.mat)*BOUND);
    cropped_bb.x=[min(BOUND(1,:)):max(BOUND(1,:))];
    cropped_bb.y=[min(BOUND(2,:)):max(BOUND(2,:))];
    cropped_bb.z=[min(BOUND(3,:)):max(BOUND(3,:))];
    % Crop the volumes around the cerebellum
    crop_image=fullfile(source_dir,['c_' Source ext]);
    crop_seg1=fullfile(source_dir,['c_' Source '_seg1' ext]);
    crop_seg2=fullfile(source_dir,['c_' Source '_seg2' ext]);
    try
        spmj_crop_vol([Sourcec ext],crop_image,cropped_bb.x,cropped_bb.y,cropped_bb.z);
        spmj_crop_vol(seg1_image,crop_seg1,cropped_bb.x,cropped_bb.y,cropped_bb.z);
        spmj_crop_vol(seg2_image,crop_seg2,cropped_bb.x,cropped_bb.y,cropped_bb.z);
    catch
        error(['Cropping failed.','Make sure that the original image is in LPI orientation.' ...
            'Otherwise call suit_isolate with option ''crop'',''reslice''']);
    end;
    fprintf('done\n');
    Source =['c_' Source];
    Sourcec =fullfile(source_dir,Source);
    seg1_image=crop_seg1;
    seg2_image=crop_seg2;
elseif (strcmp(crop,'reslice'))
    fprintf('Reslicing...');

    % determine transformation from bounding box into the original image
    % preserving the voxelsize of the image optimally
    x=[bb(1,1):voxelsize(1):bb(1,2)];
    y=[bb(2,1):voxelsize(1):bb(2,2)];
    z=[bb(3,1):voxelsize(1):bb(3,2)];
    dim=[length(x) length(y) length(z)];
    mat=spm_matrix([bb(1,1)-voxelsize(1),bb(2,1)-voxelsize(2),bb(3,1)-voxelsize(3),0,0,0,...
        voxelsize(1),voxelsize(2),voxelsize(3)]);
    mat=inv(Ms2a)*mat;

    % Approximate the affine transform without shear
    R         = mat(1:3,1:3);
    vx        = sqrt(sum(R.^2));
    R         = R * diag(1./vx);
    [U,S,V] = svd(R);
    R       = U*V';
    Rot=eye(4);Rot(1:3,1:3)=R;
    Trans=eye(4);Trans(1:3,4)=mat(1:3,4);
    Scale=diag([vx 1]);
    mat=Trans*Rot*Scale;

    % Crop the volumes around the cerebellum
    crop_image=fullfile(source_dir,['c_' Source ext]);
    crop_seg1=fullfile(source_dir,['c_' Source '_seg1' ext]);
    crop_seg2=fullfile(source_dir,['c_' Source '_seg2' ext]);
    spmj_reslice_vol([Sourcec ext],dim,mat,crop_image);
    spmj_reslice_vol(seg1_image,dim,mat,crop_seg1);
    spmj_reslice_vol(seg2_image,dim,mat,crop_seg2);
    fprintf('done\n');
    Source =['c_' Source];
    Sourcec =fullfile(source_dir,Source);
    seg1_image=crop_seg1;
    seg2_image=crop_seg2;
end;

% -----------------------------------------------------------------
% read relevant volumnes
% -----------------------------------------------------------------
% Read in anatomy and segmenttation volumes
Vana=spm_vol([Sourcec ext]);
ANA=spm_read_vols(Vana);
Vs1=spm_vol(seg1_image);
S1=spm_read_vols(Vs1);
Vs2=spm_vol(seg2_image);
S2=spm_read_vols(Vs2);
xdim=Vana.dim(1);
ydim=Vana.dim(2);
zdim=Vana.dim(3);
if (any((Vana.mat(:) - Vs1.mat(:))~=0) | any((Vana.mat(:) - Vs2.mat(:))~=0))
    error('anatomical and segmentation have to be in the same space');
end;

% -----------------------------------------------------------------
% Get the prior
% -----------------------------------------------------------------
if (~exist('prior') | isempty(prior))
    prior=fullfile(prior_dir,prior1);
end;
V_prior=spm_vol(prior);
PRIOR=zeros([xdim ydim zdim]);
M=inv(V_prior.mat)*Ms2a*Vana.mat;
[X,Y]=meshgrid(1:xdim,1:ydim);
for i=1:zdim
    Z=ones(ydim,xdim)*i;
    [x,y,z]=affine_transform(X,Y,Z,M);
    PRIOR(:,:,i)=spm_sample_vol(V_prior,x,y,z,1)';
end;

% -----------------------------------------------------------------
% make the information on range of cortical and cerebellar white matter
% -----------------------------------------------------------------
fprintf('White matter range.');
if (~exist([Sourcec '_inrange_in' ext]) | ~exist([ Sourcec '_inrange_out' ext]) | force_range)

    % Use smoothing of Prior within whitematter body only to seperate
    % Cerebellar from cortical white matter
    P_Outside=S2;
    P_Outside(P_Outside<white_cutoff)=0;
    WM=double(P_Outside>0);
    PRIORM=PRIOR.*WM;
    WM1=zeros(size(WM));
    PRIORM1=zeros(size(WM));
    prior_Sigma  = prior_FWHM/sqrt(8*log(2))./voxelsize;				% FWHM -> Gaussian parameter
    spm_smooth(WM,WM1,prior_Sigma);fprintf('.');
    spm_smooth(PRIORM,PRIORM1,prior_Sigma);fprintf('.');
    WM1(WM==0)=NaN;
    P_Inside=P_Outside.*((PRIORM./WM1).*WM>0.5);
    clear WM1 WM PRIORM1 PRIORM;

    % Reduce cerebellar white matter body to one object
    indx=find(P_Inside>0);
    [x,y,z]=ind2sub([xdim ydim zdim],indx);
    A=spm_clusters([x y z]');
    num_clusters=max(A);
    for c=1:num_clusters
        n(c)=length(find(A==c));
    end;
    [dummy,maxcl]=max(n);
    P_Inside(indx(A~=maxcl))=0;
    save_vol(P_Inside,[Sourcec '_inside' ext],Vs1);
    fprintf('.');

    % Do white matter range outside
    P_Outside=P_Outside-P_Inside;
    indx=find(P_Outside>0);
    [x,y,z]=ind2sub([xdim ydim zdim],indx);
    A=spm_clusters([x y z]');
    num_clusters=max(A);
    for c=1:num_clusters
        i=find(A==c);
        if length(i)<150
            P_Outside(indx(i))=0;
        end;
    end;
    save_vol(P_Outside,[Sourcec '_outside' ext],Vs1);
    save_vol(P_Inside,[Sourcec '_inside' ext],Vs1);

    % Now measure cerebral_range mm away from cortical white matter
    range_Sigma  = range_FWHM/sqrt(8*log(2))./voxelsize;				% FWHM -> Gaussian parameter
    kernelsize=round(cerebral_range./voxelsize.*2);
    P_OutsideS=zeros(size(P_Outside));
    for i=1:3
        x  = round(6*range_Sigma(i)); x = [-x:x];
        x  = exp(-(x).^2/(2*(range_Sigma(i)).^2));
        x=x./sum(x);
        kernel{i}=ones(1,kernelsize(i));
        kernel{i}=conv(kernel{i},x);
        lengthkernel(i)=length(kernel{i});
    end;
    spm_conv_vol(P_Outside,P_OutsideS,kernel{i},kernel{2},kernel{3},-[lengthkernel/2]);
    P_OutsideS(P_OutsideS>1)=1;
    V_OutsideS=save_vol(P_OutsideS,[Sourcec '_inrange_out' ext],Vs1);
    fprintf('.');

    % Or 5mm away from the cerebellar white matter
    kernelsize=round(cerebellar_range./voxelsize.*2);
    P_InsideS=zeros(size(P_Inside));
    for i=1:3
        x  = round(6*range_Sigma(i)); x = [-x:x];
        x  = exp(-(x).^2/(2*(range_Sigma(i)).^2));
        x=x./sum(x);
        kernel{i}=ones(1,kernelsize(i));
        kernel{i}=conv(kernel{i},x);
        lengthkernel(i)=length(kernel{i});
    end;
    spm_conv_vol(P_Inside,P_InsideS,kernel{i},kernel{2},kernel{3},-[lengthkernel/2]);
    P_InsideS(P_InsideS>1)=1;
    V_InsideS=save_vol(P_InsideS,[Sourcec '_inrange_in' ext],Vs1);
    clear P_Outside P_OutsideS P_Inside P_InsideS;
end;
V_inrange_in=spm_vol([Sourcec '_inrange_in' ext]);
P_inrange_in=spm_read_vols(V_inrange_in);
V_inrange_out=spm_vol([Sourcec '_inrange_out' ext]);
P_inrange_out=spm_read_vols(V_inrange_out);
fprintf('done\n');

% -----------------------------------------------------------------
% Use Bayesian integration on all of this
% -----------------------------------------------------------------
% Probablity of p(C|cereb & Gray)
fprintf('Integrate...');
P_post0c=(S2+S1).*PRIOR./(PRIOR+(1-PRIOR).*pTnc);

% Do the integration of cerebral white matter range
PcgR=P_post0c*pRc./(P_post0c.*pRc+(1-P_post0c).*pRnc);
PcgnR=P_post0c*pnRc./(P_post0c.*pnRc+(1-P_post0c).*pnRnc);
P_post1=(P_inrange_out.*PcgR+(1-P_inrange_out).*PcgnR);

% Do the integration of cerebellar white matter range
PcgA=P_post1*pAc./(P_post1.*pAc+(1-P_post1).*pAnc);
PcgnA=P_post1*pnAc./(P_post1.*pnAc+(1-P_post1).*pnAnc);
P_post1c=(P_inrange_in.*PcgA+(1-P_inrange_in).*PcgnA);

sP_post1c=zeros(size(P_post1c));
spm_smooth(P_post1c,sP_post1c,[3 3 3]);
V_post=save_vol(sP_post1c,[Sourcec '_pcereb' ext],Vs1);
clear P_post0c PcgR PcgnR P_post1 PcgA PcgnA P_post1c sP_post1c S1 S2;
fprintf('done\n');

% -----------------------------------------------------------------
% Redo the alignment to the new cerebellar template and redo the whole thing
% Important: Start with the alignment to MNI space, to be reasonably close.
% Bug found thanks to Jill O'Reily
% -----------------------------------------------------------------
if (iterations>1)
    fprintf('New Alignment...');
    Ms2MNI=Ms2a;
    Vana.mat=Ms2MNI*Vana.mat;  % Premultiply transform to Atlas
    V_post.mat=Ms2MNI*V_post.mat;
    affreg2.smosrc = affreg.smooth2; %
    affreg2.regtype = 'subj';
    affreg2.WF = V_post;
    affreg2.WG = spm_vol(fullfile(template_dir,template_weight2));
    affreg2.sep = affreg2.smosrc/2;
    affreg2.globnorm=0;
    affreg2.debug=0;
    V_atlas=spm_vol(fullfile(template_dir,template2));
    % smooth and scale source
    VF = spm_smoothto8bit(Vana,affreg2.smosrc);
    V_atlas.pinfo(1:2,:) = V_atlas.pinfo(1:2,:)/spm_global(V_atlas);
    VF.pinfo(1:2,:) = VF.pinfo(1:2,:)/spm_global(VF);
    spm_plot_convergence('Init','Affine Registration','Mean squared difference','Iteration');
    M         = eye(4);
    Ms2a   = spm_affreg(V_atlas, VF,affreg2, M);
    spm_plot_convergence('Clear');

    Ms2a=Ms2a*Ms2MNI;       % Postmultipy the MNIspace transform
    fprintf('done\n');
    suit_isolate([Sourcec ext],'Ms2a',Ms2a,'prior',[prior_dir '/' prior2],...
        'iterations',iterations-1,'crop','none','force_range',...
        'keeptempfiles',keeptempfiles,'banner',0);
end;

if (keeptempfiles==1)
    save_vol(PRIOR,[Sourcec 'prior' sprintf('%d',iterations) ext],Vs1);
end;

if (iterations==1)
    P_final=spm_read_vols(V_post);
    V_post_corr=save_vol(P_final>0.5,[Sourcec '_pcereb_corr' ext],Vs1);
end;
% -----------------------------------------------------------------
% remove tempfile
% -----------------------------------------------------------------
if (keeptempfiles==0)
    rm_imgfile(['temp'],ext);
    rm_imgfile([Sourcec '_seg1'],ext);
    rm_imgfile([Sourcec '_seg2'],ext);
    rm_imgfile([Sourcec '_pcerebu'],ext);
    rm_imgfile([Sourcec '_inside'],ext);
    rm_imgfile([Sourcec '_outside'],ext);
    rm_imgfile([Sourcec '_inrange_in'],ext);
    rm_imgfile([Sourcec '_inrange_out'],ext);
end;

function rm_imgfile(name,ext);
if (exist([name ext],'file'))
    delete([name ext]);
end;
if (strcmp(ext,'.img'))
    if (exist([name '.hdr'],'file'))
        delete([name '.hdr']);
    end;
    if (exist([name '.mat'],'file'))
        delete([name '.mat']);
    end;
end;

% -----------------------------------------------------------------
% Affine Transform
% -----------------------------------------------------------------
function [y1,y2,y3] = affine_transform(x1,x2,x3,M)
y1 = M(1,1)*x1 + M(1,2)*x2 + M(1,3)*x3 + M(1,4);
y2 = M(2,1)*x1 + M(2,2)*x2 + M(2,3)*x3 + M(2,4);
y3 = M(3,1)*x1 + M(3,2)*x2 + M(3,3)*x3 + M(3,4);
return;

% -----------------------------------------------------------------
% save volume
% -----------------------------------------------------------------
function VR=save_vol(DATA,name,V);
VR    = struct(		'fname',name,...
    'dim',		[V.dim(1:3)],...
    'dt',   [spm_type('float32') 0],...
    'mat',		V.mat,...
    'pinfo',	[1 0 0]',...
    'descrip',	'temp');
VR    = spm_create_vol(VR);
spm_write_vol(VR,DATA);


% -----------------------------------------------------------------
% Function save fields
%------------------------------------------------------------------------
function savefields(fnam,p)
if length(p)>1, error('Can''t save fields.'); end;
fn = fieldnames(p);
if numel(fn)==0, return; end;
for i=1:length(fn),
    eval([fn{i} '= p.' fn{i} ';']);
end;
switch spm('Ver')
    case 'SPM12b'
        if spm_check_version('matlab','7') >= 0
            save(fnam,'-V6',fn{:});
        else
            save(fnam,fn{:});
        end;
    case {'SPM5','SPM8'}
        if spm_matlab_version_chk('7') >= 0
            save(fnam,'-V6',fn{:});
        else
            save(fnam,fn{:});
        end;
end
return;