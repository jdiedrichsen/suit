function suit = tbx_cfg_suit
% SPM Configuration file for the suit toolbox
% v.3.1 
% (c) Joern Diedrichsen (j.diedrichsen@ucl.ac.uk)  
% joern.diedrichsen@googlemail.com

% $Id: tbs_cfg_rwls.m 2222 2010-05-17 11:08:47Z joern $

rev = '$Rev: 3.1$';
addpath(fullfile(spm('dir'),'toolbox','suit'));

% global defaults_suit; 
% if (~isstruct(defaults_suit))
%     suit_defaults; 
% end; 

% ---------------------------------------------------------------------
% ISOLATE 
% source Source Image
% ---------------------------------------------------------------------
source         = cfg_files;
source.tag     = 'source';
source.name    = 'Source Image';
source.help    = {'Anatomical image to perform the isolation algorithm on.  The image will be cropped and a probabistic isolation will be performed. '};
source.filter = 'image';
source.ufilter = '.*';
source.num     = [1 1];

% ---------------------------------------------------------------------
% subj Subject
% ---------------------------------------------------------------------
subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {source};
subj.help    = {'Data for this subject.  The same parameters are used within subject.'};
% % ---------------------------------------------------------------------
% % esubjs Data
% % ---------------------------------------------------------------------
esubjs         = cfg_repeat;
esubjs.tag     = 'esubjs';
esubjs.name    = 'Data';
esubjs.help    = {'List of source images to be isolated.'};
esubjs.values  = {source};
esubjs.num     = [1 Inf];
% ---------------------------------------------------------------------
% bb Bounding box
% ---------------------------------------------------------------------
bb_crop         = cfg_entry;
bb_crop.tag     = 'bb';
bb_crop.name    = 'Bounding box';
bb_crop.help    = {'The bounding box (in mm) for cropping  the anatomical image. This bounding box is defined in atlas space and and will be translated into indidual space after affine registration.'};
bb_crop.strtype = 'e';
bb_crop.num     = [3 2];
bb_crop.def     = @(val)suit_get_defaults('isolate.bb', val{:});
% ---------------------------------------------------------------------
% bb Bounding box
% ---------------------------------------------------------------------
cerebellar_range         = cfg_entry;
cerebellar_range.tag     = 'cerebellar_range';
cerebellar_range.name    = 'Cerebellar gray matter thicknes';
cerebellar_range.help    = {'The value of cerebellar gray matter thickness (in mm) determines how far from the outer edge of the cerebellar white matter gray matter, voxels are surely classified as belonging to the cerebellum. This parameter determines the boundary between cerebral and cerebellar gray matter.'};
cerebellar_range.strtype = 'e';
cerebellar_range.num     = [1 1];
cerebellar_range.def     = @(val)suit_get_defaults('isolate.cerebellar_range', val{:});

% ---------------------------------------------------------------------
% bb Bounding box
% ---------------------------------------------------------------------
cerebral_range         = cfg_entry;
cerebral_range.tag     = 'cerebral_range';
cerebral_range.name    = 'Cerebral gray matter thicknes';
cerebral_range.help    = {'The value of cerebral gray matter thickness (in mm) determines how far from the outer edge of the cerebral white matter gray matter voxels are surely classified as belonging to the cerebrum. This parameter determines the boundary between cerebral and cerebellar gray matter. '};
cerebral_range.strtype = 'e';
cerebral_range.num     = [1 1];
cerebral_range.def     = @(val)suit_get_defaults('isolate.cerebral_range', val{:});

% ---------------------------------------------------------------------
% Isolate unit
% ---------------------------------------------------------------------
isolate         = cfg_exbranch;
isolate.tag     = 'isolate';
isolate.name    = 'Isolation';
isolate.val     = {esubjs bb_crop cerebral_range cerebellar_range} ;
isolate.help    = {
    'Isolation of the cerebellum from the surrounding tissue'
    'This is the original algorithm - for better performance use Segmentation & Isolation'
    }';
isolate.prog = @suit_run_isolate;

% ---------------------------------------------------------------------
% ISOLATE_SEG  
% Segmentation and isolation 
% ---------------------------------------------------------------------
sourceSeg         = cfg_files;
sourceSeg.tag     = 'source';
sourceSeg.name    = 'Subject Image(s)';
sourceSeg.help    = {'List of channels of one subject.'
                    'When loading multiple channels, select the T1 first'
                    'All images for different channels MUST be coregistered and resliced to the T1 (first channel)'};
sourceSeg.filter = 'image';
sourceSeg.ufilter = '.*';
sourceSeg.num     = [1 inf];

% ---------------------------------------------------------------------
% subj Subject
% ---------------------------------------------------------------------
subjSeg         = cfg_branch;
subjSeg.tag     = 'subj';
subjSeg.name    = 'Subject';
subjSeg.val     = {sourceSeg};
subjSeg.help    = {'Data for this subject. The same parameters are used within subject.'};

% % ---------------------------------------------------------------------
% % Channels
% % ---------------------------------------------------------------------
schan         = cfg_repeat;
schan.tag     = 'schan';
schan.name    = 'Inputs';
schan.help    = {'Anatomical image(s) to perform the isolation algorithm on.' 
                  'Multiple images (Channels e.g. T1 ,T2, PD) are allowed and will be used together for improved segmentation.'
                  'By default only one channel is needed (T1), if more channels will be used only select the images after the T1.'
                  'All images for different channels MUST be coregistered and resliced to the T1 (first channel)'
                  'To load multiple subjects click in "New: Subject Images(s)",then add the corresponding image(s)'};
schan.values  = {sourceSeg};
schan.num     = [1 Inf];

% ---------------------------------------------------------------------
% bb Bounding box
% ---------------------------------------------------------------------
bb_seg         = cfg_entry;
bb_seg.tag     = 'bb';
bb_seg.name    = 'Bounding box';
bb_seg.help    = {'The bounding box (in mm) for cropping  the anatomical image. This bounding box is defined in MNI space and and will be translated into indidual space after affine registration.'};
bb_seg.strtype = 'e';
bb_seg.num     = [3 2];
bb_seg.def     = @(val)suit_get_defaults('isolate_seg.bb', val{:});

% ---------------------------------------------------------------------
% Probability threshold
% ---------------------------------------------------------------------
threshold_prop         = cfg_entry;
threshold_prop.tag     = 'maskp';
threshold_prop.name    = 'Probability threshold for cerebellar isolation map';
threshold_prop.help    = {'The minimal probability (cerebellar white + gray matter' 
                          'to be included into the cerebellar isolation mask (default = 0.2)'};
threshold_prop.strtype = 'e';
threshold_prop.num     = [1 1];
threshold_prop.def     = @(val)suit_get_defaults('isolate_seg.maskp', val{:});

% ---------------------------------------------------------------------
% Keep temporal files
% ---------------------------------------------------------------------
keepFiles           = cfg_entry;
keepFiles.tag       = 'keeptempfiles';
keepFiles.name      = 'Keep temporary files';
keepFiles.help      = {'Set to 1 to keep temporary files created during the isolation procedure, including:'
                       'segmetation images of cerebellar gray/white matter, cortical gray/white matter,'
                       'and transformation matrix (default = 0)'};
keepFiles.strtype = 'e';
keepFiles.num     = [1 1];
keepFiles.def     = @(val)suit_get_defaults('isolate_seg.keeptempfiles', val{:});

% ---------------------------------------------------------------------
% Isolate_seg unit
% ---------------------------------------------------------------------
isolate_seg     = cfg_exbranch; 
isolate_seg.tag = 'isolate_seg'; 
isolate_seg.name = 'Segmentation & Isolation'; 
isolate_seg.val     = {schan bb_seg threshold_prop keepFiles};
isolate_seg.help    = {
    'Isolation of the cerebellum from the surrounding tissue'
    }';
isolate_seg.prog = @suit_isolate_seg;



% ---------------------------------------------------------------------
%
% NORMALIZATION 
% 
% source Source Image
% ---------------------------------------------------------------------
sourceN         = cfg_files;
sourceN.tag     = 'source';
sourceN.name    = 'Source Image';
sourceN.help    = {'The cropped image (c_*) that should be warped to match the template(s).',...
    'The result is a set of warps, which can be applied to this image, or any other image that is in register with it.',''};
sourceN.filter = 'image';
sourceN.ufilter = '.*';
sourceN.num     = [1 1];
% ---------------------------------------------------------------------
% Mask 
% ---------------------------------------------------------------------
maskN         = cfg_files;
maskN.tag     = 'mask';
maskN.name    = 'Cerebellar Mask';
maskN.help    = {'The isolation mask derived from the suit_isolate algorithm that removes all non-cerebellar structures.',...
    'The _pcereb_corr.nii should be viewed and possibly handcorrected before normalization.',''};
maskN.filter = 'image';
maskN.ufilter = '.*';
maskN.num     = [1 1];
% ---------------------------------------------------------------------
% lesion mask 
% ---------------------------------------------------------------------
lesionMaskN         = cfg_files;
lesionMaskN.tag     = 'lesion_mask';
lesionMaskN.name    = 'Lesion Mask';
lesionMaskN.val     = {''};
lesionMaskN.help    = {'Optional lesion mask for the cerebellum. If specified, the lesion will be extended by a rim of 3 mm from the normalization process.',''};
lesionMaskN.filter = 'image';
lesionMaskN.ufilter = '.*';
lesionMaskN.num     = [0 1];
% ---------------------------------------------------------------------
% subj Subject
% ---------------------------------------------------------------------
subjN         = cfg_branch;
subjN.tag     = 'subjN';
subjN.name    = 'Subject';
subjN.val     = {sourceN maskN lesionMaskN};
subjN.help    = {'Data for this subject.  The same parameters are used for all subjects.'};
% ---------------------------------------------------------------------
% esubjs Data
% ---------------------------------------------------------------------
subjsN         = cfg_repeat;
subjsN.tag     = 'subjsN';
subjsN.name    = 'Data';
subjsN.help    = {'List of subjects for the normalization step. Images of each subject should be warped differently.',''};
subjsN.values  = {subjN};
subjsN.num     = [1 Inf];
% ---------------------------------------------------------------------
% template Template Image
% ---------------------------------------------------------------------
template         = cfg_files;
template.tag     = 'template';
template.name    = 'Template Image';
template.help    = {'Specify a template image to match the source image with (default is SUIT.nii).'};
template.filter = 'image';
template.ufilter = '.*';
template.num     = [1 1]; 
template.def     = @(val)suit_get_defaults('normalise.template', val{:});
% ---------------------------------------------------------------------
% weight Template Weighting Image
% ---------------------------------------------------------------------
weight         = cfg_files;
weight.tag     = 'template_weight';
weight.name    = 'Template Weighting Image';
weight.help    = {
                  'Applies a weighting mask to the template(s) during the parameter estimation.',''};
weight.filter = 'image';
weight.ufilter = '.*';
weight.num     = [1 1];
weight.def     = @(val)suit_get_defaults('normalise.template_weight', val{:});
% ---------------------------------------------------------------------
% smosrc Source Image Smoothing
% ---------------------------------------------------------------------
prefixN         = cfg_entry;
prefixN.tag     = 'prefix';
prefixN.name    = 'Prefix anatomical';
prefixN.help    = {'Prefix with which the resliced anantomical image will be saved.',''};
prefixN.strtype = 's';
prefixN.num     = [1 Inf];
prefixN.def     = @(val)suit_get_defaults('normalise.prefix', val{:});
% ---------------------------------------------------------------------
% smosrc Source Image Smoothing
% ---------------------------------------------------------------------
parampostfix         = cfg_entry;
parampostfix.tag     = 'param_postfix';
parampostfix.name    = 'parameter postfix';
parampostfix.help    = {'Post fix for the deformation parameter file (default: _snc)',''};
parampostfix.strtype = 's';
parampostfix.num     = [1 Inf];
parampostfix.def     = @(val)suit_get_defaults('normalise.param_postfix', val{:});

% ---------------------------------------------------------------------
% smosrc Source Image Smoothing
% ---------------------------------------------------------------------
smooth_mask         = cfg_entry;
smooth_mask.tag     = 'smooth_mask';
smooth_mask.name    = 'Mask Image Smoothing';
smooth_mask.help    = {'Smoothing to apply to the mask image, before multiplying with the target image.',''};
smooth_mask.strtype = 'e';
smooth_mask.num     = [1 1];
smooth_mask.def     = @(val)suit_get_defaults('normalise.smooth_mask', val{:});
% ---------------------------------------------------------------------
% smosrc Source Image Smoothing
% ---------------------------------------------------------------------
smosrc         = cfg_entry;
smosrc.tag     = 'smosrc';
smosrc.name    = 'Source Image Smoothing';
smosrc.help    = {'Smoothing to apply to a copy of the source image.',...
    'The template and source images should have approximately the same smoothness.',...
     'By experience a value of 2mm gives good result while preserving anatomical detail.',''}; 
smosrc.strtype = 'e';
smosrc.num     = [1 1];
smosrc.def     = @(val)suit_get_defaults('normalise.estimate.smosrc', val{:});
% ---------------------------------------------------------------------
% smoref Template Image Smoothing
% ---------------------------------------------------------------------
smoref         = cfg_entry;
smoref.tag     = 'smoref';
smoref.name    = 'Template Image Smoothing';
smoref.help    = {'Smoothing to apply to a copy of the template image.','',...
    'The template and source images should have approximately the same smoothness.',...
    'On the standard SUIT image, no smoothing of the template is necessary',''};
smoref.strtype = 'e';
smoref.num     = [1 1];
smoref.def     = @(val)suit_get_defaults('normalise.estimate.smoref', val{:});
% ---------------------------------------------------------------------
% regtype Affine Regularisation
% ---------------------------------------------------------------------
regtype         = cfg_menu;
regtype.tag     = 'regtype';
regtype.name    = 'Affine Regularisation';
regtype.help    = {'Affine registration into a standard space can be made more robust by regularisation (penalising excessive stretching or shrinking).  The best solutions can be obtained by knowing the approximate amount of stretching that is needed (e.g. ICBM templates are slightly bigger than typical brains, so greater zooms are likely to be needed). If registering to an image in ICBM/MNI space, then choose the first option.  If registering to a template that is close in size, then select the second option.  If you do not want to regularise, then choose the third.',''};
regtype.labels = {
                  'ICBM space template'
                  'Average sized template'
                  'No regularisation'
}';
regtype.values = {
                  'mni'
                  'subj'
                  'none'
}';
regtype.def     = @(val)suit_get_defaults('normalise.estimate.regtype', val{:});
% ---------------------------------------------------------------------
% cutoff Nonlinear Frequency Cutoff
% ---------------------------------------------------------------------
cutoff         = cfg_entry;
cutoff.tag     = 'cutoff';
cutoff.name    = 'Nonlinear Frequency Cutoff';
cutoff.help    = {'Cutoff of DCT bases.  Only DCT bases of periods longer than the cutoff are used to describe the warps.',...
                'SUIT normalization uses higher frequencies (1cm) than normal MNI normalization',''};
cutoff.strtype = 'e';
cutoff.num     = [1 1];
cutoff.def     = @(val)suit_get_defaults('normalise.estimate.cutoff', val{:});
% ---------------------------------------------------------------------
% nits Nonlinear Iterations
% ---------------------------------------------------------------------
nits         = cfg_entry;
nits.tag     = 'nits';
nits.name    = 'Nonlinear Iterations';
nits.help    = {'Number of iterations of nonlinear warping performed. 30 is a good compromise between convergence and speed.',''};
nits.strtype = 'w';
nits.num     = [1 1];
nits.def     = @(val)suit_get_defaults('normalise.estimate.nits', val{:});
% ---------------------------------------------------------------------
% reg Nonlinear Regularisation
% ---------------------------------------------------------------------
reg         = cfg_entry;
reg.tag     = 'reg';
reg.name    = 'Nonlinear Regularisation';
reg.help    = {'The amount of regularisation for the nonlinear part of the spatial normalisation. Pick a value around one.',...
    'However, if your normalised images appear distorted, then it may be an idea to increase the amount of regularisation.',...
    '(by an order of magnitude). This may be especially useful if you have a lesioned cerebellum.',''};
reg.strtype = 'e';
reg.num     = [1 1];
reg.def     = @(val)suit_get_defaults('normalise.estimate.reg', val{:});
% ---------------------------------------------------------------------
% eoptions Estimation Options
% ---------------------------------------------------------------------
eoptions         = cfg_branch;
eoptions.tag     = 'estimate';
eoptions.name    = 'Estimation Options';
eoptions.val     = {smosrc smoref regtype cutoff nits reg };
eoptions.help    = {'Various settings for estimating warps.'};
eoptions.expanded = false; 

% ---------------------------------------------------------------------
% preserve Preserve
% ---------------------------------------------------------------------
preserveN         = cfg_menu;
preserveN.tag     = 'preserve';
preserveN.name    = 'Preserve';
preserveN.help    = {'Preserve Concentrations: Spatially normalised images are not "modulated". The warped images preserve the intensities of the original images.',...
                    'Preserve Total: Spatially normalised images are "modulated" in order to preserve the total amount of signal in the images. Areas that are expanded during warping are correspondingly reduced in intensity. Use this option for VBM analysis.',''};
preserveN.labels = {
                   'Preserve Concentrations'
                   'Preserve Amount'
}';
preserveN.values = {0 1};
preserveN.def     = @(val)suit_get_defaults('normalise.write.preserve', val{:});
% ---------------------------------------------------------------------
% bb Bounding box
% ---------------------------------------------------------------------
bbN         = cfg_entry;
bbN.tag     = 'bb';
bbN.name    = 'Bounding box';
bbN.help    = {'The bounding box (in mm) of the volume which is to be written (relative to the anterior commissure).',...
    'The default bounding box gives you images of the size of the suit template',''};
bbN.strtype = 'e';
bbN.num     = [2 3];
bbN.def     = @(val)suit_get_defaults('normalise.write.bb', val{:});
% ---------------------------------------------------------------------
% vox Voxel sizes
% ---------------------------------------------------------------------
voxN         = cfg_entry;
voxN.tag     = 'vox';
voxN.name    = 'Voxel sizes';
voxN.help    = {'The voxel sizes (x, y & z, in mm) of the written normalised images.',...
    'Usually we use 2mm for functional data and 1mm for anatomical data.',''};
voxN.strtype = 'e';
voxN.num     = [1 3];
voxN.def     = @(val)suit_get_defaults('normalise.write.vox', val{:});
% ---------------------------------------------------------------------
% interp Interpolation
% ---------------------------------------------------------------------
interpN         = cfg_menu;
interpN.tag     = 'interp';
interpN.name    = 'Interpolation';
interpN.help    = {
                  'The method by which the images are sampled when being written in a different space.'
                  '    Nearest Neighbour:     - Fastest, but not normally recommended.'
                  '    Bilinear Interpolation:     - OK for PET, or realigned fMRI.'
                  '    B-spline Interpolation:     - Better quality (but slower) interpolation/* \cite{thevenaz00a}*/, especially       with higher degree splines.  Do not use B-splines when       there is any region of NaN or Inf in the images. '
}';
interpN.labels = {
                 'Nearest neighbour'
                 'Trilinear'
                 '2nd Degree B-spline'
                 '3rd Degree B-Spline '
                 '4th Degree B-Spline '
                 '5th Degree B-Spline'
                 '6th Degree B-Spline'
                 '7th Degree B-Spline'
}';
interpN.values = {0 1 2 3 4 5 6 7};
interpN.def     = @(val)suit_get_defaults('normalise.write.interp', val{:});
% ---------------------------------------------------------------------
% wrap Wrapping
% ---------------------------------------------------------------------
wrapN         = cfg_menu;
wrapN.tag     = 'wrap';
wrapN.name    = 'Wrapping';
wrapN.help    = {
                'These are typically:'
                '    No wrapping: for PET or images that have already                   been spatially transformed. '
                '    Wrap in  Y: for (un-resliced) MRI where phase encoding                   is in the Y direction (voxel space).'
}';
wrapN.labels = {
               'No wrap'
               'Wrap X'
               'Wrap Y'
               'Wrap X & Y'
               'Wrap Z'
               'Wrap X & Z'
               'Wrap Y & Z'
               'Wrap X, Y & Z'
}';
wrapN.values = {[0 0 0] [1 0 0] [0 1 0] [1 1 0] [0 0 1] [1 0 1] [0 1 1]...
               [1 1 1]};
wrapN.def     = @(val)suit_get_defaults('normalise.write.wrap', val{:});

% ---------------------------------------------------------------------
% roptions Writing Options
% ---------------------------------------------------------------------
roptionsN         = cfg_branch;
roptionsN.tag     = 'write';
roptionsN.name    = 'Writing Options';
roptionsN.val     = {preserveN bbN voxN interpN wrapN };
roptionsN.help    = {'Various options for writing the normalised anatomical image.'};
roptionsN.expanded = false; 
% ---------------------------------------------------------------------
% est Normalise: Estimate
% ---------------------------------------------------------------------
normalise         = cfg_exbranch;
normalise.tag     = 'normalise';
normalise.name    = 'Normalise into SUIT space';
normalise.val     = {subjsN prefixN template weight parampostfix smooth_mask eoptions roptionsN };
normalise.help    = {'Computes the warp that best registers a source image (or series of source images) to match a template, saving it to a file imagename''_sn.mat''.'};
normalise.prog = @suit_run_normalise;
normalise.vout = @vout_normalise;



% ---------------------------------------------------------------------
%
% NORMALIZATION - DARTEL (s) 
% 
% source Source Image
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Gray matter images 
% ---------------------------------------------------------------------
imageD1         = cfg_files;
imageD1.tag     = 'gray';
imageD1.name    = 'Gray matter image';
imageD1.help    = {'Select the gray matter images from the segmentation set (c1_*)'};
imageD1.filter = 'image';
imageD1.ufilter = '.*';
imageD1.num     = [1 1];
% ---------------------------------------------------------------------
% White matter images 
% ---------------------------------------------------------------------
imageD2         = cfg_files;
imageD2.tag     = 'white';
imageD2.name    = 'White matter image';
imageD2.help    = {'Select the white matter images from the segmentation set (c2_*)'};
imageD2.filter = 'image';
imageD2.ufilter = '.*';
imageD2.num     = [1 1];
% ---------------------------------------------------------------------
% White matter images 
% ---------------------------------------------------------------------
imageD3         = cfg_files;
imageD3.tag     = 'isolation';
imageD3.name    = 'Isolation mask';
imageD3.help    = {'Select the isolation mask for each subject (pcereb_corr)'};
imageD3.filter = 'image';
imageD3.ufilter = '.*';
imageD3.num     = [1 1];
% ---------------------------------------------------------------------
% subj Subject
% ---------------------------------------------------------------------
subjND         = cfg_branch;
subjND.tag     = 'subjND';
subjND.name    = 'Subject';
subjND.val     = {imageD1 imageD2 imageD3};
subjND.help    = {'Data for this subject. '};
% ---------------------------------------------------------------------
% esubjs Data
% ---------------------------------------------------------------------
subjsND         = cfg_repeat;
subjsND.tag     = 'subjsND';
subjsND.name    = 'Data';
subjsND.help    = {'List of subjects for the normalization step. Images of each subject should be warped differently.',''};
subjsND.values  = {subjND};
subjsND.num     = [1 inf];

% ---------------------------------------------------------------------
% warp1 Run Dartel (existing Templates)
% ---------------------------------------------------------------------
normalise_dartel         = cfg_exbranch;
normalise_dartel.tag     = 'normalise_dartel';
normalise_dartel.name    = 'Normalise into SUIT space (using Dartel)';
normalise_dartel.val     = {subjsND};
normalise_dartel.help    = {'Run normalisation to the Suit template using Dartel. This can lead to a more accurate normalisation than the older version of suit-normalise. '};
normalise_dartel.prog = @suit_normalize_dartel;


% ---------------------------------------------------------------------
%
% NORMALIZATION with dentate 
% 
% source Source Image
% ---------------------------------------------------------------------
dentateROI         = cfg_files;
dentateROI.tag     = 'dentateROI';
dentateROI.name    = 'dentate ROI image';
dentateROI.help    = {'An ROI drawing of the dentate nucleus aligned to the source image.'};
dentateROI.filter = 'image';
dentateROI.ufilter = '.*';
dentateROI.num     = [1 1];
% ---------------------------------------------------------------------
% subj Subject
% ---------------------------------------------------------------------
subjND         = cfg_branch;
subjND.tag     = 'subjND';
subjND.name    = 'Subject';
subjND.val     = {imageD1 imageD2 imageD3 dentateROI};
subjND.help    = {'Data for this subject.  The same parameters are used for all subjects.'};
% ---------------------------------------------------------------------
% esubjs Data
% ---------------------------------------------------------------------
subjsND         = cfg_repeat;
subjsND.tag     = 'subjsND';
subjsND.name    = 'Data';
subjsND.help    = {'List of subjects for the normalization step. Images of each subject should be warped differently.'};
subjsND.values  = {subjND};
subjsND.num     = [1 Inf];

% ---------------------------------------------------------------------
% warp1 Run Dartel (existing Templates)
% ---------------------------------------------------------------------
normalise_dentate         = cfg_exbranch;
normalise_dentate.tag     = 'normalise_dentate';
normalise_dentate.name    = 'Normalise into SUIT space with dentate ROI';
normalise_dentate.val     = {subjsND};
normalise_dentate.help    = {'Run normalisation to the Suit template with dentate nucelus ROI using Dartel.' 'Please see Diedrichsen et al. (2011) as reference.'};
normalise_dentate.prog = @suit_normalize_dentate;


% ---------------------------------------------------------------------
% 
% RESLICE
% 
% matname Parameter File
% ---------------------------------------------------------------------
paramfile         = cfg_files;
paramfile.tag     = 'paramfile';
paramfile.name    = 'Parameter File';
paramfile.help    = {'Select the ''_snc.mat'' file containing the spatial normalisation parameters for that subject.' ''};
paramfile.filter = 'mat';
paramfile.ufilter = '.mat$';
paramfile.num     = [1 1];
% ---------------------------------------------------------------------
% resample Images to Write
% ---------------------------------------------------------------------
resample         = cfg_files;
resample.tag     = 'resample';
resample.name    = 'Images to Write';
resample.help    = {'These are the images for warping according to the estimated parameters. They can be any images that are in register with the "source" image used to generate the parameters.'};
resample.filter = 'image';
resample.ufilter = '.*';
resample.num     = [1 Inf];
% ---------------------------------------------------------------------
% Mask 
% ---------------------------------------------------------------------
maskR         = cfg_files;
maskR.tag     = 'mask';
maskR.name    = 'Cerebellar Mask';
maskR.val     = {''};
maskR.help    = {'Optional mask image that will be applied to the original images before reslicing.' 
                 'Typically this is the _pcereb_corr.nii from the isolation to remove all non-cerebellar structures.'
                 ''};
maskR.filter = 'image';
maskR.ufilter = '.*';
maskR.num     = [0 1];
% ---------------------------------------------------------------------
% subj Subject
% ---------------------------------------------------------------------
subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {paramfile resample maskR};
subj.help    = {'Data for this subject.  The same parameters are used within subject.'};
% ---------------------------------------------------------------------
% wsubjs Data
% ---------------------------------------------------------------------
wsubjs         = cfg_repeat;
wsubjs.tag     = 'wsubjs';
wsubjs.name    = 'Data';
wsubjs.help    = {'List of subjects. Images of each subject should be warped differently.'};
wsubjs.values  = {subj};
wsubjs.num     = [1 Inf];
% ---------------------------------------------------------------------
% preserve Preserve
% ---------------------------------------------------------------------
preserve         = cfg_menu;
preserve.tag     = 'preserve';
preserve.name    = 'Preserve';
preserve.help    = {'Preserve Concentrations: Spatially normalised images are not "modulated". The warped images preserve the intensities of the original images.',...
                    'Preserve Total: Spatially normalised images are "modulated" in order to preserve the total amount of signal in the images. Areas that are expanded during warping are correspondingly reduced in intensity. Use this option for VBM analysis.',''};

preserve.labels = {
                   'Preserve Concentrations'
                   'Preserve Amount'
}';
preserve.values = {0 1};
preserve.def     = @(val)suit_get_defaults('reslice.preserve', val{:});
% ---------------------------------------------------------------------
% bb Bounding box
% ---------------------------------------------------------------------
bbR          = cfg_entry;
bbR.tag     = 'bb';
bbR.name    = 'Bounding box';
bbR.help    = {'The bounding box (in mm) of the volume which is to be written (relative to the anterior commissure).',...
    'The default bounding box gives you images of the size of the suit template',''};
bbR.strtype = 'e';
bbR.num     = [2 3];
bbR.def     = @(val)suit_get_defaults('reslice.bb', val{:});
% ---------------------------------------------------------------------
% vox Voxel sizes
% ---------------------------------------------------------------------
voxR         = cfg_entry;
voxR.tag     = 'vox';
voxR.name    = 'Voxel sizes';
voxR.help    = {'The voxel sizes (x, y & z, in mm) of the written normalised images.',...
    'Usually we use 2mm for functional data and 1mm for anatomical data.',''};
voxR.strtype = 'e';
voxR.num     = [1 3];
voxR.def     = @(val)suit_get_defaults('reslice.vox', val{:});
% ---------------------------------------------------------------------
% interp Interpolation
% ---------------------------------------------------------------------
interpR         = cfg_menu;
interpR.tag     = 'interp';
interpR.name    = 'Interpolation';
interpR.help    = {
                  'The method by which the images are sampled when being written in a different space.'
                  '    Nearest Neighbour:     - Recommended for reslicing Atlas label images or ROIs.'
                  '    Bilinear Interpolation:     - OK for PET, or realigned fMRI.'
                  '    B-spline Interpolation:     - Better quality (but slower) interpolation/* \cite{thevenaz00a}*/, especially       with higher degree splines.  Do not use B-splines when       there is any region of NaN or Inf in the images. '
}';
interpR.labels = {
                 'Nearest neighbour'
                 'Trilinear'
                 '2nd Degree B-spline'
                 '3rd Degree B-Spline '
                 '4th Degree B-Spline '
                 '5th Degree B-Spline'
                 '6th Degree B-Spline'
                 '7th Degree B-Spline'
}';
interpR.values = {0 1 2 3 4 5 6 7};
interpR.def     = @(val)suit_get_defaults('reslice.interp', val{:});
% ---------------------------------------------------------------------
% prefix Filename Prefix
% ---------------------------------------------------------------------
prefixR         = cfg_entry;
prefixR.tag     = 'prefix';
prefixR.name    = 'Filename Prefix';
prefixR.help    = {'Specify the string to be prepended to the filenames of the normalised image file(s). Default prefix is ''wc''.',''};
prefixR.strtype = 's';
prefixR.num     = [1 Inf];
prefixR.def     = @(val)suit_get_defaults('reslice.prefix', val{:});
% ---------------------------------------------------------------------
% smosrc Source Image Smoothing
% ---------------------------------------------------------------------
smooth_mask         = cfg_entry;
smooth_mask.tag     = 'smooth_mask';
smooth_mask.name    = 'Mask Image Smoothing';
smooth_mask.help    = {'Smoothing to apply to the mask image, before multiplying with the image to reslice',''};
smooth_mask.strtype = 'e';
smooth_mask.num     = [1 1];
smooth_mask.def     = @(val)suit_get_defaults('reslice.smooth_mask', val{:});
% ---------------------------------------------------------------------
% REslice
% ---------------------------------------------------------------------
reslice         = cfg_exbranch;
reslice.tag     = 'reslice';
reslice.name    = 'Reslice into SUIT space';
reslice.val     = {wsubjs smooth_mask preserve bbR voxR interpR prefixR};
reslice.help    = {'Allows previously estimated warps (stored in imagename''_snc.mat'' files) to be applied to series of images.'};
reslice.prog = @suit_run_reslice;

% ---------------------------------------------------------------------
% 
% RESLICE_DARTEL
% 
% 
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
% Flowfield
% ---------------------------------------------------------------------
flowfield         = cfg_files;
flowfield.tag     = 'flowfield';
flowfield.name    = 'Flowfield file';
flowfield.help    = {'Select the flowfield file (u_*) containing the nonlinear transformation for that subject.' ''};
flowfield.filter = 'image';
flowfield.ufilter = '.*';
flowfield.num     = [1 1];

% ---------------------------------------------------------------------
% 
% ---------------------------------------------------------------------
affineTr         = cfg_files;
affineTr.tag     = 'affineTr';
affineTr.name    = 'Affine transformation matrix';
affineTr.help    = {'Select the affine transformation matrix, which was generated during suit_normalise_dartel' ...
                  'i.e. (Affine_<imagename>.mat)'};
affineTr.filter = 'mat';
affineTr.ufilter = '.mat$';
affineTr.num     = [1 1];


% ---------------------------------------------------------------------
% subj Subject
% ---------------------------------------------------------------------
subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {affineTr flowfield resample maskR};
subj.help    = {'Data for this subject.  The same parameters are used within subject.'};
% ---------------------------------------------------------------------
% wsubjs Data
% ---------------------------------------------------------------------
wsubjs         = cfg_repeat;
wsubjs.tag     = 'wsubjs';
wsubjs.name    = 'Data';
wsubjs.help    = {'List of subjects. Images of each subject should be warped differently.'};
wsubjs.values  = {subj};
wsubjs.num     = [1 Inf];

% ---------------------------------------------------------------------
% jactransf Modulation
% ---------------------------------------------------------------------
jactransf         = cfg_menu;
jactransf.tag     = 'jactransf';
jactransf.name    = 'Modulation';
jactransf.val     = {0};
jactransf.help    = {'This allows the spatially normalised images to be rescaled by the Jacobian determinants of the deformations. Note that the rescaling is only approximate for deformations generated using smaller numbers of time steps.'};
jactransf.labels  = {
                    'Pres. Concentration (No "modulation")'
                    'Pres. Amount ("Modulation")'
}';
jactransf.values  = {0 1};
% ---------------------------------------------------------------------
% K Time Steps
% ---------------------------------------------------------------------
K         = cfg_menu;
K.tag     = 'K';
K.name    = 'Time Steps';
K.val     = {6};
K.help    = {'The number of time points used for solving the partial differential equations.  Note that Jacobian determinants are not very accurate for very small numbers of time steps (less than about 16).'};
K.labels  = {
            '1'
            '2'
            '4'
            '8'
            '16'
            '32'
            '64'
            '128'
            '256'
            '512'
}';
K.values  = {0 1 2 3 4 5 6 7 8 9};
% ---------------------------------------------------------------------
% interp Interpolation
% ---------------------------------------------------------------------
interp         = cfg_menu;
interp.tag     = 'interp';
interp.name    = 'Interpolation';
interp.val     = {1};
interp.help    = {
                  ['The method by which the images are sampled when ',...
                  'being written in a different space. ',...
                  '(Note that Inf or NaN values are treated as zero, ',...
                  'rather than as missing data)']
                  '    Nearest Neighbour:'
                  '      - Fastest, but not normally recommended.'
                  '    Trilinear Interpolation:'
                  '      - OK for PET, realigned fMRI, or segmentations'
                  '    B-spline Interpolation:'
                  ['      - Better quality (but slower) interpolation',...
                  '/* \cite{thevenaz00a}*/, especially with higher ',...
                  'degree splines. Can produce values outside the ',...
                  'original range (e.g. small negative values from an ',...
                  'originally all positive image).']
}';
interp.labels  = {
                 'Nearest neighbour'
                 'Trilinear'
                 '2nd Degree B-spline'
                 '3rd Degree B-Spline '
                 '4th Degree B-Spline '
                 '5th Degree B-Spline'
                 '6th Degree B-Spline'
                 '7th Degree B-Spline'
}';
interp.values  = {0 1 2 3 4 5 6 7};
% ---------------------------------------------------------------------
% crt_warped Create Warped
% ---------------------------------------------------------------------
reslice_dartel         = cfg_exbranch;
reslice_dartel.tag     = 'reslice_dartel';
reslice_dartel.name    = 'Reslice into SUIT space using Dartel flowfield';
reslice_dartel.val     = {wsubjs jactransf K bbR voxR interp prefixR};
reslice_dartel.help    = {'This reslices anatomical or functional images into suit space, using a dartel flow field.'};
reslice_dartel.prog = @suit_reslice_dartel;


% ---------------------------------------------------------------------
% 
% RESLICE_INV
% 
% 
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% resample Images to Write
% ---------------------------------------------------------------------
resample         = cfg_files;
resample.tag     = 'resample';
resample.name    = 'Images to Write';
resample.help    = {'These are the images for warping according to the estimated parameters. They can be any images that are in register with the SUIT image,'
    'such as atlas image or region of interest.' };
resample.filter = 'image';
resample.ufilter = '.*';
resample.num     = [1 Inf];
% ---------------------------------------------------------------------
% resample Images to Write
% ---------------------------------------------------------------------
reference         = cfg_files;
reference.tag     = 'reference';
reference.name    = 'Image for reference';
reference.val     = {''};
reference.help    = {'Here you can provide an optimal image as reference for orientation and voxel size. If not given, the information form the parameter file (_snc.mat) is used.' };
reference.filter = 'image';
reference.ufilter = '.*';
reference.num     = [0 1];
% ---------------------------------------------------------------------
% prefix Filename Prefix
% ---------------------------------------------------------------------
prefixR         = cfg_entry;
prefixR.tag     = 'prefix';
prefixR.name    = 'Filename Prefix';
prefixR.help    = {'Specify the string to be prepended to the filenames of the normalised image file(s). Default prefix is ''i''.'};
prefixR.strtype = 's';
prefixR.num     = [1 Inf];
prefixR.def     = @(val)suit_get_defaults('reslice_inv.prefix', val{:});

interpR.def     = @(val)suit_get_defaults('reslice_inv.interp', val{:});

% ---------------------------------------------------------------------
% Reslice_inv
% ---------------------------------------------------------------------
reslice_inv         = cfg_exbranch;
reslice_inv.tag     = 'reslice_inv';
reslice_inv.name    = 'Reslice into individual space';
reslice_inv.val     = {resample paramfile prefixR interpR reference};
reslice_inv.help    = {'Allows previously estimated warps (stored in imagename''_snc.mat'' files) to be inverted and applied to images in Atlas space, such that they are resliced into individual space.',''};
reslice_inv.prog = @suit_run_reslice_inv;



% ---------------------------------------------------------------------
% 
% Summarize
% 
% 
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% resample Images to summarize 
% ---------------------------------------------------------------------
sample         = cfg_files;
sample.tag     = 'images';
sample.name    = 'Images to Summarize';
sample.help    = {'Select the images you want to summarize with into a table.' };
sample.filter = 'image';
sample.ufilter = '.*';
sample.num     = [1 Inf];
% ---------------------------------------------------------------------
% Atlas Image
% ---------------------------------------------------------------------
atlas         = cfg_files;
atlas.tag     = 'atlas';
atlas.name    = 'Atlas image';
atlas.help    = {'This image defines the ROIs over which you want to summarize your data.' };
atlas.filter = 'image';
atlas.ufilter = '.*';
atlas.num     = [1 1];
atlas.def     = @(val)suit_get_defaults('summarize.atlas', val{:});
% ---------------------------------------------------------------------
% output Output Filename
% ---------------------------------------------------------------------
output         = cfg_entry;
output.tag     = 'output';
output.name    = 'Output Filename';
output.help    = {'The table is written into a text file.'};
output.strtype = 's';
output.num     = [1 inf];
output.val     = {'table.txt'};
% ---------------------------------------------------------------------
% outdir Output Directory
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output Directory';
outdir.val{1} = {''};
outdir.help    = {'Output directory. If no directory is given, images will be written to current working directory. '};
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [0 1];
% ---------------------------------------------------------------------
% expression Expression
% ---------------------------------------------------------------------
stats         = cfg_entry;
stats.tag     = 'stats';
stats.name    = 'Statistic to be calculated (function name)';
stats.help    = {   ''
                    '' 
}';
stats.val = {'nanmean'}; 
stats.strtype = 's';
stats.num     = [2 Inf];
% REslice_inv
% ---------------------------------------------------------------------
summarize         = cfg_exbranch;
summarize.tag     = 'ROI';
summarize.name    = 'Summarize images by ROIs (atlas)';
summarize.val     = {sample atlas output outdir stats};
summarize.help    = {'Calculate statistics on images for different ROIs defined by an atas',''};
summarize.prog = @suit_run_ROI_summarize;


% ---------------------------------------------------------------------
% 
% Flatmap
% 
% 
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Volume to map 
% ---------------------------------------------------------------------
flatVol         = cfg_files;
flatVol.tag     = 'images';
flatVol.name    = 'Volume to plot';
flatVol.help    = {'Select the image that you want to plot.' };
flatVol.filter = 'image';
flatVol.ufilter = '.*';
flatVol.num     = [1 1];
% ---------------------------------------------------------------------
% Flatmap atlas space 
% ---------------------------------------------------------------------
flatSpace         = cfg_menu;
flatSpace.tag     = 'flatSpace';
flatSpace.name    = 'Atlas space of volume data';
flatSpace.val     = {1};
flatSpace.help    = {
                  ['Determines the atlas space that was used for the' ...
                  'normalisation of the volume data'],...
                  'SUIT is normal or dartel suit normalisation',...
                  'FSL refers to FNIRT normalisation to MNI152nonlin',...
                  'SPM refers to combined Segmentation and normalisation'}';
flatSpace.labels  = {
                 'SUIT'
                 'FSL'
                 'SPM'
}';
flatSpace.values  = {1 2 3};

% ---------------------------------------------------------------------
% statistics for mapping 
% ---------------------------------------------------------------------
flatStats         = cfg_entry;
flatStats.tag     = 'stats';
flatStats.name    = 'Statistic to be calculated across the different depth';
flatStats.help    = {   'Use nanmean for general data',...
                    'minmax or max for a Glass-brain projection of t-values (too see all sign. clusters)',...
                    'mode for discrete labels'}'; 
flatStats.val = {'nanmean'}; 
flatStats.strtype = 's';
flatStats.num     = [2 Inf];

% ---------------------------------------------------------------------
% Color scale: minimal and maximal value 
% ---------------------------------------------------------------------
flatCscale         = cfg_entry;
flatCscale.tag     = 'flatCscale';
flatCscale.name    = 'Colorscale [min max]';
flatCscale.help    = {'Functional value that is assigned to the lowest color and highest color in the current colormap.',...
                       'NaNs force the scaling to thr minimum and maximim of the data'};
flatCscale.strtype = 'e';
flatCscale.num     = [1 2];
flatCscale.def     = @(x)nan(1,2);

% ---------------------------------------------------------------------
% Color scale: minimal and maximal value 
% ---------------------------------------------------------------------
flatThreshold         = cfg_entry;
flatThreshold.tag     = 'flatThreshold';
flatThreshold.name    = 'Functional threshold';
flatThreshold.help    = {'Functional threshold, values below this value are not displayed.',...
                        'NaN for unthresholded data'}; 
flatThreshold.strtype = 'e';
flatThreshold.num     = [1 1];
flatThreshold.def     = @(x)nan(1,1);

% ---------------------------------------------------------------------
% Color scale: minimal and maximal value 
% ---------------------------------------------------------------------
flatCmap         = cfg_entry;
flatCmap.tag     = 'Cmap';
flatCmap.name    = 'Colormap';
flatCmap.help    = {'Colormap to be used',...
                        'String of function name or variable'}; 
flatCmap.strtype = 's';
flatCmap.num     = [2 inf];
flatCmap.def     = @(x)'jet';

% ---------------------------------------------------------------------
% Plot Border? 
% ---------------------------------------------------------------------
flatBorder         = cfg_menu;
flatBorder.tag     = 'flatBorder';
flatBorder.name    = 'Plot Border?';
flatBorder.val     = {1};
flatBorder.help    = {'Determines whether the approximate lobular boundaries are plotted as dotted lines'}';
flatBorder.labels  = {'No','Yes'}';
flatBorder.values  = {0 1};

% Flatmap main branch 
% ---------------------------------------------------------------------
flatmap         = cfg_exbranch;
flatmap.tag     = 'flatmap';
flatmap.name    = 'Plot volume on cerebellar flatmap';
flatmap.val     = {flatVol flatSpace flatStats flatCscale flatThreshold flatCmap flatBorder};
flatmap.help    = {'Plots a volume of functional or label data on a cerebellar flatmap',''};
flatmap.prog = @suit_run_flatmap;


% ---------------------------------------------------------------------
% ISOLATE_NEO 
% source Source Image
% ---------------------------------------------------------------------
sourceNeo         = cfg_files;
sourceNeo.tag     = 'source';
sourceNeo.name    = 'Source Image';
sourceNeo.help    = {'Anatomical image (T2) to perform the isolation algorithm on.  The image will be cropped and a probabistic isolation will be performed. '};
sourceNeo.filter = 'image';
sourceNeo.ufilter = '.*';
sourceNeo.num     = [1 1];

% ---------------------------------------------------------------------
% subj Subject
% ---------------------------------------------------------------------
subjNeo         = cfg_branch;
subjNeo.tag     = 'subj';
subjNeo.name    = 'Subject';
subjNeo.val     = {sourceNeo};
subjNeo.help    = {'Data for this subject.  The same parameters are used within subject.'};
% % ---------------------------------------------------------------------
% % esubjs Data
% % ---------------------------------------------------------------------
esubjsNeo         = cfg_repeat;
esubjsNeo.tag     = 'esubjs';
esubjsNeo.name    = 'Data';
esubjsNeo.help    = {'List of source images to be isolated.'};
esubjsNeo.values  = {sourceNeo};
esubjsNeo.num     = [1 Inf];
% ---------------------------------------------------------------------
% maskp 
% ---------------------------------------------------------------------
mask_probNeo         = cfg_entry;
mask_probNeo.tag     = 'maskp';
mask_probNeo.name    = 'Mask probability';
mask_probNeo.help    = {'Higher values will result in a larger cerebellar mask.'};
mask_probNeo.strtype = 'e';
mask_probNeo.num     = [1 1];
mask_probNeo.def     = @(val)suit_get_defaults('isolateNeo.maskp', val{:});
% ---------------------------------------------------------------------
% keeptfiles 
% ---------------------------------------------------------------------
keepfNeo         = cfg_entry;
keepfNeo.tag     = 'keeptfiles';
keepfNeo.name    = 'Keep temporal files';
keepfNeo.help    = {'Keep intermidiate results including segmentations.'};
keepfNeo.strtype = 'e';
keepfNeo.num     = [1 1];
keepfNeo.def     = @(val)suit_get_defaults('isolateNeo.keeptfiles', val{:});
% ---------------------------------------------------------------------
% getV calculate volume
% ---------------------------------------------------------------------
getvNeo         = cfg_entry;
getvNeo.tag     = 'getV';
getvNeo.name    = 'Calculate volume';
getvNeo.help    = {'This option calculate the volume of the cerebelar mask. if the mask need manual correction volume should be calculated again.'};
getvNeo.strtype = 'e';
getvNeo.num     = [1 1];
getvNeo.def     = @(val)suit_get_defaults('isolateNeo.getV', val{:});

% ---------------------------------------------------------------------
% Isolate unit
% ---------------------------------------------------------------------
isolateNeo         = cfg_exbranch;
isolateNeo.tag     = 'isolateNeo';
isolateNeo.name    = 'Neonatal Isolation';
isolateNeo.val     = {esubjsNeo mask_probNeo keepfNeo getvNeo};
isolateNeo.help    = {'Isolation of the neonatal cerebellum from the surrounding tissue'}';
isolateNeo.prog = @suit_run_isolate_neo;

% ---------------------------------------------------------------------
%
% NORMALIZATION Neonate
% 
% source Source Image
% ---------------------------------------------------------------------
sourceNneo         = cfg_files;
sourceNneo.tag     = 'source';
sourceNneo.name    = 'Source Image';
sourceNneo.help    = {'The cropped image (c_*) that should be warped to match the template(s).',...
    'The result is a set of warps, which can be applied to this image, or any other image that is in register with it.',''};
sourceNneo.filter = 'image';
sourceNneo.ufilter = '.*';
sourceNneo.num     = [1 1];
% ---------------------------------------------------------------------
% Mask 
% ---------------------------------------------------------------------
maskNneo         = cfg_files;
maskNneo.tag     = 'mask';
maskNneo.name    = 'Cerebellar Mask';
maskNneo.help    = {'The isolation mask derived from the suit_isolate algorithm that removes all non-cerebellar structures.',...
    'The _pcereb_corr.nii should be viewed and possibly handcorrected before normalization.',''};
maskNneo.filter = 'image';
maskNneo.ufilter = '.*';
maskNneo.num     = [1 1];
% ---------------------------------------------------------------------
% lesion mask 
% ---------------------------------------------------------------------
lesionMaskNneo         = cfg_files;
lesionMaskNneo.tag     = 'lesion_mask';
lesionMaskNneo.name    = 'Lesion Mask';
lesionMaskNneo.val     = {''};
lesionMaskNneo.help    = {'Optional lesion mask for the cerebellum. If specified, the lesion will be extended by a rim of 3 mm from the normalization process.',''};
lesionMaskNneo.filter = 'image';
lesionMaskNneo.ufilter = '.*';
lesionMaskNneo.num     = [0 1];
% ---------------------------------------------------------------------
% subj Subject
% ---------------------------------------------------------------------
subjNneo         = cfg_branch;
subjNneo.tag     = 'subjN';
subjNneo.name    = 'Subject';
subjNneo.val     = {sourceNneo maskNneo lesionMaskNneo};
subjNneo.help    = {'Data for this subject.  The same parameters are used for all subjects.'};
% ---------------------------------------------------------------------
% esubjs Data
% ---------------------------------------------------------------------
subjsNneo         = cfg_repeat;
subjsNneo.tag     = 'subjsN';
subjsNneo.name    = 'Data';
subjsNneo.help    = {'List of subjects for the normalization step. Images of each subject should be warped differently.',''};
subjsNneo.values  = {subjNneo};
subjsNneo.num     = [1 Inf];
% ---------------------------------------------------------------------
% template Template Image
% ---------------------------------------------------------------------
templateneo         = cfg_files;
templateneo.tag     = 'template';
templateneo.name    = 'Template Image';
templateneo.help    = {'Specify a template image to match the source image with (default is SUIT.nii).'};
templateneo.filter = 'image';
templateneo.ufilter = '.*';
templateneo.num     = [1 1]; 
templateneo.def     = @(val)suit_get_defaults('normaliseN.template', val{:});
% ---------------------------------------------------------------------
% weight Template Weighting Image
% ---------------------------------------------------------------------
weightneo         = cfg_files;
weightneo.tag     = 'template_weight';
weightneo.name    = 'Template Weighting Image';
weightneo.help    = {
                  'Applies a weighting mask to the template(s) during the parameter estimation.',''};
weightneo.filter = 'image';
weightneo.ufilter = '.*';
weightneo.num     = [1 1];
weightneo.def     = @(val)suit_get_defaults('normaliseN.template_weight', val{:});
% ---------------------------------------------------------------------
% smosrc Source Image Smoothing
% ---------------------------------------------------------------------
prefixNneo         = cfg_entry;
prefixNneo.tag     = 'prefix';
prefixNneo.name    = 'Prefix anatomical';
prefixNneo.help    = {'Prefix with which the resliced anantomical image will be saved.',''};
prefixNneo.strtype = 's';
prefixNneo.num     = [1 Inf];
prefixNneo.def     = @(val)suit_get_defaults('normaliseN.prefix', val{:});
% ---------------------------------------------------------------------
% smosrc Source Image Smoothing
% ---------------------------------------------------------------------
parampostfixneo         = cfg_entry;
parampostfixneo.tag     = 'param_postfix';
parampostfixneo.name    = 'parameter postfix';
parampostfixneo.help    = {'Post fix for the deformation parameter file (default: _snc)',''};
parampostfixneo.strtype = 's';
parampostfixneo.num     = [1 Inf];
parampostfixneo.def     = @(val)suit_get_defaults('normaliseN.param_postfix', val{:});

% ---------------------------------------------------------------------
% smosrc Source Image Smoothing
% ---------------------------------------------------------------------
smooth_maskneo         = cfg_entry;
smooth_maskneo.tag     = 'smooth_mask';
smooth_maskneo.name    = 'Mask Image Smoothing';
smooth_maskneo.help    = {'Smoothing to apply to the mask image, before multiplying with the target image.',''};
smooth_maskneo.strtype = 'e';
smooth_maskneo.num     = [1 1];
smooth_maskneo.def     = @(val)suit_get_defaults('normaliseN.smooth_mask', val{:});
% ---------------------------------------------------------------------
% smosrc Source Image Smoothing
% ---------------------------------------------------------------------
smosrcneo         = cfg_entry;
smosrcneo.tag     = 'smosrc';
smosrcneo.name    = 'Source Image Smoothing';
smosrcneo.help    = {'Smoothing to apply to a copy of the source image.',...
    'The template and source images should have approximately the same smoothness.',...
     'By experience a value of 2mm gives good result while preserving anatomical detail.',''}; 
smosrcneo.strtype = 'e';
smosrcneo.num     = [1 1];
smosrcneo.def     = @(val)suit_get_defaults('normaliseN.estimate.smosrc', val{:});
% ---------------------------------------------------------------------
% smoref Template Image Smoothing
% ---------------------------------------------------------------------
smorefneo         = cfg_entry;
smorefneo.tag     = 'smoref';
smorefneo.name    = 'Template Image Smoothing';
smorefneo.help    = {'Smoothing to apply to a copy of the template image.','',...
    'The template and source images should have approximately the same smoothness.',...
    'On the standard SUIT image, no smoothing of the template is necessary',''};
smorefneo.strtype = 'e';
smorefneo.num     = [1 1];
smorefneo.def     = @(val)suit_get_defaults('normaliseN.estimate.smoref', val{:});
% ---------------------------------------------------------------------
% regtype Affine Regularisation
% ---------------------------------------------------------------------
regtypeneo         = cfg_menu;
regtypeneo.tag     = 'regtype';
regtypeneo.name    = 'Affine Regularisation';
regtypeneo.help    = {'Affine registration into a standard space can be made more robust by regularisation (penalising excessive stretching or shrinking).  The best solutions can be obtained by knowing the approximate amount of stretching that is needed (e.g. ICBM templates are slightly bigger than typical brains, so greater zooms are likely to be needed). If registering to an image in ICBM/MNI space, then choose the first option.  If registering to a template that is close in size, then select the second option.  If you do not want to regularise, then choose the third.',''};
regtypeneo.labels = {
                  'ICBM space template'
                  'Average sized template'
                  'No regularisation'
}';
regtypeneo.values = {
                  'mni'
                  'subj'
                  'none'
}';
regtypeneo.def     = @(val)suit_get_defaults('normaliseN.estimate.regtype', val{:});
% ---------------------------------------------------------------------
% cutoff Nonlinear Frequency Cutoff
% ---------------------------------------------------------------------
cutoffneo         = cfg_entry;
cutoffneo.tag     = 'cutoff';
cutoffneo.name    = 'Nonlinear Frequency Cutoff';
cutoffneo.help    = {'Cutoff of DCT bases.  Only DCT bases of periods longer than the cutoff are used to describe the warps.',...
                'SUIT normalization uses higher frequencies (1cm) than normal MNI normalization',''};
cutoffneo.strtype = 'e';
cutoffneo.num     = [1 1];
cutoffneo.def     = @(val)suit_get_defaults('normaliseN.estimate.cutoff', val{:});
% ---------------------------------------------------------------------
% nits Nonlinear Iterations
% ---------------------------------------------------------------------
nitsneo         = cfg_entry;
nitsneo.tag     = 'nits';
nitsneo.name    = 'Nonlinear Iterations';
nitsneo.help    = {'Number of iterations of nonlinear warping performed. 30 is a good compromise between convergence and speed.',''};
nitsneo.strtype = 'w';
nitsneo.num     = [1 1];
nitsneo.def     = @(val)suit_get_defaults('normaliseN.estimate.nits', val{:});
% ---------------------------------------------------------------------
% reg Nonlinear Regularisation
% ---------------------------------------------------------------------
regneo         = cfg_entry;
regneo.tag     = 'reg';
regneo.name    = 'Nonlinear Regularisation';
regneo.help    = {'The amount of regularisation for the nonlinear part of the spatial normalisation. Pick a value around one.',...
    'However, if your normalised images appear distorted, then it may be an idea to increase the amount of regularisation.',...
    '(by an order of magnitude). This may be especially useful if you have a lesioned cerebellum.',''};
regneo.strtype = 'e';
regneo.num     = [1 1];
regneo.def     = @(val)suit_get_defaults('normaliseN.estimate.reg', val{:});
% ---------------------------------------------------------------------
% eoptions Estimation Options
% ---------------------------------------------------------------------
eoptionsneo         = cfg_branch;
eoptionsneo.tag     = 'estimate';
eoptionsneo.name    = 'Estimation Options';
eoptionsneo.val     = {smosrcneo smorefneo regtypeneo cutoffneo nitsneo regneo };
eoptionsneo.help    = {'Various settings for estimating warps.'};
eoptionsneo.expanded = false; 

% ---------------------------------------------------------------------
% preserve Preserve
% ---------------------------------------------------------------------
preserveNneo         = cfg_menu;
preserveNneo.tag     = 'preserveN';
preserveNneo.name    = 'Preserve';
preserveNneo.help    = {'Preserve Concentrations: Spatially normalised images are not "modulated". The warped images preserve the intensities of the original images.',...
                    'Preserve Total: Spatially normalised images are "modulated" in order to preserve the total amount of signal in the images. Areas that are expanded during warping are correspondingly reduced in intensity. Use this option for VBM analysis.',''};
preserveNneo.labels = {
                   'Preserve Concentrations'
                   'Preserve Amount'
}';
preserveNneo.values = {0 1};
preserveNneo.def     = @(val)suit_get_defaults('normaliseN.write.preserve', val{:});
% ---------------------------------------------------------------------
% bb Bounding box
% ---------------------------------------------------------------------
bbNneo         = cfg_entry;
bbNneo.tag     = 'bb';
bbNneo.name    = 'Bounding box';
bbNneo.help    = {'The bounding box (in mm) of the volume which is to be written (relative to the anterior commissure).',...
    'The default bounding box gives you images of the size of the suit template',''};
bbNneo.strtype = 'e';
bbNneo.num     = [2 3];
bbNneo.def     = @(val)suit_get_defaults('normaliseN.write.bb', val{:});
% ---------------------------------------------------------------------
% vox Voxel sizes
% ---------------------------------------------------------------------
voxNneo         = cfg_entry;
voxNneo.tag     = 'voxN';
voxNneo.name    = 'Voxel sizes';
voxNneo.help    = {'The voxel sizes (x, y & z, in mm) of the written normalised images.',...
    'Usually we use 2mm for functional data and 1mm for anatomical data.',''};
voxNneo.strtype = 'e';
voxNneo.num     = [1 3];
voxNneo.def     = @(val)suit_get_defaults('normaliseN.write.vox', val{:});
% ---------------------------------------------------------------------
% interp Interpolation
% ---------------------------------------------------------------------
interpNneo         = cfg_menu;
interpNneo.tag     = 'interpN';
interpNneo.name    = 'Interpolation';
interpNneo.help    = {
                  'The method by which the images are sampled when being written in a different space.'
                  '    Nearest Neighbour:     - Fastest, but not normally recommended.'
                  '    Bilinear Interpolation:     - OK for PET, or realigned fMRI.'
                  '    B-spline Interpolation:     - Better quality (but slower) interpolation/* \cite{thevenaz00a}*/, especially       with higher degree splines.  Do not use B-splines when       there is any region of NaN or Inf in the images. '
}';
interpNneo.labels = {
                 'Nearest neighbour'
                 'Trilinear'
                 '2nd Degree B-spline'
                 '3rd Degree B-Spline '
                 '4th Degree B-Spline '
                 '5th Degree B-Spline'
                 '6th Degree B-Spline'
                 '7th Degree B-Spline'
}';
interpNneo.values = {0 1 2 3 4 5 6 7};
interpNneo.def     = @(val)suit_get_defaults('normaliseN.write.interp', val{:});
% ---------------------------------------------------------------------
% wrap Wrapping
% ---------------------------------------------------------------------
wrapNneo         = cfg_menu;
wrapNneo.tag     = 'wrapN';
wrapNneo.name    = 'Wrapping';
wrapNneo.help    = {
                'These are typically:'
                '    No wrapping: for PET or images that have already                   been spatially transformed. '
                '    Wrap in  Y: for (un-resliced) MRI where phase encoding                   is in the Y direction (voxel space).'
}';
wrapNneo.labels = {
               'No wrap'
               'Wrap X'
               'Wrap Y'
               'Wrap X & Y'
               'Wrap Z'
               'Wrap X & Z'
               'Wrap Y & Z'
               'Wrap X, Y & Z'
}';
wrapNneo.values = {[0 0 0] [1 0 0] [0 1 0] [1 1 0] [0 0 1] [1 0 1] [0 1 1]...
               [1 1 1]};
wrapNneo.def     = @(val)suit_get_defaults('normaliseN.write.wrap', val{:});

% ---------------------------------------------------------------------
% roptions Writing Options
% ---------------------------------------------------------------------
roptionsNneo         = cfg_branch;
roptionsNneo.tag     = 'write';
roptionsNneo.name    = 'Writing Options';
roptionsNneo.val     = {preserveNneo bbNneo voxNneo interpNneo wrapNneo };
roptionsNneo.help    = {'Various options for writing the normalised anatomical image.'};
roptionsNneo.expanded = false; 
% ---------------------------------------------------------------------
% est Normalise: Estimate
% ---------------------------------------------------------------------
normaliseNeo         = cfg_exbranch;
normaliseNeo.tag     = 'normaliseNeo';
normaliseNeo.name    = 'Normalise into SUIT-N space';
normaliseNeo.val     = {subjsNneo prefixNneo templateneo weightneo parampostfixneo smooth_maskneo eoptionsneo roptionsNneo };
normaliseNeo.help    = {'Computes the warp that best registers a source image (or series of source images) to match a template, saving it to a file imagename''_sn.mat''.'};
normaliseNeo.prog = @suit_run_normalise;
normaliseNeo.vout = @vout_normalise;






% ---------------------------------------------------------------------
% 
% RESLICE NEO
% 
% matname Parameter File
% ---------------------------------------------------------------------
paramfileneo         = cfg_files;
paramfileneo.tag     = 'paramfile';
paramfileneo.name    = 'Parameter File';
paramfileneo.help    = {'Select the ''_snc.mat'' file containing the spatial normalisation parameters for that subject.' ''};
paramfileneo.filter = 'mat';
paramfileneo.ufilter = '.mat$';
paramfileneo.num     = [1 1];
% ---------------------------------------------------------------------
% resample Images to Write
% ---------------------------------------------------------------------
resampleneo         = cfg_files;
resampleneo.tag     = 'resample';
resampleneo.name    = 'Images to Write';
resampleneo.help    = {'These are the images for warping according to the estimated parameters. They can be any images that are in register with the "source" image used to generate the parameters.'};
resampleneo.filter = 'image';
resampleneo.ufilter = '.*';
resampleneo.num     = [1 Inf];
% ---------------------------------------------------------------------
% Mask 
% ---------------------------------------------------------------------
maskRneo         = cfg_files;
maskRneo.tag     = 'mask';
maskRneo.name    = 'Cerebellar Mask';
maskRneo.val     = {''};
maskRneo.help    = {'Optional mask image that will be applied to the original images before reslicing.' 
                 'Typically this is the _pcereb_corr.nii from the isolation to remove all non-cerebellar structures.'
                 ''};
maskRneo.filter = 'image';
maskRneo.ufilter = '.*';
maskRneo.num     = [0 1];
% ---------------------------------------------------------------------
% subj Subject
% ---------------------------------------------------------------------
subjneo         = cfg_branch;
subjneo.tag     = 'subj';
subjneo.name    = 'Subject';
subjneo.val     = {paramfileneo resampleneo maskRneo};
subjneo.help    = {'Data for this subject.  The same parameters are used within subject.'};
% ---------------------------------------------------------------------
% wsubjs Data
% ---------------------------------------------------------------------
wsubjsneo         = cfg_repeat;
wsubjsneo.tag     = 'wsubjs';
wsubjsneo.name    = 'Data';
wsubjsneo.help    = {'List of subjects. Images of each subject should be warped differently.'};
wsubjsneo.values  = {subjneo};
wsubjsneo.num     = [1 Inf];
% ---------------------------------------------------------------------
% preserve Preserve
% ---------------------------------------------------------------------
preserveneo         = cfg_menu;
preserveneo.tag     = 'preserve';
preserveneo.name    = 'Preserve';
preserveneo.help    = {'Preserve Concentrations: Spatially normalised images are not "modulated". The warped images preserve the intensities of the original images.',...
                    'Preserve Total: Spatially normalised images are "modulated" in order to preserve the total amount of signal in the images. Areas that are expanded during warping are correspondingly reduced in intensity. Use this option for VBM analysis.',''};

preserveneo.labels = {
                   'Preserve Concentrations'
                   'Preserve Amount'
}';
preserveneo.values = {0 1};
preserveneo.def     = @(val)suit_get_defaults('resliceN.preserve', val{:});
% ---------------------------------------------------------------------
% bb Bounding box
% ---------------------------------------------------------------------
bbRneo          = cfg_entry;
bbRneo.tag     = 'bb';
bbRneo.name    = 'Bounding box';
bbRneo.help    = {'The bounding box (in mm) of the volume which is to be written (relative to the anterior commissure).',...
    'The default bounding box gives you images of the size of the suit template',''};
bbRneo.strtype = 'e';
bbRneo.num     = [2 3];
bbRneo.def     = @(val)suit_get_defaults('resliceN.bb', val{:});
% ---------------------------------------------------------------------
% vox Voxel sizes
% ---------------------------------------------------------------------
voxRneo         = cfg_entry;
voxRneo.tag     = 'vox';
voxRneo.name    = 'Voxel sizes';
voxRneo.help    = {'The voxel sizes (x, y & z, in mm) of the written normalised images.',...
    'Usually we use 2mm for functional data and 1mm for anatomical data.',''};
voxRneo.strtype = 'e';
voxRneo.num     = [1 3];
voxRneo.def     = @(val)suit_get_defaults('resliceN.vox', val{:});
% ---------------------------------------------------------------------
% interp Interpolation
% ---------------------------------------------------------------------
interpRneo         = cfg_menu;
interpRneo.tag     = 'interp';
interpRneo.name    = 'Interpolation';
interpRneo.help    = {
                  'The method by which the images are sampled when being written in a different space.'
                  '    Nearest Neighbour:     - Recommended for reslicing Atlas label images or ROIs.'
                  '    Bilinear Interpolation:     - OK for PET, or realigned fMRI.'
                  '    B-spline Interpolation:     - Better quality (but slower) interpolation/* \cite{thevenaz00a}*/, especially       with higher degree splines.  Do not use B-splines when       there is any region of NaN or Inf in the images. '
}';
interpRneo.labels = {
                 'Nearest neighbour'
                 'Trilinear'
                 '2nd Degree B-spline'
                 '3rd Degree B-Spline '
                 '4th Degree B-Spline '
                 '5th Degree B-Spline'
                 '6th Degree B-Spline'
                 '7th Degree B-Spline'
}';
interpRneo.values = {0 1 2 3 4 5 6 7};
interpRneo.def     = @(val)suit_get_defaults('resliceN.interp', val{:});
% ---------------------------------------------------------------------
% prefix Filename Prefix
% ---------------------------------------------------------------------
prefixRneo         = cfg_entry;
prefixRneo.tag     = 'prefix';
prefixRneo.name    = 'Filename Prefix';
prefixRneo.help    = {'Specify the string to be prepended to the filenames of the normalised image file(s). Default prefix is ''wc''.',''};
prefixRneo.strtype = 's';
prefixRneo.num     = [1 Inf];
prefixRneo.def     = @(val)suit_get_defaults('resliceN.prefix', val{:});
% ---------------------------------------------------------------------
% smosrc Source Image Smoothing
% ---------------------------------------------------------------------
smooth_maskneo         = cfg_entry;
smooth_maskneo.tag     = 'smooth_mask';
smooth_maskneo.name    = 'Mask Image Smoothing';
smooth_maskneo.help    = {'Smoothing to apply to the mask image, before multiplying with the image to reslice',''};
smooth_maskneo.strtype = 'e';
smooth_maskneo.num     = [1 1];
smooth_maskneo.def     = @(val)suit_get_defaults('resliceN.smooth_mask', val{:});
% ---------------------------------------------------------------------
% REslice neo
% ---------------------------------------------------------------------
resliceNeo         = cfg_exbranch;
resliceNeo.tag     = 'resliceNeo';
resliceNeo.name    = 'Reslice into SUIT-N space';
resliceNeo.val     = {wsubjsneo smooth_maskneo preserveneo bbRneo voxRneo interpRneo prefixRneo};
resliceNeo.help    = {'Allows previously estimated warps (stored in imagename''_snc.mat'' files) to be applied to series of images.'};
resliceNeo.prog = @suit_run_reslice;






% ---------------------------------------------------------------------
% 
% Overall menu and toolbox 
% 
% 
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% SUIT toolbox
% ---------------------------------------------------------------------
suit         = cfg_choice;
suit.tag     = 'suit';
suit.name    = 'SUIT';
suit.help    = {
    'The spatially unbiased infra-tentorial template (SUIT) toolbox provides a more accurate normalization of the cerebellum and the brainstem than other common normalization methods.'
    ''
    'For help and more information see: http://diedrichsenlab.org/imaging/suit.htm'
    ''
    'Please cite:'
    'Diedrichsen, J. (2006). A spatially unbiased atlas template of the human cerebellum. Neuroimage, 33, 1, p. 127-138.'
    'Diedrichsen, J., Balsters, J. H., Flavell, J., Cussans, E., & Ramnani, N. (2009). A probabilistic atlas of the human cerebellum. Neuroimage.'
    'Diedrichsen, J., Maderwald, S., Kuper, M., Thurling, M., Rabe, K., Gizewski, E. R., et al. (2011). ' 
    'Imaging the deep cerebellar nuclei: A probabilistic atlas and normalization procedure. Neuroimage.'
    'Diedrichsen, J. & Zotow, E. (2015). Surface-based display of volume-averaged cerebellar imaging data.'
    ''
    }';
suit.values  = {isolate_seg normalise_dartel normalise_dentate reslice_dartel summarize flatmap isolate normalise reslice reslice_inv isolateNeo normaliseNeo resliceNeo};


%------------------------------------------------------------------------
function dep = vout_normalise(job)
for k=1:numel(job.subjN)
    dep(k)            = cfg_dep;
    dep(k).sname      = sprintf('Norm Params File (Subj %d)',k);
    dep(k).src_output = substruct('()',{k},'.','params');
    dep(k).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end;

