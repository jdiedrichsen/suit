function suit_defaults
% Sets the defaults which are used by the SUIT toolbox
%
% FORMAT suit_defaults
%_______________________________________________________________________
%
% This file is intended to be customised for the site.
% Individual users can make copies which can be stored in their own
% matlab subdirectories. If ~/matlab is ahead of the SPM directory
% in the MATLABPATH, then the users own personal defaults are used.
%
% Care must be taken when modifying this file
% v.2.5: Compatibility with SPM job manager (26/08/2010) 
% v.2.5.1: SWD problem fixed from new releases of SPM 8 
%_______________________________________________________________________
% Copyright (C) 2010 
% Joern Diedrichsen (joern.diedrichsen@googlemail.com)

spm_Dir= fileparts(which('spm'));
global defaults_suit; 

% Priors
defaults_suit.isolate.template_dir=[spm_Dir '/toolbox/suit/templates'];
defaults_suit.isolate.prior_dir=[spm_Dir '/toolbox/suit/priors'];
defaults_suit.isolate.template1='T1.nii';
defaults_suit.isolate.template2='SUIT.nii';
defaults_suit.isolate.template_weight2='SUIT_weight.nii';
defaults_suit.isolate.prior1='MNIa_pcerebellum.nii';
defaults_suit.isolate.prior2='SUIT_pcerebellum.nii';
defaults_suit.isolate.gray_prior='gray_cereb.nii';
defaults_suit.isolate.white_prior='white_cereb.nii';
defaults_suit.isolate.csf_prior='csf_cereb.nii';
defaults_suit.isolate.affreg.smosrc = 8;
defaults_suit.isolate.affreg.regtype = 'mni';
defaults_suit.isolate.affreg.weight = '';
defaults_suit.isolate.affreg.WF='';
defaults_suit.isolate.affreg.smooth2=4;
defaults_suit.isolate.bb=[-76 76;-108 -6;-75 11]; % Bounding box for the cropped image in atlas space
defaults_suit.isolate.white_cutoff=0.2;
defaults_suit.isolate.prior_FWHM=[5 5 5];
defaults_suit.isolate.cerebral_range=3.5;       % distance from cortical white matter 
defaults_suit.isolate.cerebellar_range=2.5;     % distance from cerebellar white matter 
defaults_suit.isolate.range_FWHM=[3 3 3];       % smoothing on range map to approximate probability
defaults_suit.isolate.iterations=2;           % Number of iterations
defaults_suit.isolate.crop='crop';         % 'none': No cropping of image
defaults_suit.isolate.keeptempfiles=0;        % Keep temporary files (0/1)

% isolateNeo
defaults_suit.isolateNeo.getV=0;
defaults_suit.isolateNeo.keeptfiles=0;
defaults_suit.isolateNeo.maskp=0.5;
defaults_suit.isolateNeo.bb = [-39 40;-65 -6;-41 -2];
defaults_suit.isolateNeo.priors=[spm_Dir '/toolbox/suit/neo/neoTPM.nii'];

% isolate_seg
defaults_suit.isolate_seg.keeptempfiles=0;
defaults_suit.isolate_seg.maskp=0.2;
defaults_suit.isolate_seg.bb=[-76 76;-108 -6;-70 11];
defaults_suit.isolate_seg.prior_dir=[spm_Dir '/toolbox/suit/priors'];
defaults_suit.isolate_seg.priors='TPM_SUIT.nii';


defaults_suit.normalise.template_dir=[spm_Dir '/toolbox/suit/templates'];
defaults_suit.normalise.prior_dir=[spm_Dir '/toolbox/suit/priors'];
defaults_suit.normalise.template_dir=[spm_Dir '/toolbox/suit/templates'];
defaults_suit.normalise.prior_dir=[spm_Dir '/toolbox/suit/priors'];
defaults_suit.normalise.template={[defaults_suit.normalise.template_dir '/SUIT.nii']};
defaults_suit.normalise.template_weight={[defaults_suit.normalise.template_dir '/SUIT_weight.nii']};
defaults_suit.normalise.source_weight=[];
defaults_suit.normalise.lesion_mask=[];
defaults_suit.normalise.lesion_rim=2;
defaults_suit.normalise.prefix='wsuit_';
defaults_suit.normalise.param='';
defaults_suit.normalise.outfilename='';
defaults_suit.normalise.param_postfix='_snc';
defaults_suit.normalise.smooth_mask=2;

defaults_suit.normalise.estimate.smosrc=2;
defaults_suit.normalise.estimate.smoref=0;
defaults_suit.normalise.estimate.regtype='subj';
defaults_suit.normalise.estimate.cutoff=10;
defaults_suit.normalise.estimate.nits=30;
defaults_suit.normalise.estimate.reg=1;
defaults_suit.normalise.estimate.wtsrc   = 0;
defaults_suit.normalise.write.preserve=0;
defaults_suit.normalise.write.interp=1;
defaults_suit.normalise.write.vox=[1 1 1];
defaults_suit.normalise.write.bb=[[-70 -100 -75];[70 -6 11]];
defaults_suit.normalise.write.wrap = [0 0 0];

% Normalize dentate 
defaults_suit.normalise_dentate.template_dentate={[defaults_suit.normalise.template_dir '/SUIT.nii']};
defaults_suit.normalise_dentate.nucleus_weight=4;
defaults_suit.normalise_dentate.param_postfix='_snd';
defaults_suit.normalise_dentate.prefix='wsuitd_';

% Reslicing options 
defaults_suit.reslice.prefix='wc';
defaults_suit.reslice.interp = 1;
defaults_suit.reslice.preserve=0;
defaults_suit.reslice.vox=[2 2 2];
defaults_suit.reslice.bb=[[-70 -100 -75];[70 -6 11]];
defaults_suit.reslice.smooth_mask=2;

% Reslicing options for dartel 
defaults_suit.reslice_dartel.prefix='wd';
defaults_suit.reslice_dartel.K = 6;
defaults_suit.reslice_dartel.interp = 1;
defaults_suit.reslice_dartel.jactransf=0;
defaults_suit.reslice_dartel.vox=[1 1 1];
defaults_suit.reslice_dartel.bb=[[-70 -100 -75];[70 -6 11]];

% Inverse reslicing option 
defaults_suit.reslice_inv.prefix='i';
defaults_suit.reslice_inv.interp = 0;

% Lobuli Summarize function
defaults_suit.summarize.atlas_dir={[spm_Dir '/toolbox/cerebellar_atlases']};
defaults_suit.summarize.atlas={[defaults_suit.summarize.atlas_dir{1} '/Diedrichsen_2009/atl-Anatom_space-SUIT_dseg.nii']};
defaults_suit.summarize.stats={'nanmean'};

% normalize neonate
defaults_suit.normaliseN.template={[spm_Dir '/toolbox/suit/neo/SUIT-N.nii']};
defaults_suit.normaliseN.template_weight={[spm_Dir '/toolbox/suit/neo/SUIT-N_w.nii']};
defaults_suit.normaliseN.source_weight=[];
defaults_suit.normaliseN.lesion_mask=[];
defaults_suit.normaliseN.lesion_rim=2;
defaults_suit.normaliseN.prefix='wsn_';
defaults_suit.normaliseN.param='';
defaults_suit.normaliseN.outfilename='';
defaults_suit.normaliseN.param_postfix='_snc';
defaults_suit.normaliseN.smooth_mask=2;
defaults_suit.normaliseN.estimate.smosrc=2;
defaults_suit.normaliseN.estimate.smoref=0;
defaults_suit.normaliseN.estimate.regtype='subj';
defaults_suit.normaliseN.estimate.cutoff=10;
defaults_suit.normaliseN.estimate.nits=30;
defaults_suit.normaliseN.estimate.reg=1;
defaults_suit.normaliseN.estimate.wtsrc = 0;
defaults_suit.normaliseN.write.preserve=0;
defaults_suit.normaliseN.write.interp=1;
defaults_suit.normaliseN.write.vox=[1 1 1];
defaults_suit.normaliseN.write.bb=[-39 -65 -41; 40 -6 -2];
defaults_suit.normaliseN.write.wrap = [0 0 0];

% Reslicing neonate 
defaults_suit.resliceN.prefix='wn';
defaults_suit.resliceN.interp = 1;
defaults_suit.resliceN.preserve=1;
defaults_suit.resliceN.vox=[1 1 1];
defaults_suit.resliceN.bb=[-39 -65 -41;40 -6 -2];
defaults_suit.resliceN.smooth_mask=2;
