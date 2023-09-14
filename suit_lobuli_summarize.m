function D=suit_lobuli_summarize(images,varargin);
% function varargout=suit_ROI_summarize(images,varargin);
% Uses cerebellar ROIs from an atlas image to summarize cerebellar data 
%__________________________________________________________________________
% INPUT:
%   images: character or cell array of images to summarize
%__________________________________________________________________________
% OPTIONS:
%   'atlas',atlas_image: Atlas image (default: MDTB_10regions.nii)
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

global defaults_suit;
spm_dir = fileparts(which('spm'));
atlas=defaults_suit.summarize.atlas{1};
stats={'nanmean'};
regionname={'Left I_IV','Right I_IV','Left V','Right V',...
    'Left VI','Vermis VI','Right VI',...
    'Left Crus I','Vermis Crus I','Right Crus I',...
    'Left Crus II','Vermis Crus II','Right Crus II',...
    'Left VIIb','Vermis VIIb','Right VIIb',...
    'Left VIIIa','Vermis VIIIa','Right VIIIa',...
    'Left VIIIb','Vermis VIIIb','Right VIIIb',...
    'Left IX','Vermis IX','Right IX',...
    'Left X','Vermis X','Right X'};

D=suit_ROI_summarize(images,'atlas',atlas,'regionname',regionname,varargin(:)); 
