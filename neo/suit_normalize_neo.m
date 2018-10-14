function suit_normalize_neo(source,mask,varargin)
% Normalize a source image into SUIT-N template masking it with the
% resulting image from suit_isolate_neo.
%
% Input:    
%       source  T2 neonatal image
%       mask    <name>_pcereb.nii from suit_isolate_neo (corrected)
% Options:  
%       mod     Modulation, set to 1 to modulate or 0 to preserve(default)
%       vox     Voxel size of final image, default [1 1 1]
%       interp  Interpolation trilinear(default) set to 0 for NN
%       prefix  Prefix for the output
%
% Output:
%       wsn_<name>       Normalized image into SUIT-N space
%       m<name>_snc.mat  Parametric file of transformation
%
% Example
%
%       suit_normalize_neo('<T2name>.nii','c_<T2name>_pcereb.nii')
% ________________________________________________________________________
% Carlos R. Hernandez-Castillo 2018

% Set defaults
mod = 0;
vox = [1 1 1];
interp = 1;
prefix = 'wsn_';
spm_Dir= fileparts(which('spm'));
template_dir=[spm_Dir '/toolbox/suit/neo/'];

% Get the defaults from varargin
vararginoptions(varargin,{'mod','vox','interp','prefix'});

write.interp=interp;
write.vox=vox;
write.preserve=mod;
write.bb=[-39 -65 -41; 40 -6 -2];

suit_normalize(source,'mask',mask,...
    'template',[template_dir 'SUIT-N.nii'],...
    'template_weight',[template_dir 'SUIT-N_w.nii'],...
    'write',write,'prefix',prefix);

[sdir,sname]= spm_fileparts(source);
if isempty(sdir); sdir = pwd; end
[mdir,mname]= spm_fileparts(mask);
if isempty(mdir); mdir = pwd; end
delete([sdir '/m' sname '.nii']);
delete([mdir '/s' mname '.nii']);

end