function suit_reslice_neo(source,mask,param,varargin)
% Reslice a source image into SUIT-N template masking it with the
% resulting image from suit_isolate_neo and the parameter file from
% suit_normalize_neo.
%
% Input:    
%       source  Array of image(s) to reslice using the same transformation
%       mask    <name>_pcereb.nii from suit_isolate_neo (corrected)
%       param   m<name>_snc.mat from suit_normalize_neo
%
% Options:  
%       mod     Modulation, set to 1 to modulate(default) or 0 to preserve
%       vox     Voxel size of final image, default [1 1 1]
%       interp  Interpolation trilinear(default) set to 0 for NN
%       prefix  Prefix for the output
%
% Output:
%       wn_<name>   Resliced image(s) into SUIT-N space
%
% Example
%
%   suit_reslice_neo({'<name>.nii'},'c_<name>_pcereb.nii','<name>_snc.mat')
% ________________________________________________________________________
% Carlos R. Hernandez-Castillo 2018

% Set defaults
vox = [1 1 1];
prefix = 'wn';
mod = 1;
interp = 1;
% Get the defaults from varargin
vararginoptions(varargin,{'mod','vox','interp','prefix'});

% Check for multiple files 
if (iscell(source)); source=char(source); end;
for i = 1:size(source,1)
        image = source(i,:);

    [source_dir,Sname,ext]=fileparts(image);

    if (isempty(source_dir))
        source_dir=pwd;
    end;
    

    suit_reslice(image,param,'mask',mask,...
        'bb',[-39 -65 -41;40 -6 -2],...
        'vox',vox,'prefix',prefix,'preserve',mod,'interp',interp);

    delete([source_dir '/m' Sname ext]);
    [mdir,mname]= spm_fileparts(mask);
    if isempty(mdir); mdir = pwd; end
    delete([mdir '/s' mname '.nii']);
end