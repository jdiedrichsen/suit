function suit_mni2suit(input,varargin)
% Apply a deformation field to an image to transform it from MNI space 
% to suit space or the other way around
%
% INPUT:
%   input:  Image to transform
%
% VARARGIN:
%   def:    'mni2suit' (default)   Direction of the trasnformation
%           'suit2mni'
%   mask:   Set to 1 to mask image before transformation, this creates a
%           new file with prefix 'm_' to be used as input.
%   Interp: Interpolation, default is nearest neighbour (0) 
%           set to 1 for trilinear interpolation 
% 
% OUTPUT:
%           Nifti image with prefix 'Wm2s_' for transformation from MNI 
%           to SUIT or Ws2m for transformation from SUIT to MNI.
%_______________________________________________________________________
% Carlos Hernandez 2018

def     = 'mni2suit';
mask    = 0;
Interp  = 0;
vararginoptions(varargin,{'def','mask','Interp'});

% Template Location
spm_Dir = fileparts(which('spm'));
templateDir = [spm_Dir '/toolbox/suit/templates'];

switch (def)
    case 'mni2suit'
        D = [templateDir '/def_mni2suit.nii'];
        opfx = 'Wm2s_';
        Tmask = [templateDir '/maskMNI.nii'];
    case 'suit2mni'
        D=[templateDir '/def_suit2mni.nii'];
        opfx = 'Ws2m_';
        Tmask = [templateDir '/maskSUIT.nii'];
end;

% Masking
if mask == 1
    M = spm_vol(Tmask);
    M = spm_read_vols(M);
    VI = spm_vol(input);
    Vi = spm_read_vols(VI);
    mVI = Vi.*M;
    save_vol(mVI,['m_' input],VI);
    In = ['m_' input];
else
    In = input;
end;

% Read Image and deformation
VI = spm_vol(In);
VD = spm_vol(D);

% Create Output
VO = VI;
for i=1:length(VO),
	VO(i).fname    = [opfx VO(i).fname];
	VO(i).dim(1:3) = VD(1).dim(1:3);
	VO(i).mat      = VD(1).mat;
	if ~isfield(VO,'descrip'), VO(i).descrip = ''; end;
	VO(i).descrip  = ['warped ' VO(i).descrip];
end;
VO = spm_create_vol(VO(i));

% Apply deformation
for p=1:VD(1).dim(3),
	M  = spm_matrix([0 0 p]);
	x1 = spm_slice_vol(VD(1), M, VD(1).dim(1:2),Interp);
	x2 = spm_slice_vol(VD(2), M, VD(1).dim(1:2),Interp);
	x3 = spm_slice_vol(VD(3), M, VD(1).dim(1:2),Interp);
	for i=1:length(VI),
		M     = inv(VI(i).mat);
		y1    = M(1,1)*x1+M(1,2)*x2+M(1,3)*x3+M(1,4);
		y2    = M(2,1)*x1+M(2,2)*x2+M(2,3)*x3+M(2,4);
		y3    = M(3,1)*x1+M(3,2)*x2+M(3,3)*x3+M(3,4);
		img   = spm_sample_vol(VI(i),y1,y2,y3,Interp);
		VO(i) = spm_write_plane(VO(i),img,p);
	end;
end;

end

function VR=save_vol(DATA,name,V)
VR    = struct(		'fname',name,...
    'dim',		[V.dim(1:3)],...
    'dt',   [spm_type('float32') 0],...
    'mat',		V.mat,...
    'pinfo',	[1 0 0]',...
    'descrip',	'temp');
VR    = spm_create_vol(VR);
spm_write_vol(VR,DATA);
end
