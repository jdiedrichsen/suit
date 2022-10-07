function suit_save_darteldef(name,varargin)
% function suit_save_darteldef(varargin);
% Saves a dartel flowfield as a deformation map that can be used for
% subsequent normalization using SUITPy 
% _________________________________________________________________________
% INPUT: 
%   name:           name of the main file 
%      Expects to find the following two output files from suit_normalize_dartel: 
%       Affine_<name>.mat: Matrix that contains linear part of the mapping 
%       u_a_<name>.nii   : Flow field that contains the nonlinear part of the mapping
%__________________________________________________________________________
% OPTIONAL: 
%   wdir: working Directory name (default:pwd)
%__________________________________________________________________________
% OUTPUT:
%    saves the following file" 
%       y_<name>_suitdef.nii: Deformation image 
% ----------------------------------------------------------------------
% Typical usuage:
% suit_save_darteldef('sub01-anatomical','wdir','/usr/name/projectdir')
% ----------------------------------------------------------------------
% Joern Diedrichsen 2022 joern.diedrichsen@googlemail.com


% Deal with variable arguments:
wdir=[]; 
vararginoptions(varargin,{'wdir'});

if isempty(wdir)
    wdir= pwd(); 
end

affine = fullfile(wdir,['Affine_c_' name '_seg1.mat']);
flowfield = fullfile(wdir,['u_a_c_' name '_seg1.nii']); 

[Def,mat]=spmdefs_get_dartel(flowfield,affine); 
spmdefs_save_def(Def,mat,fullfile(wdir,['y_' name '_suitdef.nii'])); 
   

function spmdefs_save_def(Def,mat,fname)
% function spmdefs_save_def(Def,mat,ofname)
% Saves a deformation field as a nii-file 
% INPUT: 
%   Def: deformation field 
%   mat: transformation matrix 
%   fname: outfilename
% 

dim   = [size(Def{1},1) size(Def{1},2) size(Def{1},3) 1 3]; % was 1 3 
dtype = 'FLOAT32-BE';
off   = 0;
scale = 1;
inter = 0;
dat   = file_array(fname,dim,dtype,off,scale,inter);

N      = nifti;
N.dat  = dat;
N.mat  = mat;
N.mat0 = mat;
N.mat_intent  = 'Aligned';
N.mat0_intent = 'Aligned';
N.intent.code = 'VECTOR';
N.intent.name = 'Mapping';
N.descrip = 'Deformation field';
create(N);
N.dat(:,:,:,1,1) = Def{1}; % was 1,1 
N.dat(:,:,:,1,2) = Def{2}; % was 1,2 
N.dat(:,:,:,1,3) = Def{3}; % was 1,3 
