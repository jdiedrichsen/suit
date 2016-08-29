function suit_reslice_inv(source,param,varargin)
% function suit_reslice_inv(source,param,varname1,value1,...);
%   Reslice a set of volumes (source) from the SUIT template to the "original space"
%__________________________________________________________________________
% INPUT:
%   source: image(s) to be resliced: character/cell array 
%   param: mat file with the deformations from normalization
%   (mc_<name>_snc.mat).
%__________________________________________________________________________
% OUTPUT:
%       is written as i_<name> (or as prefix_<name> if prefix given)
%__________________________________________________________________________
% OPTIONS:
%   'interp': interpolation 0-> nearest neighbour 1-> tri-linear interpolation (default = 0 )
%   'prefix': prefix of output volume (default='i')
%   'reference': image that defines the dimension of output Volume (default: image given in _snc.mat)
%__________________________________________________________________________
% EXAMPLE:
%       suit_reslice_inv('<source>','mc_<name>_snc.mat');
% -------------------------------------------------------------------------
% Copyright (C) 2010 
% Tobias Wiestler, (Tobias.Wiestler@gmx.de) 
% Joern Diedrichsen (j.diedrichsen@ucl.ac.uk)
% Diedrichsen, J. (2006). A spatially unbiased atlas template of the human
% cerebellum. Neuroimage.

% v.2.4: Compatibility with SPM8
%        Ability to write multiple files 
% v.2.5 Job manager support in SPM8 (24/08/10)
% v.2.6 bug fixed initialized interp variable
% v.2.7 Dropped support for SPM2 and included support for SPM12b


outfilename=[];
reference=[];
interp = 0;

global defaults_suit; 
if (~isstruct(defaults_suit))
    suit_defaults;
end;
names=fieldnames(defaults_suit.reslice_inv); 
for i=1:length(names)
    eval([names{i} '= defaults_suit.reslice_inv.' names{i} ';']); 
end;


SCCSid   = '3.0';
SPMid    = spm('FnBanner',mfilename,SCCSid);

vararginoptions(varargin,{ 'interp', 'prefix','reference'});

% ----------------------------------------------------------------------
% Get the image and mask if necessary
% ----------------------------------------------------------------------
if (nargin<1 || isempty(source))
    source=spm_select([1 inf],'image','Get Source Image(s)');
end;
V_source=spm_vol(source);


if (nargin<2 || isempty(param))
    param=spm_select(1,'mat','Deformation map (*_snc.mat)');
end;

% ----------------------------------------------------------------------
% Make cells to char
% ----------------------------------------------------------------------
if (iscell(reference)) 
    reference=char(reference); 
end; 
if (iscell(param))
    param=char(param); 
end; 
if (iscell(source))
    source=char(source); 
end;

% ----------------------------------------------------------------------
% Open refence image if necessary
% ----------------------------------------------------------------------
if (isempty(reference))
    T=load(param);
    VRef=T.VF;
else
    if(isstruct(reference))
        VRef=reference;
    else
        VRef=spm_vol(reference);
    end;
end;

%-------------------------------------------------------------------
% get the deformation map in x y z. tells you the position in the original
% image as a function of position in suit space & invert 
%-------------------------------------------------------------------
[def,defMat] = spmdefs_get_sn2def(param);
[i_def,i_defMat] = spmdefs_get_inv(def,defMat,VRef); % invert the transformation


%-------------------------------------------------------------------
% apply the deformation to the source image
%-------------------------------------------------------------------
target_dir=spm_fileparts(param); 

for i=1:size(source,1)
    [source_dir,sourceName,ext] = spm_fileparts(source(i,:));

    outfilename= fullfile(target_dir,[prefix, sourceName, ext]);
    spmdefs_apply_def(i_def,i_defMat,V_source(i).fname,interp,outfilename);
end; 


%-------------------------------------------------------------------
function [Def,mat] = spmdefs_get_sn2def(matname,bb,vox)
% function [Def,mat] = spmdefs_get_sn2def(matname,bb,vox)
% Convert a SPM _sn.mat file into a deformation field, and return it.
% INPUT:
%   matname: Name of mat file
%   bb:     Bounding box (default: atlas space of deformation)
%   vox:    voxel size (default: atlas space of deformation)
% OUTPUT:
%   Def:    Nonlinear parts
%   mat:    Affine transformation
sn  = load(matname);
if nargin<3
    [bb,vox]=spmdefs_bbvox_from_V(sn.VG(1));
end;

[bb0,vox0] = spmdefs_bbvox_from_V(sn.VG(1));

if any(~isfinite(vox)), vox = vox0; end;
if any(~isfinite(bb)),  bb  = bb0;  end;
bb  = sort(bb);
vox = abs(vox);

% Adjust bounding box slightly - so it rounds to closest voxel.
bb(:,1) = round(bb(:,1)/vox(1))*vox(1);
bb(:,2) = round(bb(:,2)/vox(2))*vox(2);
bb(:,3) = round(bb(:,3)/vox(3))*vox(3);

M   = sn.VG(1).mat;
vxg = sqrt(sum(M(1:3,1:3).^2));
ogn = M\[0 0 0 1]';
ogn = ogn(1:3)';

% Convert range into range of voxels within template image
x   = (bb(1,1):vox(1):bb(2,1))/vxg(1) + ogn(1);
y   = (bb(1,2):vox(2):bb(2,2))/vxg(2) + ogn(2);
z   = (bb(1,3):vox(3):bb(2,3))/vxg(3) + ogn(3);

og  = -vxg.*ogn;
of  = -vox.*(round(-bb(1,:)./vox)+1);
M1  = [vxg(1) 0 0 og(1) ; 0 vxg(2) 0 og(2) ; 0 0 vxg(3) og(3) ; 0 0 0 1];
M2  = [vox(1) 0 0 of(1) ; 0 vox(2) 0 of(2) ; 0 0 vox(3) of(3) ; 0 0 0 1];
mat = sn.VG(1).mat*inv(M1)*M2;
% dim = [length(x) length(y) length(z)];

[X,Y] = ndgrid(x,y);

st = size(sn.Tr);

if (prod(st) == 0),
    affine_only = true;
    basX = 0;
    basY = 0;
    basZ = 0;
else
    affine_only = false;
    basX = spm_dctmtx(sn.VG(1).dim(1),st(1),x-1);
    basY = spm_dctmtx(sn.VG(1).dim(2),st(2),y-1);
    basZ = spm_dctmtx(sn.VG(1).dim(3),st(3),z-1);
end,

Def = single(0);
Def(numel(x),numel(y),numel(z)) = 0;
Def = {Def; Def; Def};

for j=1:length(z)
    if (~affine_only)
        tx = reshape( reshape(sn.Tr(:,:,:,1),st(1)*st(2),st(3)) *basZ(j,:)', st(1), st(2) );
        ty = reshape( reshape(sn.Tr(:,:,:,2),st(1)*st(2),st(3)) *basZ(j,:)', st(1), st(2) );
        tz = reshape( reshape(sn.Tr(:,:,:,3),st(1)*st(2),st(3)) *basZ(j,:)', st(1), st(2) );

        X1 = X    + basX*tx*basY';
        Y1 = Y    + basX*ty*basY';
        Z1 = z(j) + basX*tz*basY';
    end

    Mult = sn.VF.mat*sn.Affine;
    if (~affine_only)
        X2= Mult(1,1)*X1 + Mult(1,2)*Y1 + Mult(1,3)*Z1 + Mult(1,4);
        Y2= Mult(2,1)*X1 + Mult(2,2)*Y1 + Mult(2,3)*Z1 + Mult(2,4);
        Z2= Mult(3,1)*X1 + Mult(3,2)*Y1 + Mult(3,3)*Z1 + Mult(3,4);
    else
        X2= Mult(1,1)*X + Mult(1,2)*Y + (Mult(1,3)*z(j) + Mult(1,4));
        Y2= Mult(2,1)*X + Mult(2,2)*Y + (Mult(2,3)*z(j) + Mult(2,4));
        Z2= Mult(3,1)*X + Mult(3,2)*Y + (Mult(3,3)*z(j) + Mult(3,4));
    end

    Def{1}(:,:,j) = single(X2);
    Def{2}(:,:,j) = single(Y2);
    Def{3}(:,:,j) = single(Z2);
end;



%-------------------------------------------------------------------
function [Def,mat] = spmdefs_get_inv(Def0,mat0,VT)
% function [Def,mat] = spmdefs_get_inv(Def0,mat0,VT)
% Invert a deformation field (derived from a composition of deformations)
M0      = mat0;
M1      = inv(VT.mat);
M0(4,:) = [0 0 0 1];
M1(4,:) = [0 0 0 1];
[Def{1},Def{2},Def{3}]    = spm_invdef(Def0{:},VT.dim(1:3),M1,M0);
mat         = VT.mat;




%-------------------------------------------------------------------
function spmdefs_apply_def(Def,mat,fnames,intrp,ofnames)
% function spmdefs_apply_def(Def,mat,fnames,intrp,ofnames)
% Warp an image or series of images according to a deformation field
intrp = [intrp*[1 1 1], 0 0 0];

for i=1:size(fnames,1),
    V = spm_vol(fnames(i,:));
    M = inv(V.mat);
    [pth,nam,ext] = spm_fileparts(fnames(i,:));
    if (nargin<5)
        ofname = fullfile(pth,['w',nam,ext]);
    else
        ofname= deblank(ofnames(i,:));
    end;
    Vo = struct('fname',ofname,...
        'dim',[size(Def{1},1) size(Def{1},2) size(Def{1},3)],...
        'dt',V.dt,...
        'pinfo',V.pinfo,...
        'mat',mat,...
        'n',V.n,...
        'descrip',V.descrip);
    C  = spm_bsplinc(V,intrp);
    Vo = spm_create_vol(Vo);
    for j=1:size(Def{1},3)
        d0    = {double(Def{1}(:,:,j)), double(Def{2}(:,:,j)),double(Def{3}(:,:,j))};
        d{1}  = M(1,1)*d0{1}+M(1,2)*d0{2}+M(1,3)*d0{3}+M(1,4);
        d{2}  = M(2,1)*d0{1}+M(2,2)*d0{2}+M(2,3)*d0{3}+M(2,4);
        d{3}  = M(3,1)*d0{1}+M(3,2)*d0{2}+M(3,3)*d0{3}+M(3,4);
        dat   = spm_bsplins(C,d{:},intrp);
        Vo    = spm_write_plane(Vo,dat,j);
    end;
end;
return;

%-------------------------------------------------------------------
function [bb,vx] = spmdefs_bbvox_from_V(V)
% Return the default bounding box for an image volume

vx = sqrt(sum(V.mat(1:3,1:3).^2));
o  = V.mat\[0 0 0 1]';
o  = o(1:3)';
bb = [-vx.*(o-1) ; vx.*(V.dim(1:3)-o)];
return;
