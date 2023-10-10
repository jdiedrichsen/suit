function suit_reslice_dartel_inv(job)
%   Reslice a set of volumes from SUIT space to the native space using
%   a flow field and affine transformation previouly calcualted with
%   suit_normalize_dartel.
%__________________________________________________________________________
% INPUT:
%   job: Structure with fields:
% 
%       job.Affine:     Affine Transform parameter file (Affine*.mat)
%       job.flowfield:  Flowfield file (u_a*.nii)
%       job.resample:   Images that need to be resliced. If this field
%                       doesn't exist, It will reslice the SUIT Lobular 
%                       atlas into the native space.
%       job.ref         Reference image for final geometry of the ouput. 
%                       e.g. initial T1. If this field doesn't exist, final
%                       image will have the same size than SUIT template.
%__________________________________________________________________________
% OUTPUT:
%   iw_<source_name>:   The resliced source images.
%__________________________________________________________________________
% Carlos Hernandez 2018
%

% Directiores
global defaults_suit;
if isempty(defaults_suit)
    suit_defaults;
end

spmDir   =  fileparts(which('spm'));
LobAtlas =  defaults_suit.summarize.atlas{1};

% Check the resample input
if (~isfield(job,'resample'))
    if (~exist(LobAtlas,'file'))
        error(sprintf('Atlas file: %s not found. \nYou may have to download github/DiedrichsenLab/cerebellar_atlases, \n or set the location of the atlas directory in suit_defaults.m',LobAtlas));
    end
    job.resample={LobAtlas};
end

% Load the Affine Matrix
load(char(job.Affine));

% Get the dartel deformation
times  = [1 0]; 
K      = 6; 
N      = nifti(job.flowfield);
y      = spm_dartel_integrate(N.dat,times,K);
Def    = cell(3,1);
M      = single(inv(Affine)*N.mat);
mat    = N.mat0;
Def{1} = y(:,:,:,1)*M(1,1) + y(:,:,:,2)*M(1,2) + y(:,:,:,3)*M(1,3) + M(1,4);
Def{2} = y(:,:,:,1)*M(2,1) + y(:,:,:,2)*M(2,2) + y(:,:,:,3)*M(2,3) + M(2,4);
Def{3} = y(:,:,:,1)*M(3,1) + y(:,:,:,2)*M(3,2) + y(:,:,:,3)*M(3,3) + M(3,4);

% Check if reference image was provided
if (~isfield(job,'ref'));job.ref={LobAtlas};end

% Inverse deformation
VT = spm_vol(char(job.ref));
M0      = mat;
M1      = inv(VT.mat);
M0(4,:) = [0 0 0 1];
M1(4,:) = [0 0 0 1];
[invDef{1},invDef{2},invDef{3}] = spm_invdef(Def{:},VT.dim(1:3),M1,M0);
invmat         = VT.mat;

% Apply inverse transformation
for i=1:size(job.resample,1),
    V = spm_vol(char(job.resample{i}));
    M = inv(V.mat);
    [~,nam,ext] = spm_fileparts(char(job.resample{i}));
    [pth,naff,~] = spm_fileparts(char(job.flowfield));
    ofname = fullfile(pth,['iw_',nam,'_',naff,ext]);

    Vo = struct('fname',ofname,...
                'dim',[size(invDef{1},1) size(invDef{1},2) size(invDef{1},3)],...
                'dt',V.dt,...
                'pinfo',V.pinfo,...
                'mat',invmat,...
                'n',V.n,...
                'descrip',V.descrip);
    C  = spm_bsplinc(V,[0 0 0 0 0 0]);
    Vo = spm_create_vol(Vo);
    for j=1:size(invDef{1},3)
        d0    = {double(invDef{1}(:,:,j)), double(invDef{2}(:,:,j)),double(invDef{3}(:,:,j))};
        d{1}  = M(1,1)*d0{1}+M(1,2)*d0{2}+M(1,3)*d0{3}+M(1,4);
        d{2}  = M(2,1)*d0{1}+M(2,2)*d0{2}+M(2,3)*d0{3}+M(2,4);
        d{3}  = M(3,1)*d0{1}+M(3,2)*d0{2}+M(3,3)*d0{3}+M(3,4);
        dat   = spm_bsplins(C,d{:},[0 0 0 0 0 0]);
        Vo    = spm_write_plane(Vo,dat,j);
    end;
end;

end
    