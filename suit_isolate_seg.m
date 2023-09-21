function suit_isolate_seg(Source,varargin)
% Function to isolate the cerebellum and brainsteam using the spm12
% segmentation to create a cerebellum mask for subsecuent normalization to  
% a specialized template. This function is part of SUIT.
%__________________________________________________________________________
% INPUT:
%       Cell array of Whole brain anatomical scan (T1)
%       Additional channels (T2, PD, etc.) can be used to improve
%       isolation, In order to add more channels include them as a
%       different element in the cell array.
%       All channels must be corregistered and resliced to the T1.
%__________________________________________________________________________
% OUTPUT:
%   c_<Source>
%       The cropped anatomical containing the cerebellum
%   c_<Source>_pcereb
%       The binarized mask after thresholded at p=0.2 for hand correction
%   <Source>_seg1 - <Source>_seg2
%       Probabilities maps of Gray matter(seg1) or White matter(seg2)
%__________________________________________________________________________
% OPTIONS
%   maskp           Probability value for final mask higher value 
%                   returns a tighter mask. defualt 0.2
%   bb              Use a bounding box different from default feed as 
%                   [-x,x;-y,y;-z,z] in mm in MNI space.
%   keeptempfiles   set to 1 to keep temporal files
%                   If selected the probability maps will be stored as
%                   c*<Source>.nii, where c1 = cerebellar gray matter, 
%                   c2 = cerebellar white matter, c7 = cortical gray matter
%                   c8 = cortical white matter.
%__________________________________________________________________________
% Carlos Hernandez-Castillo 2016
% SUIT Copyright (C) 2010 
% Use of the SUIT software should be cited as:
% Diedrichsen, J. (2006). A spatially unbiased atlas template of the human
% cerebellum. Neuroimage.

% ------------------------------------------------------------------------
% Set the defaults.
% ------------------------------------------------------------------------

global defaults;
if (~isstruct(defaults))
    error('Start SPM to use SUIT_isolate');
end;

global defaults_suit; 
if (~isstruct(defaults_suit))
    suit_defaults;
end;

names=fieldnames(defaults_suit.isolate_seg); 
for i=1:length(names)
    eval([names{i} '= defaults_suit.isolate_seg.' names{i} ';']); 
end;

% -----------------------------------------------------------------
% check spm-version
% -----------------------------------------------------------------
spmVer=spm('Ver');
if (strcmp(spmVer,'SPM5')||strcmp(spmVer,'SPM8'))
   error('Install SPM 12 or newer version to use SUIT_isolate_seg');
end

% -----------------------------------------------------------------
% Get Source image
% -----------------------------------------------------------------
if (nargin<1 || isempty(Source))
    spmIn=spm_select(Inf,'image','Get Source Image(s)');
    for n=1:size(spmIn,1)
        Source{n}=spmIn(n,:);
    end
end;

if (~isstruct(Source)) % Check for Job structure or single file

[source_dir,Sname,ext,~]=spm_fileparts(Source{1});
    
if (isempty(source_dir))
    source_dir=pwd;
end;

% Get the defaults from varargin
vararginoptions(varargin,{'maskp','keeptempfiles','bb'});

% ------------------------------------------------------------------
% Segmentation
% ------------------------------------------------------------------
    J=[]; 
    for chan=1:length(Source)
        volumen=Source{chan};
        J.channel(chan).vols = {volumen};
        J.channel(chan).biasreg = 0.001;
        J.channel(chan).biasfwhm = 60;
        J.channel(chan).write = [0 0];
    end
    J.tissue(1).tpm = {[prior_dir,'/',priors,',1']};
    J.tissue(1).ngaus = 1;
    J.tissue(1).native = [1 0];
    J.tissue(1).warped = [0 0];
    J.tissue(2).tpm = {[prior_dir,'/',priors,',2']};
    J.tissue(2).ngaus = 1;
    J.tissue(2).native = [1 0];
    J.tissue(2).warped = [0 0];
    J.tissue(3).tpm = {[prior_dir,'/',priors,',3']};
    J.tissue(3).ngaus = 2;
    J.tissue(3).native = [1 0];
    J.tissue(3).warped = [0 0];
    J.tissue(4).tpm = {[prior_dir,'/',priors,',4']};
    J.tissue(4).ngaus = 4;
    J.tissue(4).native = [0 0];
    J.tissue(4).warped = [0 0];
    J.tissue(5).tpm = {[prior_dir,'/',priors,',5']};
    J.tissue(5).ngaus = 3;
    J.tissue(5).native = [0 0];
    J.tissue(5).warped = [0 0];
    J.tissue(6).tpm = {[prior_dir,'/',priors,',6']};
    J.tissue(6).ngaus = 3;
    J.tissue(6).native = [0 0];
    J.tissue(6).warped = [0 0];
    J.tissue(7).tpm = {[prior_dir,'/',priors,',7']};
    J.tissue(7).ngaus = 1;
    J.tissue(7).native = [1 0];
    J.tissue(7).warped = [0 0];
    J.tissue(8).tpm = {[prior_dir,'/',priors,',8']};
    J.tissue(8).ngaus = 1;
    J.tissue(8).native = [1 0];
    J.tissue(8).warped = [0 0];

    J.warp.mrf = 1;
    J.warp.cleanup = 1;
    J.warp.reg = [0 0.001 0.5 0.05 0.2];
    J.warp.affreg = 'mni';
    J.warp.fwhm = 0;
    J.warp.samp = 3;
    J.warp.write = [0 1];

    spm_preproc_run(J,'run'); 
  
% -----------------------------------------------------------------
% crop the volumes to a standard bounding box
% -----------------------------------------------------------------
    fprintf('Cropping...');
    [defy,maty]=spmdefs_get_def([source_dir,'/y_',Sname,ext]);
    bbox=[bb(1,1),bb(2,1),bb(3,1);bb(1,2),bb(2,2),bb(3,2);...
            bb(1,2),bb(2,1),bb(3,1);bb(1,1),bb(2,2),bb(3,2);...
            bb(1,2),bb(2,2),bb(3,1);bb(1,1),bb(2,2),bb(3,1);...
            bb(1,1),bb(2,1),bb(3,2);bb(1,2),bb(2,1),bb(3,2)];
    BOUND=spmdefs_transformM(defy,maty,bbox);
    V=spm_vol([source_dir,'/',Sname,ext]);
    BOUND=round((BOUND' ./ repmat(diag(V.mat(1:3,1:3)),1,8))-...
        (repmat(V.mat(1:3,4),1,8) ./ repmat(diag(V.mat(1:3,1:3)),1,8)));
    cropped_bb.x = min(BOUND(1,:)):max(BOUND(1,:));
    cropped_bb.y = min(BOUND(2,:)):max(BOUND(2,:));
    cropped_bb.z = min(BOUND(3,:)):max(BOUND(3,:));

    crop_image=[source_dir,'/c_',Sname,ext];
    crop_seg1=[source_dir,'/c_',Sname,'_seg1',ext];
    crop_seg2=[source_dir,'/c_',Sname,'_seg2',ext];
    try
        spmj_crop_vol([source_dir,'/',Sname,ext],crop_image,cropped_bb.x,cropped_bb.y,cropped_bb.z);
        spmj_crop_vol([source_dir,'/c1',Sname,ext],crop_seg1,cropped_bb.x,cropped_bb.y,cropped_bb.z);
        spmj_crop_vol([source_dir,'/c2',Sname,ext],crop_seg2,cropped_bb.x,cropped_bb.y,cropped_bb.z);
    catch
        error(['Cropping failed.','Make sure that the original image is in LPI orientation.']);
    end
    
    fprintf(' done\n');
 
% -----------------------------------------------------------------
% creating final mask
% -----------------------------------------------------------------
s1 = spm_vol(crop_seg1);
s2 = spm_vol(crop_seg2);
S1 = spm_read_vols(s1);
S2 = spm_read_vols(s2);
M = S1+S2;
M(M < maskp) = 0;
M(M > 0) = 1;
save_vol(M,[source_dir,'/c_',Sname,'_pcereb',ext],s1);
 
% -----------------------------------------------------------------
% cleaning
% -----------------------------------------------------------------
if (keeptempfiles==0)
    movefile([source_dir,'/c1',Sname,ext],[source_dir,'/',Sname,'_seg1',ext]);
    movefile([source_dir,'/c2',Sname,ext],[source_dir,'/',Sname,'_seg2',ext]);
    rm_imgfile([source_dir,'/c3',Sname],ext);
    rm_imgfile([source_dir,'/c7',Sname],ext);
    rm_imgfile([source_dir,'/c8',Sname],ext);
    rm_imgfile([source_dir,'/c_',Sname,'_seg1'],ext);
    rm_imgfile([source_dir,'/c_',Sname,'_seg2'],ext);
    rm_imgfile([source_dir,'/',Sname,'_seg8'],'.mat');
    rm_imgfile([source_dir,'/y_',Sname],ext);

end;
        

else % loop for Job strucure
    s=Source.source;
    for i=1:length(s)
        suit_isolate_seg(s{i},'maskp',Source.maskp,'keeptempfiles',Source.keeptempfiles,'bb',Source.bb); 
    end;  
end

end

function rm_imgfile(name,ext)
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

    
end

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
    