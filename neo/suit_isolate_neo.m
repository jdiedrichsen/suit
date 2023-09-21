function suit_isolate_neo(source,varargin)
% Neonatal cerebellar isolation tool. Resulting mask might require manual
% correction. 
%
% Input:            T2 weighted anatomical image. Use cell array of char
%                   for multiple images.
%
% Options:
%    maskp:         Probability threshold for mask, higher values will
%                   result in a more conservative mask. (defualt 0.5)
%    getV:          Estimates the cerebellar-, cortical- and Intracraneal
%                   volume and save it in a mat file. Volume must be
%                   recalculated after correction of cerebellar mask or any
%                   of the mask provided using the function suit_vol
%    Keeptfiles:    Set to 1 to do not erase temporal files
%
% Output:
%   c_<image>           Cropped source image
%   c_<image>_pcereb    Isolation mask image
%
%   getV option:
%   cort_<name>         Cortical mask image
%   ICV_<name>          Intracraneal volume mask image
%   vol_<name>          Volume estimate .mat file
%
% Example:
%
%   Default         suit_isolate_neo({'<name>.nii'})
%
%   Volume masks    suit_isolate_neo({'<name>.nii'},'getV',1)
% _______________________________________________________________________
% Carlos R. Hernandez-Castillo 2018

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

names=fieldnames(defaults_suit.isolateNeo); 
for i=1:length(names)
    eval([names{i} '= defaults_suit.isolateNeo.' names{i} ';']); 
end;


% -----------------------------------------------------------------
% Check spm-version
% -----------------------------------------------------------------
spmVer=spm('Ver');
if (strcmp(spmVer,'SPM5')||strcmp(spmVer,'SPM8'))
   error('Install SPM 12 or newer version to use SUIT_isolate_seg');
end

% Get the defaults from varargin
vararginoptions(varargin,{'maskp','keeptfiles','getV'});

% Check for multiple files
if (iscell(source)); source=char(source); end;
for i = 1:size(source,1)
        image = source(i,:);

    [source_dir,Sname,ext,num]=spm_fileparts(image);

    if (isempty(source_dir))
        source_dir=pwd;
    end;

    % ------------------------------------------------------------------
    % Segmentation
    % ------------------------------------------------------------------
        J.channel(1).vols = {image};
        J.channel(1).biasreg = 0.001;
        J.channel(1).biasfwhm = 60;
        J.channel(1).write = [0 0];
        J.tissue(1).tpm = {[priors ',1']};
        J.tissue(1).ngaus = 1;
        J.tissue(1).native = [1 0];
        J.tissue(1).warped = [0 0];
        J.tissue(2).tpm = {[priors ',2']};
        J.tissue(2).ngaus = 1;
        J.tissue(2).native = [1 0];
        J.tissue(2).warped = [0 0];
        J.tissue(3).tpm = {[priors ',3']};
        J.tissue(3).ngaus = 2;
        J.tissue(3).native = [1 0];
        J.tissue(3).warped = [0 0];
        J.tissue(4).tpm = {[priors ',4']};
        J.tissue(4).ngaus = 4;
        J.tissue(4).native = [0 0];
        J.tissue(4).warped = [0 0];
        J.tissue(5).tpm = {[priors ',5']};
        J.tissue(5).ngaus = 3;
        J.tissue(5).native = [0 0];
        J.tissue(5).warped = [0 0];
        J.tissue(6).tpm = {[priors ',6']};
        J.tissue(6).ngaus = 1;
        J.tissue(6).native = [1 0];
        J.tissue(6).warped = [0 0];
        J.tissue(7).tpm = {[priors ',7']};
        J.tissue(7).ngaus = 1;
        J.tissue(7).native = [1 0];
        J.tissue(7).warped = [0 0];

        J.warp.mrf = 1;
        J.warp.cleanup = 1;
        J.warp.reg = [0 0.001 0.5 0.05 0.2];
        J.warp.affreg = 'mni';
        J.warp.fwhm = 0;
        J.warp.samp = 3;
        J.warp.write = [0 1];

        spm_preproc_run(J,'run'); 

    if getV == 1
    % -----------------------------------------------------------------
    % Calculate volumes
    % -----------------------------------------------------------------
    spm_imcalc({[source_dir,'/c1',Sname,ext],[source_dir,'/c2',Sname,ext]},...
        [source_dir,'/cere_',Sname,ext],'(i1+i2)>0.5',{[],[],[],4});
    spm_imcalc({[source_dir,'/c6',Sname,ext],[source_dir,'/c7',Sname,ext]},...
        [source_dir,'/cort_',Sname,ext],'(i1+i2)>0.5',{});
    spm_imcalc({[source_dir,'/c1',Sname,ext],[source_dir,'/c2',Sname,ext],...
        [source_dir,'/c3',Sname,ext],[source_dir,'/c6',Sname,ext],...
        [source_dir,'/c7',Sname,ext]},[source_dir,'/ICV_',Sname,ext],'(i1+i2+i3+i4+i5)>0.5',{});

    V1 = spm_vol([source_dir,'/cere_',Sname,ext]);
    V2 = spm_vol([source_dir,'/cort_',Sname,ext]);
    V3 = spm_vol([source_dir,'/ICV_',Sname,ext]);
    
    SVoxSize = V1.mat(1,1) * V1.mat(2,2) * V1.mat(3,3);

    v1 = spm_read_vols(V1);
    v2 = spm_read_vols(V2);
    v3 = spm_read_vols(V3);

    cereVox = numel(find(single(v1)==1));
    cortVox = numel(find(single(v2)==1));
    ICVVox = numel(find(single(v3)==1));
    cereMM = cereVox * SVoxSize;
    cortMM = cortVox * SVoxSize;
    ICVMM = ICVVox * SVoxSize;
    delete([source_dir,'/cere_',Sname,ext]);
    save([source_dir,'/vol_',Sname,'.mat'],'cereVox','cortVox','ICVVox','cereMM','cortMM','ICVMM','SVoxSize');
    fprintf('Volumes saved in %s\nIf any mask is corrected, use suit_vol to recalculate the volume\n',...
        [source_dir,'/vol_',Sname,'.mat']);
    end

    % -----------------------------------------------------------------
    % Crop the volumes to a standard bounding box
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
    % Creating final mask
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
    % Cleaning
    % -----------------------------------------------------------------
    if (keeptfiles==0)
        rm_imgfile([source_dir,'/c1',Sname],ext);
        rm_imgfile([source_dir,'/c2',Sname],ext);
        rm_imgfile([source_dir,'/c3',Sname],ext);
        rm_imgfile([source_dir,'/c6',Sname],ext);
        rm_imgfile([source_dir,'/c7',Sname],ext);
        rm_imgfile([source_dir,'/c_',Sname,'_seg1'],ext);
        rm_imgfile([source_dir,'/c_',Sname,'_seg2'],ext);
        rm_imgfile([source_dir,'/',Sname,'_seg8'],'.mat');
        rm_imgfile([source_dir,'/y_',Sname],ext);
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
end;
end
