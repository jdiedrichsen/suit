function suit_reslice_dartel(job,varargin);
% function suit_reslice_dartel(job,varargin...);
%   Reslice a set of volumes (source) to the SUIT template using
%   normalization of a dartel flow field.
%   We recommend that you use a mask to mask out non-cerebellar structures
%   before this. For functional data from the deep cerebellar nuclei it is
%   highly recommended that you use the dentate / interposed mask to mask
%   the images.
%__________________________________________________________________________
% INPUT:
%   job:                Structure with fields:
%       job.subj:           A Structure for each subject (individual flowfields)
%       job.subj.affineTr:  Affine Transform parameter file (Affine*.mat)
%       job.subj.flowfield: Flowfield file (u_*)
%       job.subj.resample:  Images that need to be resliced
%       job.subj.mask:      Masking file that is applied before reslicing
%__________________________________________________________________________
% OUTPUT:
%   w<source_name>:     The resliced source images.
%__________________________________________________________________________
% OPTIONS:
%   job.jactransf:      Jacobian volume modulation 0:no, for fMRI, VBQ data
%                       1: yes, use for VBM
%   job.K:              Number of time steps (default 6)
%   job.bb:             Bounding box:  2x3 matrix
%   job.vox:            voxel resolution (default 1 1 1)
%   job.interp:         Interpolation: 0; nearest neighbour, 1: trilinear, etc
%   job.prefix:         prefix for resliced images (default 'wd'
%   job.subj.outname:   Outnames - otherwise uses prefix to
%                       determine outname 
%__________________________________________________________________________
% EXAMPLE:
%         for s=1
%             job.subj(s).affineTr={['Affine_<name>' subjNames{s} '_seg1.mat']};       % where subjName is a cell array of subject names
%             job.subj(s).flowfield={['u_a_<name>' subjNames{s} '_seg1.nii']};       
%             job.subj(s).resample{1}={['<name>' subjNames{s} '.nii']}; % Cell array containing all images to be resliced.
%             job.subj(s).mask={['c_struct<name>' subjNames{s} '_pcereb_corr.nii']}; % Mask to be applied 
%         end; 
%        suit_reslice_dartel(job); 
% ----------------------------------------------------------------------

% Joern Diedrichsen 8/4/2013 j.diedrichsen@ucl.ac.uk
% v.2.6: First version
% v.2.7: Bug fixes for removal of temporary files and output files  
% v.3.0: No changes 

global defaults;
global defaults_suit;
if (~isstruct(defaults_suit))
    suit_defaults;
end;


% Add parameters from the defaults if not specified in job structure
fields={'K','bb','vox','interp','prefix','jactransf'}; 
for i=1:length(fields) 
    if (~isfield(job,fields{i}))
        job.(fields{i})=defaults_suit.reslice_dartel.(fields{i}); 
    end; 
end; 


SCCSid   = '3.0';
SPMid    = spm('FnBanner',mfilename,SCCSid);

% Loop over subjects
for s=1:length(job.subj)
    S=job.subj(s);

    % --------------------------------------------------
    % Check and load the affine transformation matrix 
    if (iscell(S.affineTr))
        S.affineTr=char(S.affineTr); 
    end; 
    load(S.affineTr);
    if (~exist('Affine','var'))
        error('.mat file does not contain an affine transformation matrix'); 
    end; 
    
    % --------------------------------------------------
    % get and load the flowfield 
    if (iscell(S.flowfield))
            S.flowfield=S.flowfield{1};
    end;
    [flow_dir,flow_name,flow_ext,flow_num]=spm_fileparts(S.flowfield); 
    Vff=spm_vol(S.flowfield); 
    % Generate grid in the space of the template 
    [X,Y,Z]=ndgrid(1:Vff.dim(1),1:Vff.dim(2),1:Vff.dim(3));
    num_slice=Vff.dim(3);
    
    %---------------------------------------------------
    % If necessary load mask into the space of 
    if (~isempty(S.mask))
        if (iscell(S.mask))
            S.mask=char(S.mask); 
        end; 
        VM=spm_vol(S.mask); 
        [Xm,Ym,Zm]=spmj_affine_transform(X,Y,Z,inv(Affine*VM.mat)*Vff.mat);
        for z=1:num_slice
            MaskData(:,:,z)=spm_sample_vol(VM,Xm(:,:,z),Ym(:,:,z),Zm(:,:,z),job.interp);
        end;
        MaskData=double(MaskData>0);
    else 
        MaskData=ones(size(X)); 
    end; 
    
    %---------------------------------------------------
    % Resample all images into affine space first
    N=length(S.resample);
    for i=1:N
        [image_dir,image_name,image_ext,image_num]=spm_fileparts(S.resample{i}); 
        V(i)=spm_vol(S.resample{i}); 
        [Xm,Ym,Zm]=spmj_affine_transform(X,Y,Z,inv(Affine*V(i).mat)*Vff.mat);
        for z=1:num_slice
            Data(:,:,z)=spm_sample_vol(V(i),Xm(:,:,z),Ym(:,:,z),Zm(:,:,z),job.interp);
        end;
        Data=Data.*MaskData; % Mask the image
        O(i)=V(i); 
        O(i).dim=Vff.dim; 
        O(i).mat=Vff.mat;
        J.images{i}{1}=fullfile(image_dir,['atemp_' image_name image_ext]);
        if (~isfield(S,'outname'))
            finalname{i}=fullfile(image_dir,[job.prefix image_name image_ext]); 
        else
            finalname{i}=S.outname{i}; 
        end;
        O(i).fname=J.images{i}{1}; 
        if (O(i).dt(1)==2)
            O(i).pinfo=[1./254*max(Data(:));0;O(i).pinfo(3)]; 
        end; 
        spm_write_vol(O(i),Data); 
    end;
    
    % -----------------------------------------------------
    % Then do the nonlinear deformation using spm_dartel_norm
    J.flowfields=job.subj(s).flowfield; 
    J.interp=job.interp;
    J.K=job.K; 
    J.jactransf=job.jactransf; 
    out=spm_dartel_norm(J); 
    
    % -----------------------------------------------------
    % Now load and reslice the images into the required space
    % Generate wished grid 
    [X,Y,Z]=ndgrid(job.bb(1,1):job.vox(1):job.bb(2,1),...
                   job.bb(1,2):job.vox(2):job.bb(2,2),...
                   job.bb(1,3):job.vox(3):job.bb(2,3)); % Eucledian coordinates 
    [Xm,Ym,Zm]=spmj_affine_transform(X,Y,Z,inv(Vff.mat));% voxel coordinates in images 
    if (isempty(Xm))
        error('bb(1,i) must be smaller than bb(2,i) for all i={1..3}');
    end; 
    
    % Now get the output images into the right resolution and bounding box 
    for i=1:N  
        V=spm_vol(out.files{i});
        Data=zeros(size(Xm)); 
        for z=1:size(Xm,3); 
            Data(:,:,z)=spm_sample_vol(V,Xm(:,:,z),Ym(:,:,z),Zm(:,:,z),job.interp);
        end;
        V.fname=finalname{i};
        V.dim=size(Xm); 
        % Determine the transformation matrix from 4 points
        % This is safe, but can probably be done more elegantly 
        i1=V.dim(1); 
        i2=V.dim(2); 
        i3=V.dim(3); 
        x=[[X(1,1,1);Y(1,1,1);Z(1,1,1);1],...
            [X(i1,1,1);Y(i1,1,1);Z(i1,1,1);1],...
            [X(i1,i2,1);Y(i1,i2,1);Z(i1,i2,1);1],...
            [X(i1,i2,i3);Y(i1,i2,i3);Z(i1,i2,i3);1]];
        v=[[1;1;1;1],...
           [i1;1;1;1],...
           [i1;i2;1;1],...
           [i1;i2;i3;1]]; 
        V.mat=x*pinv(v); 
        spm_write_vol(V,Data); 
        delete(out.files{i}); 
        delete(J.images{i}{1}); 
    end; 
end;


function [y1,y2,y3] = spmj_affine_transform(x1,x2,x3,M)
% function [y1,y2,y3] = affine_transform(x1,x2,x3,M)
% -----------------------------------------------------------------
% Affine Transform for input stuff in any format (N-dim strcutures)
% -----------------------------------------------------------------
y1 = M(1,1)*x1 + M(1,2)*x2 + M(1,3)*x3 + M(1,4);
y2 = M(2,1)*x1 + M(2,2)*x2 + M(2,3)*x3 + M(2,4);
y3 = M(3,1)*x1 + M(3,2)*x2 + M(3,3)*x3 + M(3,4);


