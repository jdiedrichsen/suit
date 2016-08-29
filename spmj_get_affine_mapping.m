function M = spmj_get_affine_mapping(VF,VG,aflags)
% M = spmj_get_affine_mapping(VF,VG,aflags)
% VF: Source 
% VG: Target 
% get affine mapping to template 
% Stolen from spm_segment 
% used to get priors into the space of a anatomical image 
% The returned matrix is the affine transformation matrix of the space of
% the anatomical image to the space of the template image 
% M_atlas: tranform of Atlas
% M_source: transform of source 
% Total transformation from atlas VOXEL to source VOXEL: 
% inv(M*M_source)*M_atlas or M_atlas\(M*M_source)
if ~isempty(VG) & ischar(VG), VG = spm_vol(VG); end;

if ~isempty(VG) & isstruct(VG),
	% Affine registration so that a priori images match the image to
	% be segmented.
	%-----------------------------------------------------------------------

	VFS = spm_smoothto8bit(VF(1),aflags.smosrc);

	% Scale all images approximately equally
	% ---------------------------------------------------------------
	for i=1:length(VG),
		VG(i).pinfo(1:2,:) = VG(i).pinfo(1:2,:)/spm_global(VG(i));
	end;
	VFS(1).pinfo(1:2,:) = VFS(1).pinfo(1:2,:)/spm_global(VFS(1));

	spm_plot_convergence('Init','Affine Registration','Mean squared difference','Iteration');
	flags     = struct('sep',aflags.smosrc, 'regtype',aflags.regtype,'WG',[],'globnorm',0,'debug',0,'WF',aflags.WF);
	M         = eye(4);
	[M,scal]  = spm_affreg(VG, VFS, flags, M);

	if ~isempty(aflags.weight), flags.WG = spm_vol(aflags.weight); end;
    
    flags.sep = aflags.smosrc/2;
	M         = spm_affreg(VG, VFS, flags, M,scal);
	spm_plot_convergence('Clear');

elseif all(size(VG) == [4 4])
	% Assume that second argument is a matrix that will do the job
	%-----------------------------------------------------------------------
	M = VG;
else
	% Assume that image is normalized
	%-----------------------------------------------------------------------
	M = eye(4);
end
return;
%=======================================================================
