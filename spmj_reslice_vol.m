function spmj_reslice_vol(VS,dim,mat,name);
% spmj_reslice_vol(Source,dim,mat,name);
% Reslices source volume into a new voxel format
% Source: memory-mapped volume file
% dim: [x y z], size of the new image in 3d
% mat: 4*4 transformation matrix to indicate the location within Source
% name: Output file name
% version 2.0: can also reslice 4d-volumes
% j.diedrichsen@bangor.ac.uk
if ischar(VS)
    VS=spm_vol(VS);
end;
if (any(size(dim)~=[1 3]))
    error ('dimension has to be a 3x1 vector');
end;
if (any(size(mat)~=[4 4]))
    error ('transformation matrix has to be a 4x4 matrix');
end;

Data=zeros(dim);
spmVer=spm('Ver');

numimages=length(VS);
for i=1:numimages
    if (strcmp(spmVer,'SPM2'))
        VT    = struct(		'fname',	name,...
            'dim',		[dim VS.dim(4)],...
            'mat',		mat,...
            'pinfo',	VS.pinfo,...
            'descrip',	'resliced');
    else
        VT    = struct(		'fname',	name,...
            'dim',		dim,...
            'dt',       VS(i).dt,...
            'mat',		mat,...
            'pinfo',	VS(i).pinfo,...
            'n',        [i 1],...
            'descrip',	'resliced');
    end;


    VT = spm_create_vol(VT);
    A=inv(VS(i).mat)*mat;
    
    for j=1:dim(3)
        DATA=spm_slice_vol(VS(i),A*spm_matrix([0 0 j]),dim(1:2),1)';
        VT = spm_write_plane(VT,DATA',j);
    end;
    fprintf('.');
end;
if (strcmp(spmVer,'SPM2'))
    spm_close_vol(VT);
end;
fprintf('\n');
