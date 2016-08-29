function [DATA,VT]=spmj_crop_vol(source,target,x,y,z);
% spmj_crop_vol(source,target,x,y,z);
% spmj_crop_volume makes a cropped Image from the source image
% voxels in range x,y,z are taken over from the source image
% resolution is preserved
% x,y,z are in voxel space 
% v 2.0: Joern Diedrichsen j.diedrichsen@ucl.ac.uk
% v 2.1: Supports mutivol files
% v 2.4: Compatible with SPM8
if (isstruct(source))
    VS=source;
else
    VS=spm_vol(source);
end;
xdim=length(x);
ydim=length(y);
zdim=length(z);
spmVer=spm('Ver');
N=length(VS);
[X,Y]=meshgrid(x,y);
DATA=zeros([xdim ydim]);

for n=1:N
    switch (spmVer)
        case 'SPM2'
            VT(n,1)    = struct(		'fname',	target,...
                'dim',		[xdim ydim zdim VS(n).dim(4)],...
                'mat',		VS(n).mat*spm_matrix([x(1)-1 y(1)-1 z(1)-1]),...
                'pinfo',	VS(n).pinfo,...
                'descrip',	'cropped');
        otherwise
            VT(n,1)    = struct(		'fname',	target,...
                'dim',		[xdim ydim zdim],...
                'dt',       VS(n).dt,...
                'mat',		VS(n).mat*spm_matrix([x(1)-1 y(1)-1 z(1)-1]),...
                'pinfo',	VS(n).pinfo,...
                'n',        VS(n).n,...
                'descrip',	'cropped');
    end;
end;

VT    = spm_create_vol(VT);
for n=1:N
    for i=1:zdim
        Z=ones(ydim,xdim)*z(i);
        DATA=spm_sample_vol(VS(n),X,Y,Z,1)';
        VT(n)=spm_write_plane(VT(n),DATA,i);
    end;
end;

if (strcmp(spmVer,'SPM2'))
    spm_close_vol(VT);
end;
