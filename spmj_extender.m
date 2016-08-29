function outfilename=spmj_extender(Source, width)
% function outfilename=spmj_extender(Source, width)
% INPUT: 
%   Source: Name of source file (img or nii) 
%   Width: thickness of extension in voxels 
% OUTPUT: 
%   outfilename: name under which the new image is stored 
% copyright by Niall Lally, 2009
% V 1.2 comments, value of output voxels==1, fix bug in the extension (Tobias Wiestler)
% 
% -----------------------------------------------------------------
% Get Source image
% -----------------------------------------------------------------
if (nargin<1 || isempty(Source))
        Source=spm_get(1,'*.img','Get Source Image');
end;
[source_dir,Source,ext,num]=spm_fileparts(Source);

if (isempty(source_dir))
    source_dir=pwd;
end;

Sourcec=fullfile(source_dir,Source);
V_source=spm_vol([Sourcec ext]);

%----get informativ voxels (value>0)
X=spm_read_vols(V_source);
i=find(X>0);
ind2sub(size(X),i);
 
%----define the kernelsize:  
%----width=2 => kernelsize=5
kernelsize=width*2+1;  

%----init kernel
x=[]; 
x=ones(1,kernelsize);

%----init the new volume
Y=zeros(size(X));  

%----X->source image; Y->output image; x->fx fy fz the separable form of; 
%----offsets  [i j k] contains the x, y and z shifts to reposition the output the function in x
spm_conv_vol(X,Y,x,x,x,[-width -width -width]);
Y=(Y>=1);

%----save the extent volume withe the prefix 'e'
V_source.fname=fullfile(source_dir,['e' Source ext]);
spm_write_vol(V_source,Y);

%----give back the name of the new file 
outfilename=V_source.fname;