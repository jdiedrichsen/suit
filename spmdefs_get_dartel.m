function [Def,mat]=spmdefs_get_dartel(flowfield,Affine)
% Gets a deformation flow field from Dartel 
% INPUT: flowfield: Name of the flowfield file 
%        Affine:    
% OUTPUT: 
%        Deformation 

% If no affine transformation given, set to identity 

if (nargin<2 || isempty(Affine))
    Affine=eye(4); 
elseif (ischar(Affine))
    load(Affine); 
end; 

times  = [1 0]; 
K      = 6; 

N      = nifti(flowfield);

y      = spm_dartel_integrate(N.dat,times,K);
Def    = cell(3,1);
M      = single(inv(Affine)*N.mat);
mat    = N.mat0;
Def{1} = y(:,:,:,1)*M(1,1) + y(:,:,:,2)*M(1,2) + y(:,:,:,3)*M(1,3) + M(1,4);
Def{2} = y(:,:,:,1)*M(2,1) + y(:,:,:,2)*M(2,2) + y(:,:,:,3)*M(2,3) + M(2,4);
Def{3} = y(:,:,:,1)*M(3,1) + y(:,:,:,2)*M(3,2) + y(:,:,:,3)*M(3,3) + M(3,4);
