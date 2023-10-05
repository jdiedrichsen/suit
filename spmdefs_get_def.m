function [Def,mat] = spmdefs_get_def(fname)
% Load a deformation field saved as an image
P      = [repmat(fname,3,1), [',1,1';',1,2';',1,3']];
V      = spm_vol(P);
Def    = cell(3,1);
Def{1} = spm_load_float(V(1));
Def{2} = spm_load_float(V(2));
Def{3} = spm_load_float(V(3));
mat    = V(1).mat;
