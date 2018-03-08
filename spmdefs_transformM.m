function [xyz_new] = spmdefs_transformM(Def,mat,xyz,varargin)
% function [xyz_new] = spmdefs_transform(Def,mat,xyz,varargin)
% Transforms coordinates from VG space into VF space 
% INPUT: 
%   Def: Deformation cell array {3} 
%   mat: affine transform for Def 
%   x,y,z: source locations in world coordinates 
% VARAGRGIN
%   'interp',1: Default trialinear  

interp=1; 
vararginoptions(varargin,{'interp'}); 
[vx,vy,vz]=spmj_affine_transform(xyz(:,1),xyz(:,2),xyz(:,3),inv(mat)); 
xyz_new(:,1)=spm_sample_vol(Def{1},vx,vy,vz,interp); 
xyz_new(:,2)=spm_sample_vol(Def{2},vx,vy,vz,interp); 
xyz_new(:,3)=spm_sample_vol(Def{3},vx,vy,vz,interp); 
