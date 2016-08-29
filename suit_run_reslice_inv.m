function out = suit_run_reslice_inv(varargin)
% function out = suit_run_reslice_inv(job)
% Runs an Suit-toolbox reslice_invjob from the job manager
% 
% v.2.5: Compatibility with SPM job manager (26/08/2010) 
%_______________________________________________________________________
% Copyright (C) 2010 
% Joern Diedrichsen (j.diedrichsen@ucl.ac.uk)
job    = varargin{1};
suit_reslice_inv(job.resample,...
    job.paramfile,...
    'prefix',job.prefix,...
    'reference',job.reference,...
    'interp',job.interp);
