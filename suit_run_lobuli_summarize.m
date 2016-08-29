function out = suit_run_lobuli_summarize(varargin)
% function out = suit_run_lobuli_summarize(job)
% Runs an Suit-toolbox summarize job from the job manager
% 
% v.2.5: Compatibility with SPM job manager (26/08/2010) 
%_______________________________________________________________________
% Copyright (C) 2010 
% Joern Diedrichsen (j.diedrichsen@ucl.ac.uk)

job    = varargin{1};
suit_lobuli_summarize(job.images,...
    'atlas',job.atlas{1},...
    'outfilename',fullfile(job.outdir{1},job.output),...
    'stats',{job.stats});
