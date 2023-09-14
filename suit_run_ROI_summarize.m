function out = suit_run_ROI_summarize(varargin)
% function out = suit_run_lobuli_summarize(job)
% Runs an Suit-toolbox summarize job from the job manager
% 
% v.3.6: ROI run summarize initial version 
%_______________________________________________________________________
% Copyright (C) 2023
% Joern Diedrichsen (joern.diedrichsen@googlemail.com)

job    = varargin{1};
suit_ROI_summarize(job.images,...
    'atlas',job.atlas{1},...
    'outfilename',fullfile(job.outdir{1},job.output),...
    'stats',{job.stats});
