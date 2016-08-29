function suit_run_flatmap(varargin)
% function out = suit_run_reslice(job)
% Runs an Suit-toolbox reslice job from the job manager
% 
% v.3.0: First version of the cerebellar flatmap 
%_______________________________________________________________________
% Copyright (C) 2015 
% Joern Diedrichsen (joern.diedrichsen@googlemail.com)
flatSpace={'SUIT','FSL','SPM'}; 
job    = varargin{1};

Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph);

Data=suit_map2surf(job.images{1},...
        'stats',job.stats,...
        'space',flatSpace{job.flatSpace}); 
    
if (job.flatBorder) 
    border = 'fissures_flat.mat';
else 
    border = []; 
end; 

suit_plotflatmap(Data,'cscale',job.flatCscale,'threshold',...
    job.flatThreshold,'border',border,'cmap',eval(job.Cmap)); 


