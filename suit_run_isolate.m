function suit_run_isolate(job)
% function suit_run_isolate(job)
% Runs an Suit-toolbox Isolation job from the job manager
% 
% v.3.0: Compatibility with SPM job manager (26/08/2010) 
%_______________________________________________________________________
% Copyright (C) 2010 
% Joern Diedrichsen (j.diedrichsen@ucl.ac.uk)
s=job.source;
for i=1:length(s)
     suit_isolate(s{i}{1},'bb',job.bb,'cerebral_range',job.cerebral_range,'cerebellar_range',job.cerebellar_range); 
end; 
