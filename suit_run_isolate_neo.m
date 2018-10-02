function suit_run_isolate_neo(job)
% function suit_run_isolate_neo(job)
% Runs an Suit-toolbox Isolation job from the job manager
% 
% v.3.0: Compatibility with SPM job manager (26/08/2010) 
%_______________________________________________________________________

s=job.source;
for i=1:length(s)
     suit_isolate_neo(s{i}{1},'maskp',job.maskp,'keeptfiles',job.keeptfiles,'getV',job.getV); 
end; 
