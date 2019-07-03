clear
% delete logs/*
delete jobs/*

addpath(genpath('/m/nbe/scratch/psykoosi/scripts'))

cfg=[];
cfg.subs='BothN66'; % subjects of interest
cfg.refs='ConsN15'; % reference group
cfg.condSubs='FU'; % which data to use for subjects of interest; options baseline 'BL' or follow-up 'FU'
cfg.condRefs='BL'; % which data to use for reference subjects
cfg.res='2mm'; % voxel resolution; options '2mm', '4mm', '8mm', '16mm', or '32mm'
cfg.mask=''; % if empty, use the default MNI152 mask, e.g. "MNI152_T1_2mm_brain_mask.nii"
cfg.useMeanOverRefs=0; % take average over the reference group and use that to calculate ISC; options 1->yes or 0->no 
cfg.useNonSpatialSmoothedData=0;

% slurm parameters
cfg.partition='batch';
cfg.mem='64000';
cfg.time='00:59:59';

%% Make jobs -> one job per time window
notps=245; % time points in the original niftis
twl=30; % time window length
step=1; % step between consecutive time windows
counter=0;
for iter=1:step:notps
    counter=counter+1;
    tw=iter:(twl+iter-1);
    cfg.toi=tw; % time points of interest as vector
    if cfg.useNonSpatialSmoothedData
        cfg.outdir=[cfg.condSubs '_' num2str(twl) 'TRwin' num2str(step) 'TRstep_ref' cfg.refs '_MNI152wholeBrainMask_' cfg.res '_noSpatialSmoothing/timeInt' num2str(counter)]; % label of the output folder
    else
        cfg.outdir=[cfg.condSubs '_' num2str(twl) 'TRwin' num2str(step) 'TRstep_ref' cfg.refs '_MNI152wholeBrainMask_' cfg.res '/timeInt' num2str(counter)]; % label of the output folder
    end
    
    cfg.func='calculate_indepRef_twISC(cfg)';
    cfg.ind=1000+counter; % this is just the job (and log) index
    function_make_scripts_slurm(cfg)
    if tw(end)>=notps
        break
    end
end

%% Run the jobs

system('source slurm_run_jobs_auto.sh');
