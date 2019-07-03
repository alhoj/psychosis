clear
% delete logs/*
% delete jobs/*

cfg=[];
cfg.subs='PatsN36'; % subjects of interest: 'PatsN36', 'ConsN30', or 'BothN66'
cfg.refs='ConsN15'; % reference group: 'ConsN15', 'PatsN36', 'ConsN30', or 'BothN66'
cfg.condSubs='BL'; % which data to use for subjects of interest; options baseline 'BL' or follow-up 'FU'
cfg.condRefs='BL'; % which data to use for reference subjects
cfg.res='2mm'; % voxel resolution; options '2mm', '4mm', '8mm', '16mm', or '32mm'
cfg.mask=''; % if empty, use the default MNI152 mask, e.g. "MNI152_T1_2mm_brain_mask.nii"
cfg.useMeanOverRefs=0; % take average over the reference group and use that to calculate ISPS; options 1->yes or 0->no 
cfg.useNonSpatialSmoothedData=1; 
cfg.averageOverTimePoints=0;

%% Make jobs
if cfg.useNonSpatialSmoothedData
    cfg.outdir=[cfg.condSubs '_ref' cfg.refs '_MNI152wholeBrainMask_' cfg.res '_noSpatialSmoothing']; % label of the output folder
else
    cfg.outdir=[cfg.condSubs '_ref' cfg.refs '_MNI152wholeBrainMask_' cfg.res]; % label of the output folder
end
cfg.func='calculate_indepRef_ISPS(cfg)';
cfg.ind=1; % this is just the job (and log) index
cfg.toi=[]; % not applicable here

function_make_scripts_slurm(cfg)

%% Run the jobs

make_slurm_run_jobs
system('source slurm_run_jobs_auto.sh');