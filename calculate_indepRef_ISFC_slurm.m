clear
delete logs/*
% delete jobs/*

subs={'consN29','patsN36'};
conds={'BL','FU'};

cfg=[];
% cfg.subs='consN29'; % subjects of interest
cfg.refs='consN15'; % reference group
% cfg.cond='BL'; % which data to use for subjects of interest; options baseline 'BL' or follow-up 'FU'
cfg.res='2mm'; % voxel resolution; options '2mm', '4mm', '8mm', '16mm', or '32mm'
cfg.mask='/m/nbe/scratch/psykoosi/masks/EPIgroupMask_2mm.nii'; % if empty, use the default MNI152 mask, e.g. "MNI152_T1_2mm_brain_mask.nii"
cfg.atlas='/m/nbe/scratch/psykoosi/masks/brainnetome_atlas_w_cerebellum_v2.nii'; % define atlas; leave empty for voxel-to-voxel connectivity
% cfg.atlas=[];
cfg.roiTC='pca'; % how to extract roi timecourse; options 'pca' or 'mean'
cfg.symmetrize=0; % 0 or 1; symmetrize the correlation matrix; i.e. average between sub1-roi1 vs ref1-roi2 and sub1-roi2 vs ref1-roi1

% slurm parameters
cfg.partition='batch';
cfg.mem='120000';
cfg.time='00:29:59';

%% Make jobs

cfg.ind=5;
for s=1:length(subs)
    for c=1:length(conds)
        cfg.subs=subs{s};
        cfg.cond=conds{c};
        cfg.outdir=[cfg.cond '_' cfg.subs '_brainnetome_EPIgroupMask_patsN28ref_unsymmetrized']; % label of the output folder
        cfg.func='calculate_indepRef_ISFC(cfg)';
        function_make_scripts_slurm(cfg)
        cfg.ind=cfg.ind+1; % this is just the job (and log) index
    end
end


%% Run the jobs

system('source slurm_run_jobs_auto.sh');
