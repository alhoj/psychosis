clear

addpath(genpath('/m/nbe/scratch/psykoosi/scripts'));

cfg=[];

% Define subjects of interest here
cfg.subs = {
'EPHE602'
'EPHE606'
'EPHE607'
};

cfg.cond='BL'; % use follow-up data for subjects of interest; options baseline 'BL' or follow-up 'FU'
cfg.res='2mm'; % voxel resolution; options '2mm', '4mm', '8mm', '16mm', or '32mm'
cfg.mask=[]; % if empty, use the default MNI152 mask, e.g. "MNI152_T1_4mm_brain_mask.nii"
cfg.atlas='brainnetome_atlas_w_cerebellum_v2.nii'; % define atlas; leave empty for voxel-to-voxel connectivity
cfg.roiTC='PCA'; % how to extract the roi timecourse; options 'PCA' or 'mean'
cfg.outdir='test'; % label of the output folder
cfg.useNonSpatialSmoothedData=0;
%%
cfg=calculate_FC(cfg);
