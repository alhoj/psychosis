clear

addpath(genpath('/m/nbe/scratch/psykoosi/scripts'));

cfg=[];

% Define 4D nifti files
cfg.subs = {
'EPVE516'
'EPVE517'
'EPVE520'
'EPVE524'
'EPVE525'
'EPVE526'
'EPVE527'
'EPVE530'
'EPVE531'
'EPVE533'
'EPVE536'
'EPVE540'
'EPVE541'
'EPVE542'
'EPVE544'
'EPVE546'
'EPVE548'
'EPVE549'
'EPVE551'
'EPVE554'
'EPVE555'
'EPVE556'
'EPVE557'
'EPVE558'
'EPVE559'
'EPVE560'
'EPVE561'
'EPVE562'
'EPVE563'
};
%%
cfg.indir='/m/nbe/scratch/psykoosi/data/1PCsFromConsN15regressedOut_36pt29hc_2mm_zscored/';
cfg.mask='/m/nbe/scratch/psykoosi/masks/MNI152_T1_2mm_brain_mask.nii';
cfg.atlas='/m/nbe/scratch/psykoosi/masks/brainnetome_atlas_w_cerebellum_v2.nii'; % define atlas; leave empty for voxel-to-voxel connectivity
cfg.roiTC='PCA'; % how to extract the roi timecourse; options 'PCA' or 'mean'
cfg.outdir='consN29_1PCsFromConsN15regressedOut'; % label of the output folder

%%
cfg=calculate_FC(cfg);
