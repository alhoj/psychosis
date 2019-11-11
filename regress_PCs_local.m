clear

addpath(genpath('/m/nbe/scratch/psykoosi/scripts'));

cfg=[];

% Define subjects of interest here; PCs regressed out from this group
cfg.subs = {
'EPHE602'

};

% Define reference group here; PCs calculated from this group
cfg.refs = {
'EPVE518'
'EPVE519'
'EPVE521'
'EPVE522'
'EPVE523'
'EPVE529'
'EPVE532'
'EPVE534'
'EPVE537'
'EPVE538'
'EPVE539'
'EPVE545'
'EPVE547'
'EPVE550'
'EPVE552'
};
%%
cfg.condSubs='BL'; % baseline 'BL' or follow-up 'FU'
cfg.condRefs='BL'; % baseline 'BL' or follow-up 'FU'
cfg.res='2mm'; % voxel resolution; options '2mm', '4mm', '8mm', '16mm', or '32mm'
cfg.mask=[]; % if empty, use the default MNI152 mask, e.g. "MNI152_T1_2mm_brain_mask.nii"
cfg.NPC=1; % how many PCs to be regressed out
cfg.zscore=1; % zscore or not
cfg.useNonSpatialSmoothedData=0;
cfg.outdir=[cfg.condSubs '_' num2str(cfg.NPC) 'PCsFromConsN15regressedOut_36pt29hc_2mm_zscored']; % label of the output folder under /m/nbe/scratch/psykoosi/data/
%% run the script
regress_PCs(cfg);
