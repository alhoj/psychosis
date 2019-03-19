clear

addpath(genpath('/m/nbe/scratch/psykoosi/scripts'));

cfg=[];

% Define subjects of interest here
cfg.subs = {
'EPHE602'
'EPHE606'
'EPHE607'
};

% Define reference group here
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

cfg.condSubs='FU'; % use follow-up data for subjects of interest; options baseline 'BL' or follow-up 'FU'
cfg.condRefs='BL'; % use baseline data for reference subjects
cfg.res='8mm'; % voxel resolution; options '2mm', '4mm', '8mm', '16mm', or '32mm'
cfg.mask=[]; % if empty, use the default MNI152 mask, e.g. "MNI152_T1_4mm_brain_mask.nii"
cfg.useMeanOverRefs=0; % take average over the reference group and use that to calculate ISC; options 1->yes or 0->no 
cfg.outdir='test'; % label of the output folder
cfg.useNonSpatialSmoothedData=0;
%%
cfg=calculate_indepRef_ISC(cfg);
