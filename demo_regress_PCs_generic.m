clear

addpath(genpath('/m/nbe/scratch/psykoosi/scripts'));

cfg=[];

% Define reference group here; PCs calculated from this group
cfg.refs = {
'/m/nbe/scratch/psykoosi/data/BL/EPVE518/epi_preprocessed_8mm.nii'
'/m/nbe/scratch/psykoosi/data/BL/EPVE519/epi_preprocessed_8mm.nii'
'/m/nbe/scratch/psykoosi/data/BL/EPVE521/epi_preprocessed_8mm.nii'
'/m/nbe/scratch/psykoosi/data/BL/EPVE522/epi_preprocessed_8mm.nii'
'/m/nbe/scratch/psykoosi/data/BL/EPVE523/epi_preprocessed_8mm.nii'
'/m/nbe/scratch/psykoosi/data/BL/EPVE529/epi_preprocessed_8mm.nii'
'/m/nbe/scratch/psykoosi/data/BL/EPVE532/epi_preprocessed_8mm.nii'
'/m/nbe/scratch/psykoosi/data/BL/EPVE534/epi_preprocessed_8mm.nii'
'/m/nbe/scratch/psykoosi/data/BL/EPVE537/epi_preprocessed_8mm.nii'
'/m/nbe/scratch/psykoosi/data/BL/EPVE538/epi_preprocessed_8mm.nii'
};

cfg.mask='/m/nbe/scratch/psykoosi/masks/MNI152_T1_8mm_brain_mask.nii'; % if empty, use the default MNI152 mask, e.g. "MNI152_T1_2mm_brain_mask.nii"
cfg.NPC=10; % how many PCs to be regressed out
cfg.outID='controls5PCsRegressed'; % ID for the output niftis -> ..."_cfg.outID".nii

%%

% subjects of interest; PCs regressed out from this group
infiles = {
'/m/nbe/scratch/psykoosi/data/BL/EPHE602/epi_preprocessed_8mm.nii'
'/m/nbe/scratch/psykoosi/data/BL/EPHE606/epi_preprocessed_8mm.nii'
'/m/nbe/scratch/psykoosi/data/BL/EPHE607/epi_preprocessed_8mm.nii'
'/m/nbe/scratch/psykoosi/data/BL/EPHE608/epi_preprocessed_8mm.nii'
'/m/nbe/scratch/psykoosi/data/BL/EPHE610/epi_preprocessed_8mm.nii'
};

outfiles = {
'/m/nbe/scratch/psykoosi/data/PCsRegressed/EPHE602_8mm.nii'
'/m/nbe/scratch/psykoosi/data/PCsRegressed/EPHE606_8mm.nii'
'/m/nbe/scratch/psykoosi/data/PCsRegressed/EPHE607_8mm.nii'
'/m/nbe/scratch/psykoosi/data/PCsRegressed/EPHE608_8mm.nii'
'/m/nbe/scratch/psykoosi/data/PCsRegressed/EPHE610_8mm.nii'
};


for s=1:length(infiles)
    cfg.sub=infiles{s};
    cfg=regress_PCs(cfg);    
    save_nii(make_nii(cfg.dataReg),outfiles{s});
    temp=fixOriginator(outfiles{s},cfg.mask);
    save_nii(temp,outfiles{s});
end
