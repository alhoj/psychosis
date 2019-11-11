clear
delete logs/*
% delete jobs/*

addpath(genpath('/m/nbe/scratch/psykoosi/scripts'))

% Define subjects of interest here; PCs regressed out from this group
subs = {
'EPHE602'
'EPHE606'
'EPHE607'
'EPHE608'
'EPHE610'
'EPHE613'
'EPHE617'
'EPHE618'
'EPHE619'
'EPJO208'
'EPJO210'
'EPJO211'
'EPJO213'
'EPJO214'
'EPJO215'
'EPJO216'
'EPJO217'
'EPJO218'
'EPJO220'
'EPJO222'
'EPJO224'
'EPJO225'
'EPJO227'
'EPJO228'
'EPKK303'
'EPPE114'
'EPPE115'
'EPPE117'
'EPPE118'
'EPPE119'
'EPPE120'
'EPPE121'
'EPPE122'
'EPPE124'
'EPPE126'
'EPPE128'
'EPPE129'
'EPPE131'
'EPPE132'
'EPPE134'
'EPPE135'
'EPPE136'
'EPPE137'
'EPPE139'
'EPPE140'
'EPPE143'
'EPPE145'
'EPPE148'
'EPPE149'
'EPPE150'
'EPPE151'
'EPPE152'
'EPPE153'
'EPPE154'
'EPPE157'
'EPPE158'
'EPPK401'
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
conds={'BL'};
NPCs={1, 5, 10};

cfg=[];
% cfg.subs='BothN66'; % subjects of interest
cfg.refs='ConsN15'; % reference group; PCs calculated from this group
% cfg.condSubs='BL'; % which data to use for subjects of interest; options baseline 'BL' or follow-up 'FU'
cfg.condRefs='BL'; % which data to use for reference subjects
cfg.res='2mm'; % voxel resolution; options '2mm', '4mm', '8mm', '16mm', or '32mm'
cfg.mask=[]; % if empty, use the default MNI152 mask, e.g. "MNI152_T1_2mm_brain_mask.nii"
% cfg.NPC=1; % how many PCs to be regressed out
cfg.method='covmat'; % how to calculate the components; options: using covariance matrices as in maxcorr-function ('covmat') or directly from voxelwise BOLD-signals ('bold') 
cfg.normalization='zscore'; % 'zscore', 'detrend' or 'none'
cfg.comps='extrinsic'; % which components to regress out; options 'extrinsic' or 'intrinsic'
cfg.nonparametric=1; % 1=find common component count using permutations for circularly shifted timeseries; 0=find common component count using parametric estimate (conservative); applies only if cfg.comps is 'intrinsic'
% cfg.indir=['/m/nbe/scratch/psykoosi/data/' cfg.condSubs '_36pt29hc_2mm_butterBandpass/'];
% cfg.outdir=['/m/nbe/scratch/psykoosi/data/' cfg.condSubs '_' num2str(cfg.NPC) 'PCsFromConsN15regressedOut_' cfg.method '_36pt29hc_2mm_' cfg.normalization '_butterBandpass']; % output folder

% slurm parameters
cfg.partition='batch';
cfg.mem='32000';
cfg.time='00:09:59';

cfg.func='regress_PCs_new(cfg)';
%% make jobs
counter=1001;
for c=1:length(conds)
    for p=1:length(NPCs)
        for iter=1:length(subs)
            counter=counter+1;
            cfg.condSubs=conds{c};
            cfg.NPC=NPCs{p};
            cfg.indir=['/m/nbe/scratch/psykoosi/data/' cfg.condSubs '_' cfg.res '_butterBandpass/'];
%             cfg.outdir=['/m/nbe/scratch/psykoosi/data/' cfg.condSubs '_' num2str(cfg.NPC) 'PCsFromConsN15regressedOut_' cfg.method '_57pt29hc_' cfg.res '_' cfg.normalization '_butterBandpass'];
            cfg.outdir=['/m/nbe/scratch/psykoosi/data/' cfg.condSubs '_' num2str(cfg.NPC) 'PCsConsN15MaxCorr_nonparamCompCount_' cfg.method '_57pt29hc_' cfg.res '_' cfg.normalization '_butterBandpass'];
            cfg.ind=counter; % this is just the job (and log) index
            cfg.infile=subs{iter};
            cfg.outfile=cfg.infile;
            function_make_scripts_slurm(cfg)
        end
    end
end

%% Run the jobs

system('source slurm_run_jobs_auto.sh');

