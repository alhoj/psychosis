clear
% delete logs/*
delete jobs/*

%% INPUT DATA // An array of strings of valid files with location

cfg=[];  % this will contain the parameters of our expriment

path='/m/nbe/scratch/psykoosi/';
subs={'consN29','patsN36'};
conds={'BL','FU'};
data_labels={
'1PCsFromConsN15regressedOut'
'5PCsFromConsN15regressedOut'
'10PCsFromConsN15regressedOut'
};
methods={'covmat'};
norms={'detrend'};
res='4mm';

cfg.func='funpsy_slurm_script';

% slurm parameters
cfg.partition='batch';
cfg.mem='32000';
cfg.time='00:59:59';

cfg.ind=100;
for s=1:length(subs)
    for c=1:length(conds)
        for d=1:length(data_labels)
            for m=1:length(methods)
                for n=1:length(norms)
                    cfg.indir=[path 'data/' conds{c} '_' data_labels{d} '_' methods{m} '_36pt29hc_' res '_' norms{n} '_butterBandpass/'];
                    cfg.outdir = [path 'SBPS/' conds{c} '_' subs{s} '_' data_labels{d} '_' methods{m} '_' norms{n} '_voxel2voxel4mm'];
                    if ~exist(cfg.outdir,'dir') 
                        system(['mkdir -p ' cfg.outdir]);
                    end   
                    % NAME OF YOUR ANALYSIS SESSION
                    cfg.session_name = ['SBPS_' conds{c} '_' subs{s} '_' data_labels{d} '_' methods{m} '_' norms{n} '_voxel2voxel4mm'];
                    cfg.subs=subs{s};
                    cfg.ind=cfg.ind+1; % this is just the job (and log) index
                    function_make_scripts_slurm(cfg)
                end
            end
        end
    end
end

%% Run the jobs

system('source slurm_run_jobs_auto.sh');

