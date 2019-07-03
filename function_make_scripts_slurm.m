function function_make_scripts_slurm(cfg)


% delete logs/*
% delete jobs/*

% mkdir jobs
% mkdir logs

if ~isfield(cfg,'averageOverTimePoints')
    cfg.averageOverTimePoints=0;
end
if ~isfield(cfg,'useNonSpatialSmoothedData')
    cfg.useNonSpatialSmoothedData=0;
end
if ~isfield(cfg,'useMeanOverRefs')
    cfg.useMeanOverRefs=0;
end
if ~isfield(cfg,'partition')
    cfg.partition='batch';
end
if ~isfield(cfg,'mem')
    cfg.mem='100000';
end
if ~isfield(cfg,'time')
    cfg.time='05:59:59';
end
if ~isfield(cfg,'subs')
    cfg.subs=[];
end
if ~isfield(cfg,'refs')
    cfg.refs=[];
end
if ~isfield(cfg,'cond')
    cfg.cond=[];
end
if ~isfield(cfg,'condSubs')
    cfg.condSubs=[];
end
if ~isfield(cfg,'condRefs')
    cfg.condRefs=[];
end
if ~isfield(cfg,'group1')
    cfg.group1=[];
end
if ~isfield(cfg,'group2')
    cfg.group2=[];
end
if ~isfield(cfg,'twID')
    cfg.twID=[];
end
if ~isfield(cfg,'count')
    cfg.count=0;
end
if ~isfield(cfg,'mask')
    cfg.mask=[];
end
if ~isfield(cfg,'toi')
    cfg.toi=[];
end
if ~isfield(cfg,'res')
    cfg.res='2mm';
end
if ~isfield(cfg,'regType')
    cfg.regType=[];
end
if ~isfield(cfg,'regi')
    cfg.regi=1;
end
if ~isfield(cfg,'ti')
    cfg.ti=1;
end
if ~isfield(cfg,'infile')
    cfg.infile=[];
end
if ~isfield(cfg,'clusterstatistic')
    cfg.clusterstatistic='maxsum';
end
if ~isfield(cfg,'alpha')
    cfg.alpha=0.05;
end
if ~isfield(cfg,'cdtP')
    cfg.cdtP=0.01;
end
if ~isfield(cfg,'NPC')
    cfg.NPC=1;
end
if ~isfield(cfg,'normalization')
    cfg.normalization='none';
end
if ~isfield(cfg,'method')
    cfg.method=[];
end
if ~isfield(cfg,'comps')
    cfg.comps=[];
end
if ~isfield(cfg,'nonparametric')
    cfg.nonparametric=0;
end
if ~isfield(cfg,'indir')
    cfg.indir=[];
end
if ~isfield(cfg,'infile')
    cfg.infile=[];
end
if ~isfield(cfg,'indata')
    cfg.indata=[];
end
if ~isfield(cfg,'outdir')
    cfg.outdir=[];
end
if ~isfield(cfg,'outfile')
    cfg.outfile=cfg.infile;
end
if ~isfield(cfg,'atlas')
    cfg.atlas=[];
end
if ~isfield(cfg,'roiTC')
    cfg.roiTC='mean';
end
if ~isfield(cfg,'session_name')
    cfg.session_name=[];
end
if ~isfield(cfg,'symmetrize')
    cfg.symmetrize=0;
end


jobind=cfg.ind;
disp('Making jobs ...')

fid = fopen(['./jobs/job_' num2str(jobind) '.sh'],'w');

fprintf(fid,'#!/bin/bash\n\n');

fprintf(fid,['#SBATCH -p ' cfg.partition '\n']);
fprintf(fid,['#SBATCH -t ' cfg.time '\n']);
fprintf(fid,['#SBATCH -o ' './logs/log_' num2str(jobind) '\n']);
fprintf(fid,'#SBATCH --qos=normal\n');
% fprintf(fid,'#SBATCH --exclusive\n');
fprintf(fid,['#SBATCH --mem-per-cpu=' cfg.mem '\n\n']);

% fprintf(fid,'sleep 5s\n\n');

fprintf(fid,['matlab -nojvm -r "cd ' pwd '/; cfg.group1=''' cfg.group1 '''; cfg.group2=''' cfg.group2 '''; cfg.subs=''' cfg.subs '''; cfg.refs=''' cfg.refs '''; cfg.cond=''' cfg.cond '''; cfg.condSubs=''' cfg.condSubs '''; cfg.condRefs=''' cfg.condRefs '''; cfg.regType=''' cfg.regType '''; cfg.method=''' cfg.method '''; cfg.comps=''' cfg.comps '''; cfg.nonparametric=' num2str(cfg.nonparametric) '; cfg.symmetrize=' num2str(cfg.symmetrize) '; cfg.roiTC=''' cfg.roiTC '''; cfg.res=''' cfg.res '''; cfg.mask=''' cfg.mask '''; cfg.atlas=''' cfg.atlas '''; cfg.infile=''' cfg.infile '''; cfg.indata=''' cfg.indata '''; cfg.outfile=''' cfg.outfile '''; cfg.session_name=''' cfg.session_name '''; cfg.normalization=''' cfg.normalization '''; cfg.twID=''' cfg.twID '''; cfg.clusterstatistic=''' cfg.clusterstatistic '''; cfg.alpha=' num2str(cfg.alpha) '; cfg.cdtP=' num2str(cfg.cdtP) '; cfg.regi=' num2str(cfg.regi) '; cfg.NPC=' num2str(cfg.NPC) '; cfg.ti=' num2str(cfg.ti) '; cfg.useMeanOverRefs=' num2str(cfg.useMeanOverRefs) '; cfg.useNonSpatialSmoothedData=' num2str(cfg.useNonSpatialSmoothedData) '; cfg.averageOverTimePoints=' num2str(cfg.averageOverTimePoints) '; cfg.toi=[' num2str(cfg.toi) ']; cfg.indir=''' cfg.indir '''; cfg.outdir=''' cfg.outdir '''; ' cfg.func ' ; exit;"']);


fclose(fid);

% disp('Done making jobs!')

end
