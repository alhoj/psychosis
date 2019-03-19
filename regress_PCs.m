function cfg=regress_PCs(cfg)
% Regress out BOLD signal principal components of one group of subjects
% from another group of subjects
% 
% Usage:
%   cfg = regress_PCs(cfg);
%
%   Input:
%   cfg.subs = IDs (e.g. 'EPHE602') of subjects of interest
% 
%   cfg.refs = IDs (e.g. 'EPVE516') of reference group whose data are used 
%   to calculate the PCs 
% 
%   cfg.mask = nifti file of binary mask; default is the MNI152 whole brain
%   mask
% 
%   cfg.res = resolution of the nifti files;  can be '2mm', '4mm', '6mm',
%   '8mm', '16mm', or '32mm'
%
%   cfg.cond = data condition to use; options baseline 'BL' or follow-up 'FU'
%
%   cft.outdir = label of the output folder where the individual ISC niftis are saved

%% Input validation

if ~iscell(cfg.subs)
    error('cfg.subs should contain a cell array!')
end
if ~iscell(cfg.refs)
    error('cfg.refs should contain a cell array!')
end
if ~ismember(cfg.condSubs,{'BL','FU'})
    error('cfg.condSubs must be either ''BL'' or ''FU''!')
end
if ~ismember(cfg.condRefs,{'BL','FU'})
    error('cfg.condRefs must be either ''BL'' or ''FU''!')
end
if ~ismember(cfg.res,{'2mm', '4mm', '6mm', '8mm', '16mm', '32mm'}) 
    error('cfg.res has to be either ''2mm'', ''4mm'', ''6mm'', ''8mm'', ''16mm'', or ''32mm''')
end
if ~isempty(cfg.mask) && ~isfile(cfg.mask)
    error(['could not find mask: ' cfg.mask])
end

disp(['Regress out PCs calculated from ' num2str(length(cfg.refs)) ' references from ' num2str(length(cfg.subs)) ' subjects of interested'])
disp(['using ' cfg.condRefs ' data for references and ' cfg.condSubs ' data for subjects of interested'])
fprintf('\n')

%% Load mask

disp('Loading mask...')
if isempty(cfg.mask)
    mask=load_nii(['/m/nbe/scratch/psykoosi/masks/MNI152_T1_' num2str(cfg.res) '_brain_mask.nii']);
else
    mask=load_nii(cfg.mask);
end
inmask=find(mask.img);
fprintf('\n')

%% Load brain data
% Load an example NIFTI file to find the dimensions needed for
% preallocating the matrices (this will speed things up)
test_nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.condSubs '/' cfg.subs{1} '/epi_preprocessed_' cfg.res '.nii']);

% Put time first, so that we can acess the x y z with the 1-D indices
temp=permute(test_nii.img,[4 1 2 3]);

% Get the number of time points, voxels, subjects-of-interest, and
% reference subjects
ntps=size(temp,1);
nvox=length(inmask);
nsubs=length(cfg.subs);
nrefs=length(cfg.refs);

% Preallocate data matrices with the dimensions of time points x subjects x voxels 
allsubs=zeros(ntps,nsubs,nvox); 
allrefs=zeros(ntps,nrefs,nvox);

disp('Loading brain data of subjects of interest...')
for i=1:nsubs
    disp(['Subject ' cfg.subs{i} ' - ' num2str(i) ' out of ' num2str(nsubs)])
    if cfg.useNonSpatialSmoothedData
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.condSubs '/' cfg.subs{i} '/epi_noSpatialSmoothing_' cfg.res '.nii']);
    else
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.condSubs '/' cfg.subs{i} '/epi_preprocessed_' cfg.res '.nii']);
    end
    temp=permute(nii.img,[4 1 2 3]);
    if cfg.zscore
        allsubs(:,i,:)=zscore(temp(:,inmask));
    else
        allsubs(:,i,:)=temp(:,inmask);
    end
end
fprintf('\n')
disp('Loading brain data of reference subjects...')
for i=1:nrefs
    disp(['Subject ' cfg.refs{i} ' - ' num2str(i) ' out of ' num2str(nrefs)])
    if cfg.useNonSpatialSmoothedData
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.condRefs '/' cfg.refs{i} '/epi_noSpatialSmoothing_' cfg.res '.nii']);
    else
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.condRefs '/' cfg.refs{i} '/epi_preprocessed_' cfg.res '.nii']);
    end
    temp=permute(nii.img,[4 1 2 3]);
    if cfg.zscore
        allrefs(:,i,:)=zscore(temp(:,inmask));
    else
        allrefs(:,i,:)=temp(:,inmask);
    end
end

fprintf('\n')

%% Regress out PCs of reference group from the subjects of interest
disp(['Calculating ' num2str(cfg.NPC) ' PCs from references and regressing them out from subjects of interest...'])
dataReg=zeros(nvox,ntps,nsubs);
for voxi=1:nvox
    if mod(voxi,1000)==0 % Show status every 1000 voxels
        disp([num2str(voxi) '/' num2str(nvox) ' voxels'])
    end
    [~,PCs] = pca(allrefs(:,:,voxi)); % Get the PCs
    PCs=PCs(:,1:cfg.NPC); % Take cfg.NPC number of PCs

    % Regress out PCs
    model=[ones(ntps,1) PCs];
    b=allsubs(:,:,voxi)'/model';
    TCreg=allsubs(:,:,voxi)-model*b';
    dataReg(voxi,:,:)=TCreg;
end
fprintf('\n')

%% Save files
% Replace possible NaN values with zeros.
dataReg(isnan(dataReg))=0;

% Check if output directory exists; if not, create it
dirname=['/m/nbe/scratch/psykoosi/data/' cfg.outdir];
if ~exist(dirname,'dir') 
    system(['mkdir -p ' dirname]);
end

% Save cleaned niftis
disp(['Saving files to directory ' dirname '/'])
dataReg=permute(dataReg,[2 1 3]);
for i=1:nsubs
    newbrain=zeros(ntps,size(mask.img,1),size(mask.img,2),size(mask.img,3));
    newbrain(:,inmask)=dataReg(:,:,i);
    newbrain=permute(newbrain,[2 3 4 1]);
    filename=[dirname '/' cfg.subs{i} '.nii'];
    save_nii(make_nii(newbrain),filename);
    nii=fixOriginator(filename,mask);
    save_nii(nii,filename);
end

disp('Done!')
fprintf('\n')