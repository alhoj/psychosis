function cfg=calculate_FC(cfg)
% Calculates voxel-to-voxel functional connectivity
% 
% Usage:
%   cfg = calculate_FC(cfg);
%
%   Input:
%   cfg.subs = IDs (e.g. 'EPHE602') of subjects
% 
%   cfg.cond = data condition to use; options baseline 'BL' or follow-up 'FU'
%
%   cfg.mask = nifti file of binary mask; default is the MNI152 whole brain
%   mask
% 
%   cfg.res = resolution of the nifti files;  can be '2mm', '4mm', '6mm',
%   '8mm', '16mm', or '32mm'
% 
%   cft.outdir = label of the output folder where the individual ISC nifitis are saved

%%
%   Output:
%   cfg.cors = correlation matrix (voxels x voxels x subjects)


%% Input validation

if ~iscell(cfg.subs)
    error('cfg.subs should contain a cell array!')
end
if ~ismember(cfg.cond,{'BL','FU'})
    error('cfg.condSubs must be either ''BL'' or ''FU''!')
end
if ~ismember(cfg.res,{'2mm', '4mm', '6mm', '8mm', '16mm', '32mm'}) 
    error('cfg.res has to be either ''2mm'', ''4mm'', ''6mm'', ''8mm'', ''16mm'', or ''32mm''')
end
if ~isempty(cfg.mask) && ~isfile(cfg.mask)
    error(['could not find mask: ' cfg.mask])
end

disp(['Calculating functional connectivity matrices for ' num2str(length(cfg.subs)) ' subjects based on Pearson''s correlation'])
disp(['using ' cfg.cond ' data'])
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
test_nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.cond '/' cfg.subs{1} '/epi_preprocessed_' cfg.res '.nii']);

% Put time first, so that we can acess the x y z with the 1-D indices
temp=permute(test_nii.img,[4 1 2 3]);

% Get the number of time points, voxels, subjects-of-interest, and
% reference subjects
ntps=size(temp(:,inmask),1);
nvox=size(temp(:,inmask),2);
nsub=length(cfg.subs);

% Preallocate data matrices with the dimensions of time points x subjects x voxels 
allsubs=zeros(ntps,nsub,nvox); 

disp('Loading brain data of subjects of interest...')
for i=1:nsub
    disp(['Subject ' cfg.subs{i} ' - ' num2str(i) ' out of ' num2str(nsub)])
    if cfg.useNonSpatialSmoothedData
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.cond '/' cfg.subs{i} '/epi_noSpatialSmoothing_' cfg.res '.nii']);
    else
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.cond '/' cfg.subs{i} '/epi_preprocessed_' cfg.res '.nii']);
    end
    temp=permute(nii.img,[4 1 2 3]);
    allsubs(:,i,:)=zscore(temp(:,inmask));
%     allsubs(:,i,:)=temp(:,inmask);
end
fprintf('\n')

%% Calculate FC
% Preallocate the correlation matrix
cors=zeros(nvox,nvox,nsub);
    
disp('Calculating FC...')
for subi=1:nsub
    disp(subi)
    % Calculate functional correlation over time between all voxels
    cors(:,:,subi)=corr(squeeze(allsubs(:,subi,:)));
  
end

% Replace possible NaN values with zeros.
cors(isnan(cors))=0;

fprintf('\n')

%% Average to obtain indivual ISC maps

% Take the average across subjects (first Fisher's
% z-tranform, then averaging, then back to correlation scale)
avg_cors=squeeze(tanh(mean(atanh(cors),3)));
% or without the Fisher's z-transform
%     avg_cors=squeeze(mean(cors,3));

% Replace possible NaN values with zeros.
avg_cors(isnan(avg_cors))=0;

% Check if output directory exists; if not, create it
dirname=['/m/nbe/scratch/psykoosi/FC/' cfg.outdir];
if ~exist(dirname,'dir') 
    system(['mkdir ' dirname]);
end

disp(['Saving files to directory ' dirname '/'])
for i=1:nsub
    filename=[dirname '/' cfg.subs{i} '.mat'];
    data=cors(:,:,i);
    save(filename,'data')
end
% save group average
filename=[dirname '/groupAverage.mat'];
save(filename,'avg_cors')

cfg.cors=cors;
disp('Done!')
fprintf('\n')