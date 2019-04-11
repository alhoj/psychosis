function cfg=calculate_FC_generic(cfg)
% Calculates voxel-to-voxel functional connectivity
% 
% Usage:
%   cfg = calculate_FC(cfg);
%
%   Input:
%   cfg.subs = subject IDs
% 
%   cfg.indir = path to 4D niftis of subjects defined in cfg.subs
%
%   cfg.mask = nifti file of binary mask; default is the MNI152 whole brain
%   mask
% 
%   cfg.res = resolution of the nifti files;  can be '2mm', '4mm', '6mm',
%   '8mm', '16mm', or '32mm'
%
%   cfg.atlas = ROI atlas to be used; leave empty for voxel-to-voxel connectivity
%
%   cfg.roiTC = how to extract the roi timecourse; options 'PCA' or 'mean'
% 
%   cft.outdir = label of the output folder where the individual ISC nifitis are saved

%%
%   Output:
%   cfg.cors = correlation matrix (voxels x voxels x subjects)


%% Input validation

if ~iscell(cfg.subs)
    error('cfg.subs should contain a cell array!')
end
if ~isempty(cfg.mask) && ~isfile(cfg.mask)
    error(['could not find mask: ' cfg.mask])
end
if ~isempty(cfg.atlas) && ~isfile(cfg.atlas)
    error(['could not find atlas: ' cfg.atlas])
end
if ~ismember(cfg.roiTC,{'mean','PCA'})
    error('cfg.roiTC must be either ''mean'' or ''PCA''!')
end

if isempty(cfg.atlas)
    disp(['Calculating voxel-to-voxel functional connectivity for ' num2str(length(cfg.subs)) ' subjects based on Pearson''s correlation'])
else
    disp(['Calculating functional connectivity for ' num2str(length(cfg.subs)) ' subjects between ROIs from atlas ''' cfg.atlas ''' based on Pearson''s correlation'])
end
fprintf('\n')

%% Load mask and atlas

disp('Loading mask...')
mask=load_nii(cfg.mask);
inmask=find(mask.img);

if ~isempty(cfg.atlas)
    disp('Loading atlas...')
    atlas=load_nii(cfg.atlas);
    assert(isequal(length(atlas.img(:)),length(mask.img(:))),'mask and atlas voxel resolutions do not match!')
    atlas=double(atlas.img).*double(mask.img);
    rois=nonzeros(unique(atlas));
    nroi=length(rois);
end

fprintf('\n')
%% Load brain data
% Load an example NIFTI file to find the dimensions needed for
% preallocating the matrices (this will speed things up)
%test_nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.cond '/' cfg.subs{1} '/epi_preprocessed_' cfg.res '.nii']);
test_nii=load_nii([cfg.indir '/' cfg.subs{1} '.nii']);

% Put time first, so that we can acess the x y z with the 1-D indices
temp=permute(test_nii.img,[4 1 2 3]);

% Get the number of time points, voxels, subjects-of-interest, and
% reference subjects
ntps=size(temp(:,inmask),1);
nvox=size(temp(:,inmask),2);
nsub=length(cfg.subs);

% Preallocate data matrix with the dimensions of time points x voxels x subjects  
TCvox=zeros(ntps,nvox,nsub); 

disp('Loading brain data...')
for subi=1:nsub
    disp(['File ' cfg.subs{subi} ' - ' num2str(subi) ' out of ' num2str(nsub)])
%     nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.cond '/' cfg.subs{subi} '/epi_preprocessed_' cfg.res '.nii']);
    nii=load_nii([cfg.indir '/' cfg.subs{subi} '.nii']);
    temp=permute(nii.img,[4 1 2 3]);
    TCvox(:,:,subi)=zscore(temp(:,inmask));
    
    for j=1:nroi
        roi=rois(j); % all rois of the atlas might not be within the mask, hence this
        if isequal(cfg.roiTC,'PCA')
            % Principal component timecourses for ROIs
            [~,PCTC,~] = pca(TCvox(:,atlas(atlas>0)==roi));
            TCroi(:,roi,subi)=PCTC(:,1);
        else
            % Mean timecourses for ROIs
            TCroi(:,roi,subi)=mean(TCvox(:,atlas(atlas>0)==roi),2);
        end
    end
end
fprintf('\n')

%% Calculate FC
    
disp('Calculating FC...')
for subi=1:nsub
    disp(subi)
    if isempty(cfg.atlas)       
        % Calculate functional correlation over time between all voxels
        cors(:,:,subi)=corr(squeeze(TCvox(:,:,subi))); 
    else
        % Calculate functional correlation over time between all rois
        cors(:,:,subi)=corr(squeeze(TCroi(:,:,subi)));
    end
end
    
% Replace possible NaN values with zeros.
cors(isnan(cors))=0;

fprintf('\n')

%% Save subject-wise and group FC matrices

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
for subi=1:nsub
    filename=[dirname '/' cfg.subs{subi} '.mat'];
    data=cors(:,:,subi);
    save(filename,'data')
end
% save group average
filename=[dirname '/groupAverage.mat'];
save(filename,'avg_cors')

cfg.cors=cors;
disp('Done!')
fprintf('\n')