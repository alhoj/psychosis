function cfg=calculate_ISFC(cfg)
% Calculates voxel-to-voxel functional connectivity
% 
% Usage:
%   cfg = calculate_ISFC(cfg);
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
if ~ismember(cfg.roiTC,{'mean','pca'})
    error('cfg.roiTC must be either ''mean'' or ''pca''!')
end

disp(['Calculating inter-subject functional connectivity for ' num2str(length(cfg.subs)) ' subjects based on Pearson''s correlation'])
fprintf('\n')

%% Read IDs for the subject groups

codes=importdata('subject_codes.txt'); % Import the subjects code text to split them into patients, controls and reference
% Split codes into patient codes, control codes and reference codes

mode=0; %  1 for patients, 2 for controls, 3 for reference
sample=0; % Increasing index for each subject in the same group
for codei=1:length(codes)
    
    if ~strcmp(codes{codei}(1:2),'EP') % If it is not a subject code
        codes{codei};
        sample=0; % Set the indexing number to zero
        mode=mode+1; % Increase the mode by 1 to get to the next category
    else
        sample=sample+1;  % Increase the index of the subject in this category
        if mode==1
            patsN36{sample}=codes{codei};
        elseif mode==2
            consN29{sample}=codes{codei};
        elseif mode==3
            consN15{sample}=codes{codei};
        end
    end
end

% Set IDs for reference subjects
switch cfg.refs
    case 'PatsN36'
        cfg.refs=patsN36;
    case 'ConsN29'
        cfg.refs=consN29;
    case 'ConsN15'
        cfg.refs=consN15;
end


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
        if isequal(cfg.roiTC,'pca')
            % Principal component timecourses for ROIs
%             [~,PCTC,~] = pca(TCvox(:,atlas(atlas>0)==roi,subi));
            [~,PCTC,~] = pca(TCvox(:,atlas(inmask)==roi,subi));
            TCroi(:,roi,subi) = PCTC(:,1);
        else
            % Mean timecourses for ROIs
%             TCroi(:,roi,subi) = mean(TCvox(:,atlas(atlas>0)==roi,subi),2);
            TCroi(:,roi,subi) = mean(TCvox(:,atlas(inmask)==roi,subi),2);
        end
    end
end
fprintf('\n')

%% Calculate ISFC
    
disp('Calculating ISFC...')
if isempty(cfg.atlas)
    cors=zeros(size(TCvox,2),size(TCvox,2),nsub*(nsub-1));
else
    cors=zeros(size(TCroi,2),size(TCroi,2),nsub*(nsub-1));
end
idx=0;
for s1=1:nsub
    disp(s1)
    for s2=1:nsub
        if s1~=s2
            idx=idx+1;
%             disp(s2)
            if isempty(cfg.atlas)       
                % Calculate functional correlation over time between all voxels
                cors(:,:,idx)=corr(TCvox(:,:,s1),TCvox(:,:,s2));
            else
                % Calculate functional correlation over time between all rois
                cors(:,:,idx)=corr(TCroi(:,:,s1),TCroi(:,:,s2));
            end
        end
    end
%     fprintf('done\n');
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
dirname=['/m/nbe/scratch/psykoosi/ISFC/' cfg.outdir];
if ~exist(dirname,'dir') 
    system(['mkdir ' dirname]);
end

% disp(['Saving files to directory ' dirname '/'])
% for subi=1:nsub
%     filename=[dirname '/' cfg.subs{subi} '.mat'];
%     data=cors(:,:,subi);
%     save(filename,'data')
% end
% save group average
filename=[dirname '/groupAverage.mat'];
save(filename,'avg_cors')

cfg.cors=cors;
disp('Done!')
fprintf('\n')