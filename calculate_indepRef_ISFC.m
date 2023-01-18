function cfg=calculate_indepRef_ISFC(cfg)
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

if ~isempty(cfg.mask) && ~isfile(cfg.mask)
    error(['could not find mask: ' cfg.mask])
end
if ~isempty(cfg.atlas) && ~isfile(cfg.atlas)
    error(['could not find atlas: ' cfg.atlas])
end
if ~ismember(cfg.roiTC,{'mean','pca'})
    error('cfg.roiTC must be either ''mean'' or ''pca''!')
end

fprintf('\n')
%% Read IDs for the subject groups

codes=importdata('subject_codes.txt'); % Import the subjects code text to split them into patients, controls and reference
% Split codes into patient codes, control codes and reference codes

mode=0; %  1 for patientsN36, 2 for controlsN29, 3 for referenceN15, 4 for patientsN28 (not in FU), and 5 for larger sample patientsN57 (only in BL)
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
        elseif mode==4
            patsN28{sample}=codes{codei};
        elseif mode==5
            patsN77{sample}=codes{codei};
        elseif mode==6
            patsN54{sample}=codes{codei};
        end
    end
end

% Set IDs for subjects of interest
switch cfg.subs
    case 'patsN36'
        cfg.subs=patsN36;
    case 'consN29'
        cfg.subs=consN29;
    case 'bothN65'
        cfg.subs=[patsN36 consN29];
    case 'consN15'
        cfg.subs=consN15;
    case 'patsN28'
        cfg.subs=patsN28;
    case 'allN108'
        cfg.subs=[patsN36 consN29 consN15 patsN28];
    case 'patsN77'
        cfg.subs=patsN77;
    case 'patsN54'
        cfg.subs=patsN54;
end

% Set IDs for reference subjects
refs=cfg.refs;
switch cfg.refs
    case 'patsN36'
        cfg.refs=patsN36;
    case 'consN29'
        cfg.refs=consN29;
    case 'bothN65'
        cfg.refs=[patsN36 consN29];
    case 'consN15'
        cfg.refs=consN15;
    case 'patsN28'
        cfg.refs=patsN28;
end

disp(['Calculating inter-subject functional connectivity for ' num2str(length(cfg.subs)) ' subjects with respect to ' num2str(length(cfg.refs)) ' references based on Pearson''s correlation'])
fprintf('\n')

%% Load mask and atlas
addpath(genpath('/m/nbe/scratch/psykoosi/scripts'))

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

test_nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.cond '/' cfg.subs{1} '/epi_preprocessed_' cfg.res '.nii']);

% Put time first, so that we can acess the x y z with the 1-D indices
temp=permute(test_nii.img,[4 1 2 3]);

% Get the number of time points, voxels, subjects-of-interest, and reference subjects
if isempty(cfg.toi)
    ntps=size(temp(:,inmask),1);
    cfg.toi=1:ntps;   
else
    ntps=length(cfg.toi);
end
nvox=size(temp(:,inmask),2);
nsub=length(cfg.subs);
nref=length(cfg.refs);

% Preallocate data matrix with the dimensions of time points x voxels x subjects  
TCvox_subs=zeros(ntps,nvox,nsub); 
TCvox_refs=zeros(ntps,nvox,nref); 

disp('Loading brain data...')
disp('Subjects of interest')
for subi=1:nsub   
    disp([cfg.subs{subi} ' - ' num2str(subi) ' out of ' num2str(nsub)])
    nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.cond '/' cfg.subs{subi} '/epi_preprocessed_' cfg.res '.nii']);
    temp=permute(nii.img,[4 1 2 3]);
    TCvox_subs(:,:,subi)=zscore(temp(cfg.toi,inmask));
    
    if ~isempty(cfg.atlas)
        for j=1:nroi
            roi=rois(j); % all rois of the atlas might not be within the mask, hence this
            if isequal(cfg.roiTC,'pca')
                % Principal component timecourses for ROIs
                [~,PCTC,~] = pca(TCvox_subs(:,atlas(inmask)==roi,subi));
                TCroi_subs(:,j,subi) = PCTC(:,1);
            else
                % Mean timecourses for ROIs
                TCroi_subs(:,j,subi) = mean(TCvox_subs(:,atlas(inmask)==roi,subi),2);
            end
        end
    end
end
fprintf('\n')

disp('Reference subjects')
for refi=1:nref
    disp([cfg.refs{refi} ' - ' num2str(refi) ' out of ' num2str(nref)])
    if isequal(refs,'consN15') || isequal(refs,'patsN28')
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/BL/' cfg.refs{refi} '/epi_preprocessed_' cfg.res '.nii']);
    else
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.cond '/' cfg.refs{refi} '/epi_preprocessed_' cfg.res '.nii']);
    end
    temp=permute(nii.img,[4 1 2 3]);
    TCvox_refs(:,:,refi)=zscore(temp(cfg.toi,inmask));
    
    if ~isempty(cfg.atlas)
        for j=1:nroi
            roi=rois(j); % all rois of the atlas might not be within the mask, hence this
            if isequal(cfg.roiTC,'pca')
                % Principal component timecourses for ROIs
                [~,PCTC,~] = pca(TCvox_refs(:,atlas(inmask)==roi,refi));
                TCroi_refs(:,j,refi) = PCTC(:,1);
            else
                % Mean timecourses for ROIs
                TCroi_refs(:,j,refi) = mean(TCvox_refs(:,atlas(inmask)==roi,refi),2);
            end
        end
    end
end
fprintf('\n')

%% Calculate ISFC
    
disp('Calculating ISFC...')
if isempty(cfg.atlas)
    cors=zeros(nvox,nvox,nsub,nref);
else
    nroi=size(TCroi_subs,2);
    cors=zeros(nroi,nroi,nsub,nref);
end

for s1=1:nsub
    disp(s1)
    for s2=1:nref
%             disp(s2)
        if isempty(cfg.atlas)       
            % Calculate functional correlation over time between all voxels
            cors(:,:,s1,s2)=corr(TCvox_subs(:,:,s1),TCvox_refs(:,:,s2));
        else
            % Calculate functional correlation over time between all rois
            cors(:,:,s1,s2)=corr(TCroi_subs(:,:,s1),TCroi_refs(:,:,s2));
        end
    end
end

if cfg.symmetrize
    cors=(cors+permute(cors,[2 1 3 4]))/2;
end

% Replace possible NaN values with zeros.
cors(isnan(cors))=0;

fprintf('\n')

%% Save subject-wise and group ISFC matrices

% Take the average across references (first Fisher's
% z-tranform, then averaging, then back to correlation scale)
cors=squeeze(tanh(mean(atanh(cors),4)));

% Take group average
avg_cors=squeeze(tanh(mean(atanh(cors),3)));
% Replace possible NaN values with zeros.
avg_cors(isnan(avg_cors))=0;

% Check if output directory exists; if not, create it
dirname=['/m/nbe/scratch/psykoosi/ISFC/' cfg.outdir];
if ~exist(dirname,'dir') 
    system(['mkdir -p ' dirname]);
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