function calculate_indepRef_twISC(cfg)
% Calculates individual ISC maps with respect to a reference group
% 
% Usage:
%   calculate_twISC(cfg);
%
%   Input:
%   cfg.subs = IDs (e.g. 'EPHE602') of subjects of interest
% 
%   cfg.refs = IDs (e.g. 'EPHE602') of reference subjects
% 
%   cfg.mask = nifti file of binary mask; default is the MNI152 whole brain
%   mask
% 
%   cfg.res = resolution of the nifti files;  can be '2mm', '4mm', '6mm',
%   '8mm', '16mm', or '32mm'
%
%   cfg.condSubs = data condition to use for subjects of interest; options baseline 'BL' or follow-up 'FU'
%
%   cfg.condRefs = data condition to use for reference subjects; options baseline 'BL' or follow-up 'FU'
%
%   cfg.toi = time points of interest in a vector
% 
%   cft.outdir = label of the output folder where the individual ISC nifitis are saved

%%
%   Output:
%   
%%
% for usage example, see /m/nbe/scratch/psykoosi/scripts/demo_calculate_twISC.m

%% Input validation

if ~ismember(cfg.subs,{'ConsN30','PatsN36','BothN66'})
    error('cfg.subs should be either ''ConsN30'', ''PatsN36'' or ''BothN66''!')
end
if ~ismember(cfg.refs,{'ConsN15','ConsN30','PatsN36'})
    error('cfg.refs should be either ''ConsN15'', ''ConsN30'' or ''PatsN36''!')
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
    error(['Could not find mask: ' cfg.mask])
end
if ~isvector(cfg.toi)
    error('cfg.toi must be a vector of numbers, e.g. [1 2 3 4 5]!')
end
if ~ismember(cfg.useMeanOverRefs,[0 1])
    error('cfg.useMeanOverRefs must be either 0 or 1')
end
if ~ismember(cfg.useNonSpatialSmoothedData,[0 1])
    error('cfg.useNonSpatialSmoothedData must be either 0 or 1')
end

cfg.toi=sort(cfg.toi);


disp(['Calculating time windowed ISC for ' cfg.subs ' with respect to ' cfg.refs])
disp(['using ' cfg.condSubs ' data for ' cfg.subs ' and ' cfg.condRefs ' data for ' cfg.refs])
disp(['using time points: ' num2str(cfg.toi)])
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
            consN30{sample}=codes{codei};
        elseif mode==3
            consN15{sample}=codes{codei};
        end
    end
end

% Set IDs for subjects of interest and reference subjects
switch cfg.subs
    case 'PatsN36'
        cfg.subs=patsN36;
    case 'ConsN30'
        cfg.subs=consN30;
    case 'BothN66'
        cfg.subs=[patsN36 consN30];
end

switch cfg.refs
    case 'PatsN36'
        cfg.refs=patsN36;
    case 'ConsN30'
        cfg.refs=consN30;
    case 'ConsN15'
        cfg.refs=consN15;
end

%% Load mask
addpath(genpath('/m/nbe/scratch/psykoosi/scripts'));

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

% Get the number of time points, voxels, subjects
notps=size(temp(:,inmask),1);

assert(cfg.toi(end)<=notps)
novox=size(temp(:,inmask),2);
nosubs=length(cfg.subs);
norefs=length(cfg.refs);
notoi=length(cfg.toi);

% Preallocate data matrices with the dimensions of time points of interest x subjects x voxels 
allsubs=zeros(notoi,nosubs,novox); 
allrefs=zeros(notoi,norefs,novox);

disp('Loading brain data of subjects of interest...')
for iter=1:nosubs
    disp(['Subject ' cfg.subs{iter} ' - ' num2str(iter) ' out of ' num2str(nosubs)])
    if cfg.useNonSpatialSmoothedData
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.condSubs '/' cfg.subs{iter} '/epi_noSpatialSmoothing_' cfg.res '.nii']);
    else
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.condSubs '/' cfg.subs{iter} '/epi_preprocessed_' cfg.res '.nii']);
    end
    temp=permute(nii.img,[4 1 2 3]);
    allsubs(:,iter,:)=zscore(temp(cfg.toi,inmask));
%     allsubs(:,i,:)=temp(:,inmask);
end
fprintf('\n')

disp('Loading brain data of reference subjects...')
for iter=1:norefs
    disp(['Subject ' cfg.refs{iter} ' - ' num2str(iter) ' out of ' num2str(norefs)])
    if cfg.useNonSpatialSmoothedData
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.condRefs '/' cfg.refs{iter} '/epi_noSpatialSmoothing_' cfg.res '.nii']);
    else
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.condRefs '/' cfg.refs{iter} '/epi_preprocessed_' cfg.res '.nii']);
    end
    temp=permute(nii.img,[4 1 2 3]);
    allrefs(:,iter,:)=zscore(temp(cfg.toi,inmask));
%     allrefs(:,i,:)=temp(:,inmask);
end
fprintf('\n')

%% Calculate ISC

% Preallocate the correlation matrix
cors=zeros(novox,nosubs,norefs);
    
disp('Calculating ISC...')
for voxi=1:novox
    if mod(voxi,1000)==0 % Show the status every 1000 voxels
        disp([num2str(voxi) '/' num2str(novox) ' voxels'])
    end
       
    % Calculate the correlation over time of all subjects of interest to all reference subjects
    cors(voxi,:,:)=corr(squeeze(allsubs(:,:,voxi)),squeeze(allrefs(:,:,voxi)));
  
end
fprintf('\n')

%% Average to obtain indivual ISC maps and save them

% Check if output directory exists; if not, create it
dirname=['/m/nbe/scratch/psykoosi/ISC/indepRef_ISC/' cfg.outdir];
if ~exist(dirname,'dir') 
    system(['mkdir -p ' dirname]);
end

% Replace possible NaN values with zeros.
cors(isnan(cors))=0;

if cfg.useMeanOverRefs
    disp(['Saving files to directory ' dirname '/'])
    for iter=1:nosubs
        % Check if any subject of interest was also a reference subject; 
        % if so, the cors matrix has three dimensions, otherwise only two
        if size(cors,3)>0
            newbrain=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));
            % We calculated an individual average reference for each
            % subject of interest, hence cors(:,iter,iter)
            newbrain(inmask)=cors(:,iter,iter);
            filename=[dirname '/' cfg.subs{iter} '.nii'];
            save_nii(make_nii(newbrain),filename);
            nii=fixOriginator(filename,mask);
            save_nii(nii,filename);   
        else
            newbrain=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));
            newbrain(inmask)=cors(:,iter);
            filename=[dirname '/' cfg.subs{iter} '.nii'];
            save_nii(make_nii(newbrain),filename);
            nii=fixOriginator(filename,mask);
            save_nii(nii,filename);
        end
    end
else
    % Average and save individual ISC maps as nifti files
    disp(['Averaging and saving files to directory ' dirname '/'])
    for iter=1:nosubs
        % Check if the subject of interest is also a reference subject
        if ismember(cfg.subs{iter},cfg.refs)
            disp(['Subject ' cfg.subs{iter} ' - ' num2str(iter) ' out of ' num2str(nosubs)])
            % Find out the indices of all other but the subject of interest in the reference subject list
            [~,~,inds] = setxor(cfg.subs{iter},cfg.refs);
            % Take the average across all reference subjects (first Fisher's
            % z-tranform, then averaging, then back to correlation scale)
            avg_cors=squeeze(tanh(mean(atanh(cors(:,iter,inds)),3)));
            newbrain=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));
            newbrain(inmask)=avg_cors;
            filename=[dirname '/' cfg.subs{iter} '.nii'];
            save_nii(make_nii(newbrain),filename);
            nii=fixOriginator(filename,mask);
            save_nii(nii,filename);
        else
            disp(['Subject ' cfg.subs{iter} ' - ' num2str(iter) ' out of ' num2str(nosubs)])
            % Take the average across all reference subjects (first Fisher's
            % z-tranform, then averaging, then back to correlation scale)
            avg_cors=squeeze(tanh(mean(atanh(cors(:,iter,:)),3)));
            newbrain=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));
            newbrain(inmask)=avg_cors;
            filename=[dirname '/' cfg.subs{iter} '.nii'];
            save_nii(make_nii(newbrain),filename);
            nii=fixOriginator(filename,mask);
            save_nii(nii,filename);
        end
    end
end

disp('Done!')
fprintf('\n')