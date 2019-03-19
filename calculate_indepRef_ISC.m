function cfg=calculate_indepRef_ISC(cfg)
% Calculates individual ISC maps with respect to a reference group
% 
% Usage:
%   cfg = calculate_indepRef_ISC(cfg);
%
%   Input:
%   cfg.subs = IDs (e.g. 'EPHE602') of subjects of interest (these are the
%   subjects the individual maps are calculated for)
% 
%   cfg.refs = IDs of subjects in the reference group (i.e. subjects with
%   respect to which the indivual ISC maps are calculated)
% 
%   cfg.condSubs = data condition to use for subjects of interest; options baseline 'BL' or follow-up 'FU'
%
%   cfg.condRefs = data condition to use for reference subjects; options baseline 'BL' or follow-up 'FU'
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
%   cfg.cors = temporal correlation of all subjects of interest to all
%   reference subjects (matrix with dimensions of voxels x subjects-of-interest x references

%%
% for usage example, see /m/nbe/scratch/psykoosi/scripts/demo_calculate_indepRef_ISC.m

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
if ~ismember(cfg.useMeanOverRefs,[0 1])
    error('cfg.useMeanOverRefs must be either 0 or 1!')
end

disp(['Calculating individual ISC maps for ' num2str(length(cfg.subs)) ' subjects of interest with respect to ' num2str(length(cfg.refs)) ' reference subjects'])
disp(['using ' cfg.condSubs ' data for subjects of interest and ' cfg.condRefs ' data for reference subjects'])
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
notps=size(temp(:,inmask),1);
novox=size(temp(:,inmask),2);
nosubs=length(cfg.subs);
norefs=length(cfg.refs);

% Preallocate data matrices with the dimensions of time points x subjects x voxels 
allsubs=zeros(notps,nosubs,novox); 
allrefs=zeros(notps,norefs,novox);

disp('Loading brain data of subjects of interest...')
for i=1:nosubs
    disp(['Subject ' cfg.subs{i} ' - ' num2str(i) ' out of ' num2str(nosubs)])
    if cfg.useNonSpatialSmoothedData
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.condSubs '/' cfg.subs{i} '/epi_noSpatialSmoothing_' cfg.res '.nii']);
    else
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.condSubs '/' cfg.subs{i} '/epi_preprocessed_' cfg.res '.nii']);
    end
    temp=permute(nii.img,[4 1 2 3]);
    allsubs(:,i,:)=zscore(temp(:,inmask));
%     allsubs(:,i,:)=temp(:,inmask);
end
fprintf('\n')
disp('Loading brain data of reference subjects...')
for i=1:norefs
    disp(['Subject ' cfg.refs{i} ' - ' num2str(i) ' out of ' num2str(norefs)])
    if cfg.useNonSpatialSmoothedData
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.condRefs '/' cfg.refs{i} '/epi_noSpatialSmoothing_' cfg.res '.nii']);
    else
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.condRefs '/' cfg.refs{i} '/epi_preprocessed_' cfg.res '.nii']);
    end
    temp=permute(nii.img,[4 1 2 3]);
    allrefs(:,i,:)=zscore(temp(:,inmask));
%     allrefs(:,i,:)=temp(:,inmask);
end

fprintf('\n')

%% Take average over the reference group if cfg.useMeanOverRefs==1
if cfg.useMeanOverRefs
    clear temp
    disp('Taking average over the reference group and using that to calculate ISC');
    % Check if any subject of interest is also a reference subject
    if any(ismember(cfg.subs,cfg.refs))
        for i=1:nosubs
            % Check if the subject of interest ia also a reference subject
            if ismember(cfg.subs{i},cfg.refs)
                % Find out the indices of all other but the subject of interest in the reference subject list
                [~,~,inds] = setxor(cfg.subs{i},cfg.refs);
                % Average reference groups so that the subject of interest is not a member
                temp(:,i,:)=mean(allrefs(:,inds,:),2);
            else
                temp(:,i,:)=mean(allrefs,2);
            end
        end
        allrefs=temp;
        norefs=nosubs;
    else
        % If none of the subjects of interest is a reference subject,
        % we can just average over the whole group
        allrefs=mean(allrefs,2);
        norefs=1;
    end
    fprintf('\n')
end

%% Calculate ISC
% Preallocate the correlation matrix
cors=zeros(novox,nosubs,norefs);
    
disp('Calculating ISC...')
for voxi=1:novox
    if mod(voxi,1000)==0 % Show the status every 1000 voxels
        disp([num2str(voxi) '/' num2str(novox) ' voxels'])
    end
       
    % Calculate the correlation over time of all subjects-of-interest to all reference subjects
    cors(voxi,:,:)=corr(squeeze(allsubs(:,:,voxi)),squeeze(allrefs(:,:,voxi)));
  
end

fprintf('\n')

%% Average to obtain indivual ISC maps

% Take the average across all reference subjects (first Fisher's
% z-tranform, then averaging, then back to correlation scale)
% avg_cors=squeeze(tanh(mean(atanh(cors),3)));
% or without the Fisher's z-transform
%     avg_cors=squeeze(mean(cors,3));

% Replace possible NaN values with zeros.
cors(isnan(cors))=0;

% Check if output directory exists; if not, create it
dirname=['/m/nbe/scratch/psykoosi/ISC/indepRef_ISC/' cfg.outdir];
if ~exist(dirname,'dir') 
    system(['mkdir ' dirname]);
end

if cfg.useMeanOverRefs
    disp(['Saving files to directory ' dirname '/'])
    for i=1:nosubs
        % Check if any subject of interest was also a reference subject; 
        % if so, the cors matrix has three dimensions, otherwise only two
        if any(ismember(cfg.subs,cfg.refs))
            newbrain=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));
            % We calculated an individual average reference for each
            % subject of interest, hence the cors(:,i,i)
            newbrain(inmask)=cors(:,i,i);
            filename=[dirname '/' cfg.subs{i} '.nii'];
            save_nii(make_nii(newbrain),filename);
            nii=fixOriginator(filename,mask);
            save_nii(nii,filename);   
        else
            newbrain=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));
            newbrain(inmask)=cors(:,i);
            filename=[dirname '/' cfg.subs{i} '.nii'];
            save_nii(make_nii(newbrain),filename);
            nii=fixOriginator(filename,mask);
            save_nii(nii,filename);
        end
    end
else
    % Average and save individual ISC maps as nifti files
    disp(['Averaging and saving files to directory ' dirname '/'])
    for i=1:nosubs
        % Check if the subject of interest is also a reference subject
        if ismember(cfg.subs{i},cfg.refs)
            disp(['Subject ' cfg.subs{i} ' - ' num2str(i) ' out of ' num2str(nosubs)])
            % Find out the indices of all other but the subject of interest in the reference subject list
            [~,~,inds] = setxor(cfg.subs{i},cfg.refs);
            % Take the average across all reference subjects (first Fisher's
            % z-tranform, then averaging, then back to correlation scale)
            avg_cors=squeeze(tanh(mean(atanh(cors(:,i,inds)),3)));
            newbrain=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));
            newbrain(inmask)=avg_cors;
            filename=[dirname '/' cfg.subs{i} '.nii'];
            save_nii(make_nii(newbrain),filename);
            nii=fixOriginator(filename,mask);
            save_nii(nii,filename);
        else
            disp(['Subject ' cfg.subs{i} ' - ' num2str(i) ' out of ' num2str(nosubs)])
            % Take the average across all reference subjects (first Fisher's
            % z-tranform, then averaging, then back to correlation scale)
            avg_cors=squeeze(tanh(mean(atanh(cors(:,i,:)),3)));
            newbrain=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));
            newbrain(inmask)=avg_cors;
            filename=[dirname '/' cfg.subs{i} '.nii'];
            save_nii(make_nii(newbrain),filename);
            nii=fixOriginator(filename,mask);
            save_nii(nii,filename);
        end
    end
end

cfg.cors=cors;
disp('Done!')
fprintf('\n')