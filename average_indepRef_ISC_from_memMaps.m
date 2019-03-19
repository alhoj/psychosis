function cfg=average_indepRef_ISC_from_memMaps(cfg)
% calculates individual ISC maps with respect to a reference group from ISC
% toolbox memMaps
% 
% Usage:
%   cfg = average_indepRef_ISC_from_memMaps(cfg);
%
%   Input:
%   cfg.subs = IDs (e.g. 'EPHE602') of all the subjects that were input to
%   ISC toolbox
% 
%   cfg.refs = IDs of subjects in the reference group (i.e. subjects with
%   respect to which the indivual ISC maps are calculated)
% 
%   cfg.indir = name of the folder where the memMaps.mat structure is located  
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
%   cfg.data = ISC correlation matrix from memMaps

%%
% for usage example, see /m/nbe/scratch/psykoosi/scripts/demo_average_indepRef_ISC_from_memMaps.m

%% Input validation

if ~iscell(cfg.subs)
    error('cfg.subs should contain a cell array!')
end
if ~iscell(cfg.subsoi)
    error('cfg.subsoi should contain a cell array!')
end
if ~iscell(cfg.refs)
    error('cfg.refs should contain a cell array!')
end
if ~all(ismember(cfg.subsoi,cfg.subs))
    error('subjects in cfg.subsoi must be also in cfg.subs!')
end
if ~all(ismember(cfg.refs,cfg.subs))
    error('subjects in cfg.refs must be also in cfg.subs!')
end
if ~ismember(cfg.res,{'2mm', '4mm', '6mm', '8mm', '16mm', '32mm'}) 
    error('cfg.res has to be either ''2mm'', ''4mm'', ''6mm'', ''8mm'', ''16mm'', or ''32mm''')
end
if ~isempty(cfg.mask) && ~isfile(cfg.mask)
    error(['could not find mask: ' cfg.mask])
end


%% Determine the number of subjects (total, subjects of interest, and reference subjects)
nosubs=length(cfg.subs); % total number of subjects
norefs=length(cfg.refs); % number of reference subjects
nosubsoi=length(cfg.subsoi); % number of subjects of interest
disp(['Averaging individual ISC maps for ' num2str(nosubsoi) ' subjects with respect to ' num2str(norefs) ' reference subjects'])
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

%% Load ISC data from memMaps

disp('Loading ISC data...')
load([cfg.indir 'memMaps.mat'])
temp=double(memMaps.cormatMap.whole.band0.Session1.cor.Data.xyzc);
temp=permute(temp,[4 1 2 3]);
data=temp(:,inmask)';
if (nchoosek(length(cfg.subs),2) ~= size(data,2))
    error('Different number of subjects in the specified subject list and in ISC data from memMaps!')
end
fprintf('\n')

%% Average to obtain indivual ISC maps

% Check if output directory exists; if not, create it
dirname=['/m/nbe/scratch/psykoosi/ISC/indepRef_ISC/' cfg.outdir];
if ~exist(dirname,'dir') 
    system(['mkdir ' dirname]);
end

% Find indices of reference subjects and subjects of interest
inds_refs=find(contains(cfg.subs,cfg.refs))';
inds_subsoi=find(contains(cfg.subs,cfg.subsoi))';

count=1;
for i=inds_subsoi
    disp(['Averaging ISC of subject ' cfg.subs{i} ' - ' num2str(count) ' out of ' num2str(nosubsoi)])
    % Create a matrix of zeros corresponding to the memMaps correlation
    % matrix and set ones for the subject of interest
    subjectmat=zeros(nosubs); subjectmat(i,:)=1; subjectmat(:,i)=1;
    % Find all indices of the subject of interest in the upper triangular part of the matrix
    matInds_suboi=find(subjectmat(find(triu(ones(nosubs),1)))==1);
    % Find the indices the subject of interest shares with the reference subjects
    matInds_refs=[matInds_suboi(inds_refs(find(inds_refs<i)));matInds_suboi(inds_refs(find(inds_refs>i))-1)];
    % Take the average across all reference subjects (first Fisher's
    % z-tranform, then averaging, then back to correlation scale)
    avg_ISC=tanh(mean(atanh(data(:,matInds_refs)),2));
    % Replace possible NaN values with zeros.
    avg_ISC(isnan(avg_ISC))=0;
    newbrain=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));
    newbrain(inmask)=avg_ISC;
    % Save individual ISC as nifti
    filename=[dirname '/' cfg.subs{i} '.nii'];
    save_nii(make_nii(newbrain),filename);
    nii=fixOriginator(filename,mask);
    save_nii(nii,filename);
    count=count+1;
end

cfg.data=data;
disp('Done!')
