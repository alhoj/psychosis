function cfg=calculate_BL_FU_intrasubject_corr(cfg)
% Calculates individual ISC maps with respect to a reference group
% 
% Usage:
%   cfg = calculate_BL_FU_intrasubject_corr(cfg);
%
%   Input:
%   cfg.subs = IDs (e.g. 'EPHE602') of subjects
% 
%   cfg.mask = nifti file of binary mask; default is the MNI152 whole brain
%   mask
% 
%   cfg.res = resolution of the nifti files;  can be '2mm', '4mm', '6mm',
%   '8mm', '16mm', or '32mm'
% 
%   cft.outdir = label of the output folder where the correlation map nifitis are saved

%%
%   Output:
%   cfg.cors = 

%%
% for usage example, see /m/nbe/scratch/psykoosi/scripts/

%% Input validation

if ~iscell(cfg.subs)
    error('cfg.subs should contain a cell array!')
end
if ~ismember(cfg.res,{'2mm', '4mm', '6mm', '8mm', '16mm', '32mm'}) 
    error('cfg.res has to be either ''2mm'', ''4mm'', ''6mm'', ''8mm'', ''16mm'', or ''32mm''')
end
if ~isempty(cfg.mask) && ~isfile(cfg.mask)
    error(['could not find mask: ' cfg.mask])
end

disp(['Calculating BL vs FU intrasubject correlation maps for ' num2str(length(cfg.subs)) ' subjects'])
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

%% Calculate correlation

% Check if output directory exists; if not, create it
dirname=['/m/nbe/scratch/psykoosi/ISC/' cfg.outdir];
if ~exist(dirname,'dir') 
    system(['mkdir ' dirname]);
end

nosubs=length(cfg.subs);
novox=length(inmask);

disp('Calculating correlation...')
for i=1:nosubs 
    disp(['Subject ' cfg.subs{i} ' - ' num2str(i) ' out of ' num2str(nosubs)])
    niiBL=load_nii(['/m/nbe/scratch/psykoosi/data/BL/' cfg.subs{i} '/epi_preprocessed_' cfg.res '.nii']);
    temp=permute(niiBL.img,[4 1 2 3]);
    dataBL=zscore(temp(:,inmask));
    niiFU=load_nii(['/m/nbe/scratch/psykoosi/data/FU/' cfg.subs{i} '/epi_preprocessed_' cfg.res '.nii']);
    temp=permute(niiFU.img,[4 1 2 3]);
    dataFU=zscore(temp(:,inmask));
    % Preallocate correlation matrix
    cors=zeros(novox,1);
    % Calculate correlation for each voxel
    for voxi=1:novox
        if mod(voxi,1000)==0 % Show the status every 1000 voxels
            disp([num2str(voxi) '/' num2str(novox) ' voxels'])
        end
        cors(voxi)=corr(dataBL(:,voxi),dataFU(:,voxi));
    end
    newbrain=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));
    newbrain(inmask)=cors;
    filename=[dirname '/' cfg.subs{i} '.nii'];
    save_nii(make_nii(newbrain),filename);
    nii=fixOriginator(filename,mask);
    save_nii(nii,filename);
end

fprintf('\n')
