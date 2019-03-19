function cfg=regress_PCs_generic(cfg)
% Regress out BOLD signal principal components of one group of subjects
% from another group of subjects
% 
% Usage:
%   cfg = regress_PCs(cfg);
%
%   Input:
%   cfg.sub = data paths of subject of interest
% 
%   cfg.refs = data paths of the subjects whose data are used 
%   to calculate the PCs 
%
%   cfg.mask = path to nifti file of binary mask
%
%   cft.outdir = path of the output folder where the cleaned niftis are saved

%%
%   Output:
%   cfg.dataReg = 4D array (x,y,z,t) with PCs regressed out
%
%% Input validation

if ~ischar(cfg.sub)
    error('cfg.sub should contain a path to nifti file!')
end
if ~iscell(cfg.refs)
    error('cfg.refs should contain a cell array!')
end
if isempty(cfg.mask)
    error(['could not find mask: ' cfg.mask])
end

fprintf('\n')

%% Load mask

disp('Loading mask...')
mask=load_nii(cfg.mask);
inmask=find(mask.img);
fprintf('\n')

%% Load brain data
% Load subject NIFTI file
temp=load_nii(cfg.sub);
temp=permute(temp.img,[4 1 2 3]);
sub=zscore(temp(:,inmask));

% Get the number of time points, voxels, subjects-of-interest, and
% reference subjects
ntps=size(temp,1);
nvox=length(inmask);
nrefs=length(cfg.refs);

% Preallocate data matrices with the dimensions of time points x subjects x voxels 

allrefs=zeros(ntps,nrefs,nvox);

fprintf('\n')
disp('Loading brain data of reference subjects...')
for i=1:nrefs
    disp([cfg.refs{i} ' - ' num2str(i) ' out of ' num2str(nrefs)])
    
    temp=load_nii(cfg.refs{i});

    temp=permute(temp.img,[4 1 2 3]);
    allrefs(:,i,:)=zscore(temp(:,inmask));
end

fprintf('\n')

%% Regress out PCs of reference group from the subjects of interest
disp(['Calculating ' num2str(cfg.NPC) ' PCs from reference subjects data and regressing them out from subjects of interest data'])
dataReg=zeros(nvox,ntps);
for voxi=1:nvox
    if mod(voxi,1000)==0 % Show the status every 1000 voxels
        disp([num2str(voxi) '/' num2str(nvox) ' voxels'])
    end
    [~,PCs] = pca(allrefs(:,:,voxi)); % Get the principal comps of white matter and csf and regress those out
    PCs=PCs(:,1:cfg.NPC); % take the cfg.NPC number of components

    % Regress out PCs
    model=[ones(ntps,1) PCs];
    b=sub(:,voxi)'/model';
    TCreg=sub(:,voxi)-model*b';
    dataReg(voxi,:)=TCreg;
end

% Replace possible NaN values with zeros.
dataReg(isnan(dataReg))=0;

temp=zeros(ntps,size(mask.img,1),size(mask.img,2),size(mask.img,3));
temp(:,inmask)=dataReg';
cfg.dataReg=permute(temp,[2 3 4 1]);

disp('Done!')
fprintf('\n')