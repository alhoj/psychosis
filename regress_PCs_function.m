function cfg=regress_PCs_function(cfg)
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
if ~ismember(cfg.method,{'covmat','bold'})
    error('cfg.method must be either ''covmat'' or ''bold''!')
end

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

%% Load mask
addpath(genpath('/m/nbe/scratch/psykoosi/scripts'));

disp('Loading mask...')
if isempty(cfg.mask)
    mask=load_nii(['/m/nbe/scratch/psykoosi/masks/MNI152_T1_' cfg.res '_brain_mask.nii']);
else
    mask=load_nii(cfg.mask);
end
inmask=find(mask.img);
fprintf('\n')

%% Load brain data
% Load an example NIFTI file to find the dimensions needed for
% preallocating the matrices (this will speed things up)
test_nii=load_nii([cfg.indir '/' cfg.infile '.nii']);

% Put time first, so that we can acess the x y z with the 1-D indices
temp=permute(test_nii.img,[4 1 2 3]);

% Get the number of time points, voxels, subjects-of-interest, and
% reference subjects
ntps=size(temp,1);
nvox=length(inmask);
nref=length(cfg.refs);

% Preallocate data matrix with the dimensions of time points x subjects x voxels 
allrefs=zeros(ntps,nref,nvox);

disp('Loading brain data...')
disp(['Subject ' cfg.infile])
nii=load_nii([cfg.indir '/' cfg.infile '.nii']);
temp=permute(nii.img,[4 1 2 3]);
switch cfg.normalization
    case 'zscore'
        data=zscore(temp(:,inmask));
    case 'detrend'
        data=detrend(temp(:,inmask));
    case 'none'     
        data=temp(:,inmask);
end

for i=1:nref
    disp(['Reference ' cfg.refs{i} ' - ' num2str(i) ' out of ' num2str(nref)])    
%     nii=load_nii([cfg.indir '/' cfg.refs{i} '.nii']);
    nii=load_nii(['/m/nbe/scratch/psykoosi/data/BL_' cfg.res '_butterBandpass/' cfg.refs{i} '.nii']);
    temp=permute(nii.img,[4 1 2 3]);
    switch cfg.normalization
        case 'zscore'
            allrefs(:,i,:)=zscore(temp(:,inmask));
        case 'detrend'
            allrefs(:,i,:)=detrend(temp(:,inmask));
        case 'none'
            allrefs(:,i,:)=temp(:,inmask);
    end
end

fprintf('\n')

%% Regress out PCs of reference group from the subject of interest

if isequal(cfg.method,'covmat')

    disp(['Calculating ' num2str(cfg.NPC) ' PCs from references using covariance matrices'])

    % sum of covariance matrices of reference subjects
    covmat_refs_sum=zeros(ntps);
    for i=1:nref
        covmat_refs{i}=cov(squeeze(allrefs(:,i,:))');
        covmat_refs_sum=covmat_refs_sum+covmat_refs{i};
    end
    
    Urefs=[];Srefs=[];
    if isequal(cfg.comps,'intrinsic')
        disp(['Regressing intrinsic/individual components out from ' cfg.infile])
        
        if cfg.nonparametric
            disp('Using nonparametric estimation for determining common component count')
            
            % find group/shared components
            [Urefs,Srefs]=svd(covmat_refs_sum);
            Srefs = diag(Srefs);
            
            % find common component count using permutations for circularly shifted timeseries
            for j = 1:1000
%                 disp(j)
                covmat_refs_null=zeros(ntps);
                for i=1:nref
                    perm = circshift(1:ntps,[0,randi(ntps-1)]);
                    % permuting covariance matrix rows & cols equals to permuting original
                    % timeseries
%                     temp=cov(squeeze(allrefs(:,i,:))');
                    covmat_refs_null=covmat_refs_null+covmat_refs{i}(perm,perm);
                end
                [~,s]=svds(covmat_refs_null,1);
                nullvals(j)=s;
            end
            % retain components that surpass 1% of top null values
            Urefs=Urefs(:,find(Srefs>prctile(nullvals,1)));
            
            % design matrix of extrinsic processing i.e. shared components
            X=[ones(ntps,1) Urefs];
        
            % subject's covariance matrix
            covmat_sub=cov(data');
            % regress group components from subject's covariance matrix
            sub_intrinsic_comps = covmat_sub - X*pinv(X)*covmat_sub;
        else
            disp('Using parametric estimation for determining common component count (default in maxCorr)')
            addpath(genpath('/m/nbe/scratch/braindata/jaalho/toolboxes/maxCorr/'));
            % find common component count using parametric estimate (conservative)
            % create maxcorr-object
            mc = maxCorr(cat(3,data,permute(allrefs,[1 3 2])));
            w = -ones(1,mc.N);
            % set value 1 for the subject-of-interest
            w(1) = 1;
            % get the individual noise components
            [sub_intrinsic_comps,~] = mc.separate(w);
            sub_intrinsic_comps=double(sub_intrinsic_comps);
        end
        
         % regress out individual components from fMRI data via pseudoinverse
         % X=design matrix of intrinsic processing i.e. individual components
        X=[ones(ntps,1) sub_intrinsic_comps(:,1:cfg.NPC)];
        clean_data = data - X*pinv(X)*data;
        
    elseif isequal(cfg.comps,'extrinsic')  
        disp(['Regressing extrinsic/shared components out from ' cfg.infile])
        % find group/shared components
        [Urefs,~]=svd(covmat_refs_sum);
        % design matrix of extrinsic processing i.e. shared components
        X=[ones(ntps,1) Urefs(:,1:cfg.NPC)];
        % regress extrinsic components out from fMRI data via pseudoinverse
        clean_data = data - X*pinv(X)*data;
    end

    fprintf('\n')

elseif isequal(cfg.method,'bold')

    disp(['Calculating ' num2str(cfg.NPC) ' PCs from references using voxelwise BOLD signals'])
    clean_data=zeros(ntps,nvox);
    for voxi=1:nvox
        if mod(voxi,1000)==0 % Show status every 1000 voxels
            disp([num2str(voxi) '/' num2str(nvox) ' voxels'])
        end
        [~,PCs] = pca(allrefs(:,:,voxi)); % Get the PCs
        PCs=PCs(:,1:cfg.NPC); % Take cfg.NPC number of PCs

        disp(['Regressing the components out from ' cfg.infile])
        % design matrix of shared (extrinsic stimulus processing) components
        X=[ones(ntps,1) PCs];
        % betas
        b=data(:,voxi)'/X';
        % regress out components
        clean_data(:,voxi)=data(:,voxi)-X*b';
    end
    fprintf('\n')
end

%% Save files
% Replace possible NaN values with zeros.
clean_data(isnan(clean_data))=0;

% Check if output directory exists; if not, create it
dirname=cfg.outdir;
if ~exist(dirname,'dir') 
    system(['mkdir -p ' dirname]);
end

% Save cleaned niftis
disp(['Saving files to directory ' dirname '/'])

newbrain=zeros(ntps,size(mask.img,1),size(mask.img,2),size(mask.img,3));
newbrain(:,inmask)=clean_data;
newbrain=permute(newbrain,[2 3 4 1]);
filename=[dirname '/' cfg.outfile '.nii'];
save_nii(make_nii(newbrain),filename);
nii=fixOriginator(filename,mask);
save_nii(nii,filename);

disp('Done!')
fprintf('\n')
