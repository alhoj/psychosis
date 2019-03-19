function cfg=calculate_indepRef_ISPS(cfg)
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

%% Input validation

if ~ismember(cfg.subs,{'ConsN30','PatsN36','BothN66'})
    error('cfg.subs should be either ''ConsN30'' or ''PatsN36''!')
end
if ~ismember(cfg.refs,{'ConsN15','ConsN30','PatsN36','BothN66'})
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
    error(['could not find mask: ' cfg.mask])
end
if ~ismember(cfg.useMeanOverRefs,[0 1])
    error('cfg.useMeanOverRefs must be either 0 or 1!')
end
if ~ismember(cfg.useNonSpatialSmoothedData,[0 1])
    error('cfg.useNonSpatialSmoothedData must be either 0 or 1!')
end
if ~ismember(cfg.averageOverTimePoints,[0 1])
    error('cfg.averageOverTimePoints must be either 0 or 1!')
end

disp(['Calculating ISPS for ' cfg.subs ' with respect to ' cfg.refs])
disp(['using ' cfg.condSubs ' data for ' cfg.subs ' and ' cfg.condRefs ' data for ' cfg.refs])
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
        cfg.subs=[consN30 patsN36];
end

switch cfg.refs
    case 'PatsN36'
        cfg.refs=patsN36;
    case 'ConsN30'
        cfg.refs=consN30;
    case 'ConsN15'
        cfg.refs=consN15;
    case 'BothN66'
        cfg.refs=[consN30 patsN36];
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

% Get the number of time points, voxels, subjects-of-interest, and
% reference subjects
notps=size(temp(:,inmask),1);
novox=size(temp(:,inmask),2);
nosubs=length(cfg.subs);
norefs=length(cfg.refs);

% Preallocate data matrices with the dimensions of time points x subjects x voxels 
allsubsTC=zeros(notps,nosubs,novox); 
allrefsTC=zeros(notps,norefs,novox);

% Filter frequencies
Fhp=.04; % high-pass
Flp=.07; % low-pass
TR=1.8;

disp('Loading brain data of subjects of interest...')
for i=1:nosubs
    disp([cfg.subs{i} ' - ' num2str(i) ' out of ' num2str(nosubs)])
    if cfg.useNonSpatialSmoothedData
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.condSubs '/' cfg.subs{i} '/epi_noSpatialSmoothing_' cfg.res '.nii']);
    else
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.condSubs '/' cfg.subs{i} '/epi_preprocessed_' cfg.res '.nii']);
    end
    temp=permute(nii.img,[4 1 2 3]);
    TCfull=zscore(temp(:,inmask));
    TCfilt=filterPMarray(TCfull,Fhp,Flp,TR);
    allsubsTC(:,i,:)=TCfilt;
end
fprintf('\n')
disp('Loading brain data of reference subjects...')
for i=1:norefs
    disp([cfg.refs{i} ' - ' num2str(i) ' out of ' num2str(norefs)])
    if cfg.useNonSpatialSmoothedData
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.condRefs '/' cfg.refs{i} '/epi_noSpatialSmoothing_' cfg.res '.nii']);
    else
        nii=load_nii(['/m/nbe/scratch/psykoosi/data/' cfg.condRefs '/' cfg.refs{i} '/epi_preprocessed_' cfg.res '.nii']);
    end
    temp=permute(nii.img,[4 1 2 3]);
    TCfull=zscore(temp(:,inmask));
    TCfilt=filterPMarray(TCfull,Fhp,Flp,TR);
    allrefsTC(:,i,:)=TCfilt;
end
fprintf('\n')

%% Calculate ISPS
disp('Calculating ISPS...')

% Hilbert transform
allsubsTCanal=hilbert(allsubsTC);
allrefsTCanal=hilbert(allrefsTC);
allsubsTCanal=permute(allsubsTCanal,[3 1 2]);
allrefsTCanal=permute(allrefsTCanal,[3 1 2]);
% allsubsTCanal=permute(allsubsTCanal,[1 3 2]);
% allrefsTCanal=permute(allrefsTCanal,[1 3 2]);

for s1=1:nosubs
    % Preallocate matrix
    angDif=zeros(novox,notps,norefs);
%     angDif=zeros(notps,novox,norefs);
    disp(['Subject ' cfg.subs{s1} ' - ' num2str(s1) ' out of ' num2str(nosubs)])
    for s2=1:norefs
        angDif(:,:,s2)=abs(angle(exp(1i*(angle(allsubsTCanal(:,:,s1))-angle(allrefsTCanal(:,:,s2)))))); % pairwise distances
    end

    % Average and save individual ISPS maps as nifti files
    % Check if the subject of interest is also a reference subject
    if ismember(cfg.subs{s1},cfg.refs)
        % Find out the indices of all other but the subject of interest in the reference subject list
        [~,~,inds] = setxor(cfg.subs{s1},cfg.refs);
        if cfg.averageOverTimePoints
            D=squeeze(angle(mean(mean(exp(1i*angDif(:,:,inds)),3),2))); % average across time and reference group
            isps=angle(exp(1i*(pi-2*D)))/pi;
            % Replace possible NaN values with zeros.
            isps(isnan(isps))=0;
            dirname=['/m/nbe/scratch/psykoosi/ISPS/indepRef_ISPS/' cfg.outdir '/tpsAveraged'];
            if ~exist(dirname,'dir') 
                system(['mkdir -p ' dirname]);
            end
            newbrain=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));
            newbrain(inmask)=isps;
            filename=[dirname '/' cfg.subs{s1} '.nii'];
            save_nii(make_nii(newbrain),filename);
            nii=fixOriginator(filename,mask);
            save_nii(nii,filename);
        else
            D=squeeze(angle(mean(exp(1i*angDif(:,:,inds)),3))); % average across reference group
            isps=angle(exp(1i*(pi-2*D)))/pi;
            % Replace possible NaN values with zeros.
            isps(isnan(isps))=0;
            for t=1:notps
                % Check if output directory exists; if not, create it
                dirname=['/m/nbe/scratch/psykoosi/ISPS/indepRef_ISPS/' cfg.outdir '/tp' num2str(t)];
                if ~exist(dirname,'dir') 
                    system(['mkdir -p ' dirname]);
                end
                newbrain=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));
                newbrain(inmask)=isps(:,t);
                filename=[dirname '/' cfg.subs{s1} '.nii'];
                save_nii(make_nii(newbrain),filename);
                nii=fixOriginator(filename,mask);
                save_nii(nii,filename);
            end
        end
    else
        if cfg.averageOverTimePoints
            D=squeeze(angle(mean(mean(exp(1i*angDif),3),2))); % average across time and reference group
            isps=angle(exp(1i*(pi-2*D)))/pi;
            % Replace possible NaN values with zeros.
            isps(isnan(isps))=0;
            dirname=['/m/nbe/scratch/psykoosi/ISPS/indepRef_ISPS/' cfg.outdir '/tpsAveraged'];
            if ~exist(dirname,'dir') 
                system(['mkdir -p ' dirname]);
            end
            newbrain=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));
            newbrain(inmask)=isps;
            filename=[dirname '/' cfg.subs{s1} '.nii'];
            save_nii(make_nii(newbrain),filename);
            nii=fixOriginator(filename,mask);
            save_nii(nii,filename);
        else
            D=squeeze(angle(mean(exp(1i*angDif),3))); % average across reference group
            isps=angle(exp(1i*(pi-2*D)))/pi;
            % Replace possible NaN values with zeros.
            isps(isnan(isps))=0;
            for t=1:notps
                % Check if output directory exists; if not, create it
                dirname=['/m/nbe/scratch/psykoosi/ISPS/indepRef_ISPS/' cfg.outdir '/tp' num2str(t)];
                if ~exist(dirname,'dir') 
                    system(['mkdir -p ' dirname]);
                end
                newbrain=zeros(size(mask.img,1),size(mask.img,2),size(mask.img,3));
                newbrain(inmask)=isps(:,t);
                filename=[dirname '/' cfg.subs{s1} '.nii'];
                save_nii(make_nii(newbrain),filename);
                nii=fixOriginator(filename,mask);
                save_nii(nii,filename);
            end
        end
    end
end

disp('Done!')
fprintf('\n')