
addpath(genpath('/m/nbe/scratch/psykoosi/scripts/'));

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

% Set IDs for subjects of interest
switch cfg.subs
    case 'patsN36'
        cfg.subs=patsN36;
    case 'consN29'
        cfg.subs=consN29;
    case 'consN15'
        cfg.subs=consN15;
end
%%
% INPUT DATA // An array of strings of valid files with location
for s=1:length(cfg.subs)
    cfg.indata{s}=[cfg.indir '/' cfg.subs{s} '.nii'];
end

% OUTPUT FOLDER // A folder where the software will write all data. 
cfg.outpath=cfg.outdir;   % Please, create the folder before using it

% SAMPLING RATE // In Hertz (1/TR)
cfg.Fs=1/1.8; 

% FILTER BAND // In Hertz, it has to be a four element vector
cfg.F=[0.025 0.04 0.07 0.09];

% FILTER DEVIATION // see 'help firpmord'
cfg.DEV=[0.05 0.01 0.05];

% BRAIN MASKS // masks that the software will use
cfg.coregistered_mask='funpsy/atlases/masks/MNI152_T1_2mm_brain_mask.nii';
cfg.compute_group_mask = 1;     % if = 1, it computes a group mask based on the power of each voxel
cfg.compute_spectrum = 0;       % if = 1, it computes a group frequency spectrum for each voxel

cfg.overwrite = 1; % overwrite if session already exists

% RUN PRE-ANALYSIS
sessionfile=funpsy_makepsess(cfg);  % validates input parameters
                                    % sessionfile will be the only variable needed in next function calls
                                    % it's a string with the path to the matlab file with all informations about
                                    % this phase analysis session


%% MAKE THE DATA // creates analytic signal

cfg=[];
cfg.sessionfile = sessionfile;
% cfg.compute_group_mask=1;     % overrides session settings
% cfg.compute_spectrum=1;       % overrides session settings
out = funpsy_makedata(cfg);     % filters, compute masks and creates the analytic signal

%% Pairwise ROIs analysis

cfg=[];
cfg.roimask='/m/nbe/scratch/psykoosi/masks/brainnetome_atlas_w_cerebellum_v2.nii';
atlas=load_nii(cfg.roimask);
mask=load_nii('/m/nbe/scratch/psykoosi/masks/MNI152_T1_2mm_brain_mask.nii');
rois=nonzeros(unique(double(atlas.img).*double(mask.img)));
for i=1:length(rois)
    cfg.labels{i}=rois(i);
end
cfg.overwrite = 1;
cfg.rois = bramila_makeRoiStruct(cfg);
cfg.sessionfile = sessionfile;
cfg.usemean = 0; 				% set usemean = 1 if you want to just do a mean of the voxels in the region. Default is usemean =0, which uses first principal component
out = funpsy_makeroidata(cfg); % extract roi time series based on the 1st principal component or the mean

% SBPS (Seed Based Phase Synchrony)
cfg=[];
cfg.sessionfile=sessionfile;
%cfg.useppc=1;					% Pairwise phase consistency is not implemented for SBPS since it needs testing
out = funpsy_sbps(cfg);         % takes a list of seeds/rois and computes full differential functional phase synchrony
                                % between each pair of seeds/rois.
                                % results stored in out.results.sbps

                        
%% Statistics
% SBPS (pairwise ROI analysis)
cfg = [];
cfg.sessionfile = sessionfile;

out = funpsy_avgsbps(cfg);    % calculate avgsbps

cfg.nonparam = 1;     % recommended. If 0 uses parametric tests. 
cfg.parallel = 1;   % Experimental feature - uses parallel computing
cfg.perm = 1000;     % for each ROI pair, does a non parametric test
cfg.blacklist = [];   % NOTE: the rois are already specified in data creation
                    % To modify them you need to rerun the analysis
cfg.statstype = 'sbps';
out = funpsy_stats(cfg);    % results in out.sbps_stats