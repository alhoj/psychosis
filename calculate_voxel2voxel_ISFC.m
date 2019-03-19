%% Load fMRI timecourse data for ROIs
clear

addpath(genpath('/m/nbe/scratch/psykoosi/scripts'))

save_path='/m/nbe/scratch/psykoosi/ISFC/voxels_to_voxels/';


pats = {
    
%patients N=36, main sample
'/m/nbe/scratch/psykoosi/data/BL/EPHE602/'
'/m/nbe/scratch/psykoosi/data/BL/EPHE606/'
'/m/nbe/scratch/psykoosi/data/BL/EPHE607/'
'/m/nbe/scratch/psykoosi/data/BL/EPHE608/'
'/m/nbe/scratch/psykoosi/data/BL/EPHE610/'
'/m/nbe/scratch/psykoosi/data/BL/EPHE617/'
'/m/nbe/scratch/psykoosi/data/BL/EPHE618/'
'/m/nbe/scratch/psykoosi/data/BL/EPJO208/'
'/m/nbe/scratch/psykoosi/data/BL/EPJO209/'
'/m/nbe/scratch/psykoosi/data/BL/EPJO213/'
'/m/nbe/scratch/psykoosi/data/BL/EPJO214/'
'/m/nbe/scratch/psykoosi/data/BL/EPJO216/'
'/m/nbe/scratch/psykoosi/data/BL/EPJO217/'
'/m/nbe/scratch/psykoosi/data/BL/EPJO218/'
'/m/nbe/scratch/psykoosi/data/BL/EPJO220/'
'/m/nbe/scratch/psykoosi/data/BL/EPJO222/'
'/m/nbe/scratch/psykoosi/data/BL/EPJO224/'
'/m/nbe/scratch/psykoosi/data/BL/EPJO225/'
'/m/nbe/scratch/psykoosi/data/BL/EPJO227/'
'/m/nbe/scratch/psykoosi/data/BL/EPJO228/'
'/m/nbe/scratch/psykoosi/data/BL/EPPE113/'
'/m/nbe/scratch/psykoosi/data/BL/EPPE114/'
'/m/nbe/scratch/psykoosi/data/BL/EPPE115/'
'/m/nbe/scratch/psykoosi/data/BL/EPPE116/'
'/m/nbe/scratch/psykoosi/data/BL/EPPE117/'
'/m/nbe/scratch/psykoosi/data/BL/EPPE119/'
'/m/nbe/scratch/psykoosi/data/BL/EPPE121/'
'/m/nbe/scratch/psykoosi/data/BL/EPPE122/'
'/m/nbe/scratch/psykoosi/data/BL/EPPE128/'
'/m/nbe/scratch/psykoosi/data/BL/EPPE132/'
'/m/nbe/scratch/psykoosi/data/BL/EPPE138/'
'/m/nbe/scratch/psykoosi/data/BL/EPPE145/'
'/m/nbe/scratch/psykoosi/data/BL/EPPE149/'
'/m/nbe/scratch/psykoosi/data/BL/EPPE151/'
'/m/nbe/scratch/psykoosi/data/BL/EPPE154/'
'/m/nbe/scratch/psykoosi/data/BL/EPPK401/'

%patients N=28, not in FU
% '/m/nbe/scratch/psykoosi/data/BL/EPHE609/'
% '/m/nbe/scratch/psykoosi/data/BL/EPHE613/'
% '/m/nbe/scratch/psykoosi/data/BL/EPHE614/'
% '/m/nbe/scratch/psykoosi/data/BL/EPHE619/'
% '/m/nbe/scratch/psykoosi/data/BL/EPJO210/'
% '/m/nbe/scratch/psykoosi/data/BL/EPJO211/'
% '/m/nbe/scratch/psykoosi/data/BL/EPJO215/'
% '/m/nbe/scratch/psykoosi/data/BL/EPJO226/'
% '/m/nbe/scratch/psykoosi/data/BL/EPKK303/'
% '/m/nbe/scratch/psykoosi/data/BL/EPPE118/'
% '/m/nbe/scratch/psykoosi/data/BL/EPPE120/'
% '/m/nbe/scratch/psykoosi/data/BL/EPPE124/'
% '/m/nbe/scratch/psykoosi/data/BL/EPPE126/'
% '/m/nbe/scratch/psykoosi/data/BL/EPPE129/'
% '/m/nbe/scratch/psykoosi/data/BL/EPPE131/'
% '/m/nbe/scratch/psykoosi/data/BL/EPPE134/'
% '/m/nbe/scratch/psykoosi/data/BL/EPPE135/'
% '/m/nbe/scratch/psykoosi/data/BL/EPPE136/'
% '/m/nbe/scratch/psykoosi/data/BL/EPPE137/'
% '/m/nbe/scratch/psykoosi/data/BL/EPPE139/'
% '/m/nbe/scratch/psykoosi/data/BL/EPPE140/'
% '/m/nbe/scratch/psykoosi/data/BL/EPPE143/'
% '/m/nbe/scratch/psykoosi/data/BL/EPPE148/'
% '/m/nbe/scratch/psykoosi/data/BL/EPPE150/'
% '/m/nbe/scratch/psykoosi/data/BL/EPPE152/'
% '/m/nbe/scratch/psykoosi/data/BL/EPPE153/'
% '/m/nbe/scratch/psykoosi/data/BL/EPPE157/'
% '/m/nbe/scratch/psykoosi/data/BL/EPPE158/'
};

cons = {
%controls N=30, main sample
'/m/nbe/scratch/psykoosi/data/BL/EPVE516/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE517/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE520/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE524/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE525/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE526/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE527/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE530/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE531/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE533/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE536/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE540/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE541/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE542/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE543/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE544/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE546/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE548/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE549/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE551/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE554/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE555/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE556/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE557/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE558/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE559/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE560/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE561/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE562/'
'/m/nbe/scratch/psykoosi/data/BL/EPVE563/'

%controls N=15 (reference)
% '/m/nbe/scratch/psykoosi/data/BL/EPVE518/'
% '/m/nbe/scratch/psykoosi/data/BL/EPVE519/'
% '/m/nbe/scratch/psykoosi/data/BL/EPVE521/'
% '/m/nbe/scratch/psykoosi/data/BL/EPVE522/'
% '/m/nbe/scratch/psykoosi/data/BL/EPVE523/'
% '/m/nbe/scratch/psykoosi/data/BL/EPVE529/'
% '/m/nbe/scratch/psykoosi/data/BL/EPVE532/'
% '/m/nbe/scratch/psykoosi/data/BL/EPVE534/'
% '/m/nbe/scratch/psykoosi/data/BL/EPVE537/'
% '/m/nbe/scratch/psykoosi/data/BL/EPVE538/'
% '/m/nbe/scratch/psykoosi/data/BL/EPVE539/'
% '/m/nbe/scratch/psykoosi/data/BL/EPVE545/'
% '/m/nbe/scratch/psykoosi/data/BL/EPVE547/'
% '/m/nbe/scratch/psykoosi/data/BL/EPVE550/'
% '/m/nbe/scratch/psykoosi/data/BL/EPVE552/'
};

res='16mm';

% load mask
mask=load_nii(['/m/nbe/scratch/psykoosi/masks/MNI152_T1_' res '_brain_mask.nii']);
mask=logical(mask.img);

%% Read voxel timecourses
fprintf('Loading data:\n');

% for preallocation to speed things up
hdr=load_nii_hdr([pats{1} 'epi_preprocessed_' res '.nii']);
ntps=hdr.dime.dim(5);
nvox=length(find(mask));
npat=length(pats);
ncon=length(cons);

TCvoxels_pats=zeros(ntps,nvox,npat);
TCvoxels_cons=zeros(ntps,nvox,ncon);

for pati=1:npat
    disp(['Patient ' num2str(pati) '/' num2str(length(pats))]);
    
    nii=load_nii([pats{pati} 'epi_preprocessed_' res '.nii']);
    temp=reshape(nii.img,[],ntps);

    TCvoxels_pats(:,:,pati)=zscore(temp(mask,:)');

end
fprintf('\n');
for coni=1:ncon
    disp(['Control ' num2str(coni) '/' num2str(length(cons))]);
    
    nii=load_nii([cons{coni} 'epi_preprocessed_' res '.nii']);
    temp=reshape(nii.img,[],ntps);

    TCvoxels_cons(:,:,coni)=zscore(temp(mask,:)');

end
fprintf('\n');
% disp('saving...')
% save([save_path 'TCvoxels' res '_BL_patients'],'TCvoxels_pats', '-v7.3');
% save([save_path 'TCvoxels' res '_BL_controls'],'TCvoxels_cons', '-v7.3');
disp('done!')
%% Create correlation matrices

% note! ISFC is not symmetric
cMatPatient2Patient=zeros(nvox,nvox,npat*(npat-1));
cMatControl2Control=zeros(nvox,nvox,ncon*(ncon-1));

%
fprintf('Calculating patients ISFC:\n');
idx=0;
for s1=1:npat
    disp(num2str(s1));
    for s2=1:npat
        if s1~=s2
            idx=idx+1;
            disp(num2str(s2));
            cMatPatient2Patient(:,:,idx)=corr(TCvoxels_pats(:,:,s1),TCvoxels_pats(:,:,s2));
        end
    end
    fprintf('done\n');
end

%
fprintf('Calculating controls ISFC:\n');
idx=0;
for s1=1:ncon
    disp(num2str(s1));
    for s2=1:ncon
        if s1~=s2
            idx=idx+1;
            disp(num2str(s2));
            cMatControl2Control(:,:,idx)=corr(TCvoxels_cons(:,:,s1),TCvoxels_cons(:,:,s2));
        end
    end
    fprintf('done\n');
end

save([save_path 'ISFC_voxels_to_voxels'],'cMatPatient2Patient','cMatControl2Control', '-v7.3');
