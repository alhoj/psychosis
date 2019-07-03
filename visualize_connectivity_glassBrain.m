clear
close all

addpath('/m/nbe/scratch/braindata/shared/toolboxes/NIFTI');

% connectivity matrix
load('/m/nbe/scratch/psykoosi/ISFC/anova2x2interaction_EPIgroupMask_NBSresults_cdt001.mat')
conMat=nbs.NBS.con_mat{1}+nbs.NBS.con_mat{1}';
conMat=full(conMat);

% background template
bg=load_nii('/m/nbe/scratch/psykoosi/scripts/templates/MNI152_T1_2mm_brain.nii');

% mask and atlas
mask='/m/nbe/scratch/psykoosi/masks/EPIgroupMask_2mm.nii';
mask=load_nii(mask);
atlas='/m/nbe/scratch/psykoosi/masks/brainnetome_atlas_w_cerebellum_v2.nii';
atlas=load_nii(atlas);
assert(isequal(length(atlas.img(:)),length(mask.img(:))),'mask and atlas voxel resolutions do not match!')
atlas.img=double(atlas.img).*double(mask.img);
lastROI=max(atlas.img(:));
missing=setxor(1:lastROI,nonzeros(unique(atlas.img)));

for i=length(missing):-1:1
    atlas.img(atlas.img>missing(i))=atlas.img(atlas.img>missing(i))-1;
end

clusts=atlas.img;

thr=0;
filePrefix='/m/nbe/scratch/psykoosi/figures/brainnetome';

visualizeBrainNetwork(conMat,thr,atlas.img,bg,filePrefix)

