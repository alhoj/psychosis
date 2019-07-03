% turns 3D-nifti into binary mask with ones assigned for n voxels with 
% the highest values 
clear all; close all;

% define input nifti filename 
infile='something.nii';

% how many highest values?
n=1000;

nii=niftiread(infile);
hdr=niftiinfo(infile);
[~,inds]=sort(nii(:),'descend');
mask=zeros(length(inds),1);
mask(inds(1:n))=1;
mask=reshape(mask,size(nii,1),size(nii,2),size(nii,3));

% define output mask filename
outfile='something_mask.nii';
niftiwrite(mask,outfile,hdr)
disp('done!')
