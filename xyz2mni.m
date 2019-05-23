function MNI = xyz2mni(xyz,ref)
% Converts matrix indices [x,y,z] to MNI coordinates
% ref = nifti file of data with the same (isotropic) voxel size
%
% Example: read MNI brain coordinates from 8mm standard mask
% ref='/m/nbe/scratch/psykoosi/masks/MNI152_T1_8mm_brain_mask.nii';
% nii=niftiread(ref);
% [x,y,z]=ind2sub(size(nii),find(nii));
% MNI = xyz2mni([x,y,z],ref)

xyz = double(xyz);
hdr = niftiinfo(ref);
res = hdr.PixelDimensions(1);

% transformation coefficients
T = [abs(hdr.Transform.T(4,1))/res abs(hdr.Transform.T(4,2))/res abs(hdr.Transform.T(4,3))/res] + [1 1 1];
T = repmat(T,[size(xyz,1),1]);

MNI = res*(xyz - T);

end

