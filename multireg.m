y = rand(245,1);
motionparam = load('epi_MCF.nii.par', '-ascii');

X = [ones(size(y)) motionparam];
[~,~,~,~,stats] = regress(y,X);
r2 = stats(1);