%% regress covariates

cond = 'BL';

load([cond '_brainnetome_EPIgroupMask_connMat.mat'])
data = reshape(Mat,size(Mat,1)*size(Mat,2),size(Mat,3));
data = data';

load(['main_symptoms_' cond '.mat'])
% load('inds_remisN24.mat')
% scores = scores(inds_remisN24,:);

% design matrix for FD and cpz
X = [ones(size(Mat,3),1) scores{:,5:6}];

% regress out via pseudoinverse
% clean_data = data - X*pinv(X)*data;
% clean_data = reshape(clean_data',size(Mat,1),size(Mat,2),size(Mat,3));

% or alternatively
clean_data = nan(size(data));
for i=1:size(data,2)
    [b,~,clean_data(:,i)]=regress(data(:,i),X);
    % add the intercept to the residuals
    clean_data(:,i) = clean_data(:,i)+b(1);
end
clean_data = reshape(clean_data',size(Mat,1),size(Mat,2),size(Mat,3));

Mat = clean_data;
save([cond '_brainnetome_EPIgroupMask_connMat_FDcpzRegressed.mat'],'Mat','-v7.3')
% save('main_symptoms_BL_remisN24.mat','scores')
