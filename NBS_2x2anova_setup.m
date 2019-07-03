clear

% let's model a 2x2 repeated measures
% data is ordered according to the following rows (2 groups, 2 time points,
% 29 subjects in the 1st group, 36 in the 2nd group)

ngroup1=29;
ngroup2=36;

ntime=2;

conmat=repmat(eye(ngroup1),ntime,1);
patmat=repmat(eye(ngroup2),ntime,1);

%% design matrix
% add main effect of time
% design=[ones(ngroup1,1); -1*ones(ngroup1,1); ones(ngroup2,1); -1*ones(ngroup2,1)];
% add group x time interaction
% design=[design [ones(ngroup1,1); -1*ones(ngroup1,1); -1*ones(ngroup2,1); ones(ngroup2,1)]];
design=[-1* ones(ngroup1,1); ones(ngroup1,1); ones(ngroup2,1); -1* ones(ngroup2,1)];
% add individual means of each subject
design=[design [conmat zeros(size(conmat,1),size(patmat,2)); zeros(size(patmat,1),size(conmat,2)) patmat]];
% design=[conmat zeros(size(conmat,1),size(patmat,2)); zeros(size(patmat,1),size(conmat,2)) patmat];
% design=[design [ones(ngroup1,1); -1*ones(ngroup1,1); -1*ones(ngroup2,1); ones(ngroup2,1)]];
%% contrasts
% contrast for the main effect of time 
main_effect_time=[1 zeros(1,size(design,2)-1)];
% contrast for the group x time interaction
interaction_groupXtime=[0 1 zeros(1,size(design,2)-2)];
% [1 zeros(1,65)]

%% the exhange block for NBS
exchange_block=[1:ngroup1 1:ngroup1 ngroup1+1:ngroup1+ngroup2 ngroup1+1:ngroup1+ngroup2]';
% [1:29 1:29 30:65 30:65]';
