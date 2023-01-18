clear

load('anova2x2interaction_EPIgroupMask_NBSresults_FDcpzNuisance_cdt001.mat');
% load('pairedSamples_patsN36_BL_vs_FU_EPIgroupMask_NBSresults_FDcpzNuisance_cdt001')

group = 'PC'; % 'FEP', 'PC', 'all'
cond = 'diff'; % 'BL', 'FU', 'diff'
connections = 'signif'; % 'all'
corr_type = 'basic'; % 'semipartial', 'basic'
score_labels = {'POS','NEG','Delusions','Hallucinations',...
                'bprs1011','bprs1215','Identification','Reality'}; % 'anx','Depression','Suicidality',
n_scores = length(score_labels);
covar_labels = {'FD', 'cpz'}; % 'FEP'

% adjacency matrix
adj = nbs.NBS.con_mat{1}+nbs.NBS.con_mat{1}';
adj = full(adj);

% node/ROI labels
load('roi_labels.mat')
nbs.NBS.node_label = roi_labels;
n_rois = length(roi_labels);

% indices of connections
if strcmp(connections,'signif')
    [i,j] = find(nbs.NBS.con_mat{1}); 
    inds = find(nbs.NBS.con_mat{1});
elseif strcmp(connections,'all')
    inds = find(triu(ones(nbs.STATS.N),1)); % whole brain
end

if length(inds)==65
    network_id = 'interactionNetwork';
else
    network_id = 'FEPchangeNetwork';
end
out_label = sprintf('twISFC_%s_%s_%s_%sConns_%sCorr', group, cond, ...
                    network_id, connections, corr_type);
%% setting up & loading data  
n_PC = 29;
n_FEP = 36;
n_subs = n_PC+n_FEP;
PC_inds = 1:n_PC;
FEP_inds = n_PC+1:n_subs;

% load brain data
dataID = '10TRwin1TRstep_refconsN15_brainnetome_EPIgroupMask_2mm';
brains_PC_BL = load(['BL_consN' num2str(n_PC) '_' dataID '_corrs.mat']);
brains_FEP_BL = load(['BL_patsN' num2str(n_FEP) '_' dataID '_corrs.mat']);
brains_PC_FU = load(['FU_consN' num2str(n_PC) '_' dataID '_corrs.mat']);
brains_FEP_FU = load(['FU_patsN' num2str(n_FEP) '_' dataID '_corrs.mat']);

% load behavioral data
behav_BL = load('main_symptoms_BL.mat');
behav_FU = load('main_symptoms_FU.mat');

%% process data in time windows

scores_BL = behav_BL.scores{:,score_labels};
scores_FU = behav_FU.scores{:,score_labels};

% design matrix for covariates
X_BL = [ones(n_subs,1) behav_BL.scores{:,covar_labels}];
X_FU = [ones(n_subs,1) behav_FU.scores{:,covar_labels}];
if strcmp(corr_type,'partial')
    % regress out from behavioral data via pseudoinverse
    fprintf('\nRegressing out %s\n', strjoin(string(covar_labels),','))
    scores_BL = scores_BL - X_BL*pinv(X_BL)*scores_BL;
    scores_FU = scores_FU - X_FU*pinv(X_FU)*scores_FU;
end
scores_diff = scores_FU-scores_BL;

switch group
    case 'all'
        sub_inds = 1:n_subs; % all
    case 'PC'
        sub_inds = PC_inds; % controls
    case 'FEP'
        sub_inds = FEP_inds; % patients
end

n_tws = length(brains_PC_BL.results); % number of time windows
r = zeros(n_scores, n_tws);
p = zeros(n_scores, n_tws);
conns_mean = zeros(n_subs, n_tws);
for twi=1:n_tws
    disp(['time window ' num2str(twi)])
    % create connectivity matrices
    conmat_BL = zeros(n_subs,n_rois,n_rois);
    conmat_FU = zeros(n_subs,n_rois,n_rois);
    triu_inds = find(triu(ones(n_rois),1));

    conmat_BL(PC_inds,triu_inds) = brains_PC_BL.results(twi).corrs';
    conmat_BL(FEP_inds,triu_inds) = brains_FEP_BL.results(twi).corrs';
    conmat_BL = permute(conmat_BL,[2 3 1]);
    conmat_BL = conmat_BL+permute(conmat_BL,[2 1 3])+repmat(eye(n_rois),[1,1,n_subs]);
    conmat_BL = permute(conmat_BL,[3 1 2]);
    
    conmat_FU(PC_inds,triu_inds) = brains_PC_FU.results(twi).corrs';
    conmat_FU(FEP_inds,triu_inds) = brains_FEP_FU.results(twi).corrs';
    conmat_FU = permute(conmat_FU,[2 3 1]);
    conmat_FU = conmat_FU+permute(conmat_FU,[2 1 3])+repmat(eye(n_rois),[1,1,n_subs]);
    conmat_FU = permute(conmat_FU,[3 1 2]);
    
    % ISFC in the network
    conns_BL = conmat_BL(:,inds); % BL: subjects x connections
    conns_FU = conmat_FU(:,inds); % FU: subjects x connections
    
    if strcmp(corr_type,'partial') || strcmp(corr_type,'semipartial')    
        % regress out from brain data via pseudoinverse
        conns_BL = conns_BL - X_BL*pinv(X_BL)*conns_BL;
        conns_FU = conns_FU - X_FU*pinv(X_FU)*conns_FU;
    end
    conns_diff = conns_FU-conns_BL; % BL-FU difference
    
    %   
    switch cond
        case 'BL'
            conns = conns_BL;
            scores = scores_BL;
        case 'FU'
            conns = conns_FU;
            scores = scores_FU;
        case 'diff'
            conns = conns_diff;
            scores = scores_diff;
    end

    % calculate average connectivity in the network
    conns_mean(:,twi) = mean(conns,2);
    
    % calculate correlations
    for scorei=1:size(scores,2)
        score_label = score_labels{scorei};
        score = scores(:,scorei);
        if isnumeric(score) && nnz(score)>0
            [r(scorei,twi),p(scorei,twi)] = corr(mean(conns(sub_inds,:),2), score(sub_inds));
        end
    end
end
 
out.r = r;
out.p = p;
out.score_labels = score_labels;
out.covar_labels = covar_labels;
out.roi_labels = roi_labels;
out.cond = cond;
out.group = group;
out.connections = connections;
out.corr_type = corr_type;
out.sub_inds = sub_inds;
out.conn_inds = inds;
out.n_tws = n_tws;
out.conns_mean = conns_mean;

save([out_label '.mat'], 'out')

%% correlate twISFC with ratings
clear
load('twISFC_symptomCorr_FEP_diff_signifConns_partialCorr.mat')
load('ratings.mat')
conns_mean = out.conns_mean;

for rati=1:length(labels)
    rating = ratings_convHRF_zscored(6:end-5,rati);
    r_twISFC = corr(mean(conns_mean,1)',rating);
    DF = bramila_autocorr(mean(conns_mean,1)',rating);
    p_twISFC = pvalPearson('b', r_twISFC, DF+2);
    if p_twISFC < 0.05
        rating_label = labels{rati};
        disp(['twISFC vs. ' labels{rati} ': r=' num2str(r_twISFC) ', p=' num2str(p_twISFC)])
        fig = figure();
        plot(rating,'LineWidth',2)
        hold on
        plot(zscore(mean(conns_mean,1)),'LineWidth',2)
        plot_label = sprintf('%s_%s_%s_%s', group, cond, rating_label);
        title([strrep(plot_label,'_',' - ') ' r=' num2str(rR) ',p=' num2str(pR)])
        % Enlarge figure to full screen.
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0,0,1,1]);
        print(fig, ['twISFC_rating_corr_' plot_label] ,'-r500','-dpng');
        close all
    end  
end
%% correlate twISFC vs symptoms correlation time series with ratings
group = 'all';
load(['twISFC_symptomCorr_' group '_diff_signifConns_partialCorr.mat'])
load('ratings.mat')
conns_mean = out.conns_mean;
r = out.r;
score_labels = out.score_labels;
group = out.group;
if group=='all'
    group_color = [0.694 0.184 0.756]; % [0.494 0.184 0.556]
else
    group_color = [0.85 0.325 0.098];
end
cond = out.cond;
lw = 3; % line width
rating_color = [0.25 0.25 0.25]; % [0.629 0.694 0.125] [0.466 0.674 0.188]

for sympti=1:size(r,1)
    for rati=1:length(labels)
        rating = ratings_convHRF_zscored(6:end-5,rati);
%         rating=ratings_zscored(6:end-5,rati);
        rR = corr(r(sympti,:)',rating);
        DF = bramila_autocorr(r(sympti,:)',rating);
        pR = pvalPearson('b', rR, DF+2);
        if pR < 0.05
            score_label = strrep(score_labels{sympti},'_',' - ');
            rating_label = labels{rati};
            disp(['symptom: ' score_label ', rating: ' labels{rati} ', r=' num2str(rR) ', p=' num2str(pR)])
            fig = figure();
            plot(rating,'color',rating_color,'LineWidth',lw)
            hold on
            plot(zscore(r(sympti,:)),'color',group_color,'LineWidth',lw)
            plot_label = sprintf('%s_%s_%s_%s', group, cond, score_label, rating_label);
            title([strrep(plot_label,'_',' - ') ' r=' num2str(rR) ',p=' num2str(pR)])
            % Enlarge figure to full screen.
%             set(gcf, 'Units', 'Normalized', 'OuterPosition', [0,0,1,1]);
            set(gcf,'position',[0,0,800,300])
            print(fig, ['twISFC_symptom_rating_corr_' plot_label] ,'-r500','-dpng');
            close all
        end
    end
end