clear

network = 'FEP_change'; % 'interaction' or 'FEP_change'

if strcmp(network, 'interaction')
    load('anova2x2interaction_EPIgroupMask_NBSresults_FDcpzNuisance_cdt001.mat');
else
    load('pairedSamples_patsN36_BL_vs_FU_EPIgroupMask_NBSresults_FDcpzNuisance_cdt001')
end
    
group = 'FEP'; % 'FEP', 'PC', 'all'
cond = 'diff'; % 'BL', 'FU', 'diff'
connections = 'significant'; % 'all'
corr_type = 'partial'; % 'partial', 'semipartial', 'basic'
score_labels = {'POS','NEG','Delusions','Hallucinations',...
                'bprs1011','bprs1215','Identification',...
                'Reality', 'dPrime'}; % 'anx','Depression','Suicidality',
covar_labels = {'FD', 'cpz'}; % 'FEP'
spearman = 0; % spearman or pearson
zscoring = 0;
nonparametric = 1; % use nonparametric Wilcoxon tests or parametric t-tests

% adjacency matrix
adj = nbs.NBS.con_mat{1}+nbs.NBS.con_mat{1}';
adj = full(adj);

% node/ROI labels
load('roi_labels.mat')
nbs.NBS.node_label = roi_labels;
n_roi = length(roi_labels);

% indices of connections
if strcmp(connections,'significant')
    [i,j] = find(nbs.NBS.con_mat{1}); 
    inds = find(nbs.NBS.con_mat{1});
elseif strcmp(connections,'all')
    inds = find(triu(ones(nbs.STATS.N),1)); % whole brain
end
n_conn = length(inds);

% setting up   
n_PC = 29;
n_FEP = 36;
PC_inds = 1:n_PC;
FEP_inds = n_PC+1:n_PC+n_FEP;
PC_color = [0 0.447 0.741];
FEP_color = [0.85 0.325 0.098];

% load brain data
brains_BL = load('BL_brainnetome_EPIgroupMask_connMat.mat'); 
brains_BL = brains_BL.Mat;
brains_BL = permute(brains_BL,[3 1 2]);
brains_FU = load('FU_brainnetome_EPIgroupMask_connMat.mat'); 
brains_FU = brains_FU.Mat;
brains_FU = permute(brains_FU,[3 1 2]);

% load behavioral data
behav_BL = load('main_symptoms_BL.mat');
behav_FU = load('main_symptoms_FU.mat');

n_sub = size(brains_BL,1);
conns_BL = brains_BL(:,inds); % BL: subjects x connections
conns_FU = brains_FU(:,inds); % FU: subjects x connections
conns_diff = conns_FU-conns_BL;

%
if strcmp(corr_type,'partial') || strcmp(corr_type,'semipartial')    
    % design matrix for covariates
    X_BL = [ones(n_sub,1) behav_BL.scores{:,covar_labels}];
    X_FU = [ones(n_sub,1) behav_FU.scores{:,covar_labels}];

    % regress out from brain data via pseudoinverse
    conns_BL = conns_BL - X_BL*pinv(X_BL)*conns_BL;
    conns_FU = conns_FU - X_FU*pinv(X_FU)*conns_FU;
%     conns_diff = conns_FU-conns_BL;
end

scores_BL = behav_BL.scores{:,score_labels};
scores_FU = behav_FU.scores{:,score_labels};
if strcmp(corr_type,'partial')
    % regress out from behavioral data via pseudoinverse
    scores_BL = scores_BL - X_BL*pinv(X_BL)*scores_BL;
    scores_FU = scores_FU - X_FU*pinv(X_FU)*scores_FU;
end
scores_diff = scores_FU-scores_BL;

switch group
    case 'all'
        sub_inds = 1:n_PC+n_FEP; % all
    case 'PC'
        sub_inds = PC_inds; % controls
    case 'FEP'
        sub_inds = FEP_inds; % patients
end

switch cond
    case 'BL'
        conns = conns_BL;
        conns_mean = mean(atanh(conns_BL),2);
        conns_median = median(conns_BL,2);
        scores = scores_BL;
    case 'FU'
        conns = conns_FU;
        conns_mean = mean(atanh(conns_FU),2);
        conns_median = median(conns_FU,2);
        scores = scores_FU;
    case 'diff'
        conns = conns_FU-conns_BL;
        conns_mean = mean(atanh(conns_FU),2)-mean(atanh(conns_BL),2);
        conns_median = median(conns_FU,2)-median(conns_BL,2);
        scores = scores_diff;
end

data_mean = [mean(conns_BL,2); mean(conns_FU,2)];
data_median = [median(conns_BL,2); median(conns_FU,2)];

%% connectivity matrices of negative and positive connections in each condition

labels = {'PC_BL', 'FEP_BL', 'PC_FU', 'FEP_FU'};
for i=1:length(labels)
    if contains(labels{i}, 'PC')
        this_inds = PC_inds;
    else
        this_inds = FEP_inds;
    end
    conmat = zeros(n_roi);
    if contains(labels{i}, 'BL')
        conmat(inds(find(mean(conns_BL(this_inds,:),1)<0))) = -1;
        conmat(inds(find(mean(conns_BL(this_inds,:),1)>0))) = 1;
    else
        conmat(inds(find(mean(conns_FU(this_inds,:),1)<0))) = -1;
        conmat(inds(find(mean(conns_FU(this_inds,:),1)>0))) = 1;
    end
    conmat = conmat + conmat';
    save(sprintf('%s_network_%s_posNeg.mat', network, labels{i}), 'conmat')
end

%% Violin plots

data = data_mean;
addpath('Violinplot-Matlab-master')
if zscoring
    data = zscore(data);
end
cats = [repmat({'BL-PC'},n_PC,1); repmat({'BL-FEP'},n_FEP,1); ...
        repmat({'FU-PC'},n_PC,1); repmat({'FU-FEP'},n_FEP,1)];
v = violinplot(data, cats, 'GroupOrder', {'BL-PC','BL-FEP','FU-PC','FU-FEP'});
v(1).ViolinColor = [0, 0.4470, 0.7410];
v(2).ViolinColor = [0.8500, 0.3250, 0.0980];
v(3).ViolinColor = [0, 0.4470, 0.7410];
v(4).ViolinColor = [0.8500, 0.3250, 0.0980];
print(gcf,'groupXcond_median.png','-dpng','-r500')

%% tests
data = data_mean;
% paired samples test
if nonparametric
    [p,~,stats] = signrank(data(PC_inds), data(n_sub+PC_inds));
    r = stats.zval/sqrt(2*n_PC); % wilcoxon effect size
    fprintf('\nWilcoxon signed rank test BL vs. FU in PC: z=%.2f, p=%.6f, r=%.2f', ...
            stats.zval, p, r)
else
    [~,p,~,stats] = ttest(data(PC_inds), data(n_sub+PC_inds));
    d = computeCohen_d(data(PC_inds), data(n_sub+PC_inds), 'paired');
    fprintf('\npaired t-test BL vs. FU in PC: t=%.2f, p=%.6f, d=%.2f', ...
            stats.tstat, p, d)
end

if nonparametric
    [p,~,stats] = signrank(data(FEP_inds), data(n_sub+FEP_inds));
    r = stats.zval/sqrt(2*n_FEP); % wilcoxon effect size
    fprintf('\nWilcoxon signed rank test BL vs. FU in FEP: z=%.2f, p=%.6f, r=%.2f', ...
            stats.zval, p, r)
else
    [~,p,~,stats] = ttest(data(FEP_inds), data(n_sub+FEP_inds));
    d = computeCohen_d(data(FEP_inds), data(n_sub+FEP_inds), 'paired');
    fprintf('\npaired t-test BL vs. FU in FEP: t=%.2f, p=%.6f, d=%.2f', ...
            stats.tstat, p, d)
end

% independent samples test
if nonparametric
    [p,~,stats] = ranksum(data(PC_inds), data(FEP_inds));
    r = stats.zval/sqrt(n_PC+n_FEP); % wilcoxon effect size
    fprintf('\nWilcoxon rank sum test PC vs. FEP in BL: z=%.2f, p=%.6f, r=%.2f', ...
            stats.zval, p, r)
else
    [~,p,~,stats] = ttest2(data(PC_inds), data(FEP_inds));
    d = computeCohen_d(data(PC_inds), data(FEP_inds), 'independent');
    fprintf('\ntwo-sample t-test PC vs. FEP in BL: t=%.2f, p=%.6f, d=%.2f', ...
            stats.tstat, p, d)
end

if nonparametric
    [p,~,stats] = ranksum(data(n_sub+PC_inds), data(n_sub+FEP_inds));
    r = stats.zval/sqrt(n_PC+n_FEP); % wilcoxon effect size
    fprintf('\nWilcoxon rank sum test PC vs. FEP in FU: z=%.2f, p=%.6f, r=%.2f\n', ...
            stats.zval, p, r)
else
    [~,p,~,stats] = ttest2(data(n_sub+PC_inds), data(n_sub+FEP_inds));
    d = computeCohen_d(data(n_sub+PC_inds), data(n_sub+FEP_inds), 'independent');
    fprintf('\ntwo-sample t-test PC vs. FEP in FU: t=%.2f, p=%.6f, d=%.2f\n', ...
            stats.tstat, p, d)
end

close all
%% calculate behav/symptom correlations for mean network connectivity

type = 'mean'; % 'mean' or 'median'

corr_label = sprintf('%s - %s', group, cond);
fprintf('%s:\n', corr_label)
for scorei=1:size(scores,2)
    score_label = strrep(score_labels{scorei},'_',' ');
    score = scores(sub_inds,scorei);
%     conn = mean(conns(sub_inds,:),2);

if strcmp(type, 'mean')
    conn = conns_mean(sub_inds);
else
    conn = conns_median(sub_inds);
end
    if isnumeric(score) && nnz(score)>0
        if spearman
            [r,p] = corr(conn, score,'Type','Spearman');
        else
            [r,p] = corr(conn, score);
        end         
    end
    disp([score_label ': r=' num2str(r) ', p=' num2str(p)])
    
    if p<0.1
        if zscoring
            conn = zscore(conn);
            score = zscore(score);
        end
        
        if strcmp(group,'all')
            scatter(score(PC_inds),conn(PC_inds),80,'filled',...
                    'MarkerFaceColor',PC_color,'MarkerEdgeColor','k')
            hold on
            scatter(score(FEP_inds),conn(FEP_inds),80,'filled',...
                    'MarkerFaceColor',FEP_color,'MarkerEdgeColor','k');
        else
            scatter(score,conn,80,'filled','MarkerFaceColor',FEP_color,...
                    'MarkerEdgeColor','k');
        end
        
        % plot regression line
        coef = polyfit(score, conn, 1);
        h = refline(coef(1), coef(2)); h.Color = 'k'; h.LineWidth = 2;        
        title([score_label ' - ' corr_label ', r=' num2str(r) ...
                ', p=' num2str(p)])
%         legend('PC', 'FEP')
        print(gcf,[score_label '.png'],'-dpng','-r500')
        close all
    end
end

%% calculate behav/symptom correlations for individual connections

fdr_correct = 1;

corr_label = sprintf('%s - %s', group, cond);
fprintf('%s:\n', corr_label)
for scorei=1:size(scores,2)
    score_label = strrep(score_labels{scorei},'_',' ');
    score = scores(sub_inds,scorei);    
    if isnumeric(score) && nnz(score)>0
        if spearman
            [r,p] = corr(conns(sub_inds,:), score,'Type','Spearman');
        else
            [r,p] = corr(conns(sub_inds,:), score);
        end         
    end
    if fdr_correct
        p = mafdr(p,'BHFDR',true);
    end
    signif = find(p<0.05);
    if ~isempty(signif)
        disp(score_label)
        this_conmat = zeros(n_roi);
        for s=1:length(signif)
            i_lab = nbs.NBS.node_label{i(signif(s))};
            j_lab = nbs.NBS.node_label{j(signif(s))};
            fprintf('%s to %s. r=%0.2f, p=%0.4f\n', i_lab, j_lab, r(signif(s)), p(signif(s)));
            this_conmat(i(signif(s)), j(signif(s))) = 1;
        end
        % symmetrize
        this_conmat = this_conmat+this_conmat';
        % save connectivity matrix of significant connections
        save([score_label '_indivConnCorrFDR_conmat.mat'], 'this_conmat');
    end
end

%% plot correlations for mean network connectivity for both groups

score_label = 'dPrime';

scorei = find(strcmp(score_labels,score_label));
score = scores(:, scorei);
conn = conns_mean;
if ~strcmp(corr_type, 'basic')
    score = zscore(score);
    conn = zscore(conn);
end
scatter(score(PC_inds),conn(PC_inds),80,'filled',...
        'MarkerFaceColor',PC_color,'MarkerEdgeColor','k');
hold on
scatter(score(FEP_inds),conn(FEP_inds),80,'filled',...
        'MarkerFaceColor',FEP_color,'MarkerEdgeColor','k');
    
coef_PC = polyfit(score(PC_inds), conn(PC_inds), 1);
h = refline(coef_PC(1), coef_PC(2)); h.Color = PC_color; h.LineWidth = 2;
coef_FEP = polyfit(score(FEP_inds), conn(FEP_inds), 1);
h = refline(coef_FEP(1), coef_FEP(2)); h.Color = FEP_color; h.LineWidth = 2;
xpad = (max(score)-min(score)) * 0.05;
ypad = (max(conn)-min(conn)) * 0.05;
axis([min(score)-xpad max(score)+xpad min(conn)-ypad max(conn)+ypad])
% axis([-2 102 -0.022 0.044])
h = get(gca,'Children');
set(gca,'Children',[h(4) h(3) h(2) h(1)])

title(sprintf('%s - %s', cond, score_label))
% legend('PC', 'FEP')
print(gcf,sprintf('%s_%s_%s.png', cond, score_label, corr_type),'-dpng','-r500')
close all