clear

n_FEP = 36;
n_PC = 29;

% load('anova2x2interaction_EPIgroupMask_NBSresults_FDcpzNuisance_cdt001.mat');
% load('/m/nbe/scratch/psykoosi/ISFC/anova2x2interaction_EPIgroupMask_FDRresults_p05_50000perms.mat')
% load('pairedSamples_patsN36_BL_vs_FU_EPIgroupMask_NBSresults_cdt001.mat')
% load('pairedSamples_patsN36_BL_vs_FU_EPIgroupMask_NBSresults_FDcpzNuisance_cdt001')

% adjacency matrix
n_roi = length(nbs.NBS.con_mat{1});
adj=nbs.NBS.con_mat{1}+nbs.NBS.con_mat{1}';
adj=full(adj);

fprintf('\n')
load('roi_labels.mat')
nbs.NBS.node_label=roi_labels;
[i,j]=find(nbs.NBS.con_mat{1}); 
for roii=1:n_roi
    iter = 0;
    stat_sum = 0;
    if roii<100
        roi = [' ',num2str(roii)];
    else
        roi = num2str(roii);
    end
    for n=1:length(i)
        i_lab=nbs.NBS.node_label{i(n)};
        j_lab=nbs.NBS.node_label{j(n)};
        stat=nbs.NBS.test_stat(i(n),j(n));

        if isequal(i_lab,roi) || isequal(j_lab,roi)
%             fprintf('%s to %s. Test stat: %0.2f\n',i_lab,j_lab,stat);
            iter=iter+1;
            inds_roi(iter)=sub2ind(size(adj),i(n),j(n));
            stat_sum = stat_sum + stat;
        end

    %     if stat>10        
    %         fprintf('%s to %s. Test stat: %0.2f\n',i_lab,j_lab,stat);
    %         iter=iter+1;
    %         inds_th(iter)=sub2ind(size(adj),i(n),j(n));
    %     end
    end
    if stat_sum
        fprintf('Region %s: %d connections; sum of stats=%0.1f\n',roi,iter,stat_sum);
    end
    inds=find(nbs.NBS.con_mat{1});
end
% inds=find(triu(ones(272),1)); % whole brain

%% print stats of a ROI

roi='218';
th = 1; % print connection exceeding this statistic threshold
for n=1:length(i)
    i_lab=nbs.NBS.node_label{i(n)};
    j_lab=nbs.NBS.node_label{j(n)};
    stat=nbs.NBS.test_stat(i(n),j(n));

    if isequal(i_lab,roi) || isequal(j_lab,roi)
            fprintf('%s to %s. Test stat: %0.2f\n',i_lab,j_lab,stat);
    end

%     if stat>th
%         fprintf('%s to %s. Test stat: %0.2f\n',i_lab,j_lab,stat);
%         iter=iter+1;
%         inds_th(iter)=sub2ind(size(adj),i(n),j(n));
%     end
end


%% symptom correlation
% addpath(genpath('/m/nbe/scratch/braindata/shared/toolboxes/bramila/bramila/')) 

cond = 'both'; % 'BL', 'FU', 'diff', or 'both'
group = 'FEP'; % 'FEP', 'PC', 'all', or 'remis'
connections = 'ROI'; % 'signif', 'ROI', 'th', or 'all'
partialCorr = 'semi'; % 'full' -> partial, 'semi' semipartial with FD and cpz as covars

if strcmp(cond,'BL')
    if strcmp(partialCorr,'semi')
        load('BL_brainnetome_EPIgroupMask_connMat_FDcpzRegressed.mat')
    else
        load('BL_brainnetome_EPIgroupMask_connMat.mat')
    end
    Mat=permute(Mat,[3 1 2]);
    load('main_symptoms_BL.mat')
elseif strcmp(cond,'FU')
    if strcmp(partialCorr,'semi')
        load('FU_brainnetome_EPIgroupMask_connMat_FDcpzRegressed.mat')
    else
        load('FU_brainnetome_EPIgroupMask_connMat.mat')
    end
    Mat=permute(Mat,[3 1 2]);
    load('main_symptoms_FU.mat')
elseif strcmp(cond,'diff')
    if strcmp(partialCorr,'semi')
        BL = load('BL_brainnetome_EPIgroupMask_connMat_FDcpzRegressed.mat');
        FU = load('FU_brainnetome_EPIgroupMask_connMat_FDcpzRegressed.mat');
        Mat = BL.Mat-FU.Mat;
    else
        load('diffBLFU_brainnetome_EPIgroupMask_connMat.mat')
    end
    Mat=permute(Mat,[3 1 2]);
    load('main_symptoms_diffBLFU.mat')
elseif strcmp(cond,'both')
    if strcmp(partialCorr,'semi')
        Mat_BL = load('BL_brainnetome_EPIgroupMask_connMat_FDcpzRegressed.mat'); 
        Mat_BL = Mat_BL.Mat;
        Mat_FU = load('FU_brainnetome_EPIgroupMask_connMat_FDcpzRegressed.mat');
        Mat_FU = Mat_FU.Mat;
    else
        Mat_BL = load('BL_brainnetome_EPIgroupMask_connMat.mat'); Mat_BL=Mat_BL.Mat;
        Mat_FU = load('FU_brainnetome_EPIgroupMask_connMat.mat'); Mat_FU=Mat_FU.Mat;
    end
    Mat_BL=permute(Mat_BL,[3 1 2]);
    Mat_FU=permute(Mat_FU,[3 1 2]);
    scores_BL = load('main_symptoms_BL.mat'); scores_BL=scores_BL.scores;
    scores_FU = load('main_symptoms_FU.mat'); scores_FU=scores_FU.scores;
    scores_diff = load('main_symptoms_diffBLFU.mat'); scores_diff=scores_diff.scores;
    scores = scores_BL;
end

if strcmp(group,'FEP')
    subs=n_PC+1:n_PC+n_FEP; % patients
    n = n_FEP;
elseif strcmp(group,'PC')
    subs=1:n_PC; % controls
    n = n_PC;
elseif strcmp(group,'all')
    subs=1:n_PC+n_FEP; % all
elseif strcmp(group,'remis')
    load('inds_remisN24.mat')
    subs=inds_remisN24;
end

% load('diffBLFU_remisN24_brainnetome_EPIgroupMask_connMat.mat')
% load('/m/nbe/scratch/psykoosi/ISFC/diffBLFU_patsN36_brainnetome_EPIgroupMask_connMat.mat')
% load('/m/nbe/scratch/psykoosi/ISFC/BL_patsN54_brainnetome_EPIgroupMask_connMat.mat')
% load('inds_remisN24.mat')
% scores=scores(inds_remisN24,:);
% load('/m/nbe/scratch/psykoosi/clinical_data_BL_patsN54.mat')
% load('/m/nbe/scratch/psykoosi/clinical_data_diffBLFU_patsN36.mat')

if strcmp(cond,'both')
    if strcmp(connections,'ROI')
        conns_BL=Mat_BL(subs,inds_roi);
        conns_FU=Mat_FU(subs,inds_roi);
    elseif strcmp(connections,'th')
        conns_BL=Mat_BL(subs,inds_th);
        conns_FU=Mat_FU(subs,inds_th);
    elseif strcmp(connections,'signif')
        conns_BL=Mat_BL(subs,inds);
        conns_FU=Mat_FU(subs,inds);
    elseif strcmp(connections,'all')
        conns_BL=Mat_BL(subs,:);
        conns_FU=Mat_FU(subs,:);
    end
else
    if strcmp(connections,'ROI')
        conns=Mat(subs,inds_roi);
    elseif strcmp(connections,'th')
        conns=Mat(subs,inds_th);
    elseif strcmp(connections,'signif')
        conns=Mat(subs,inds);
    elseif strcmp(connections,'all')
        conns=Mat(subs,:);
    end
end

% adjust for group difference
% conns(1:29,:)=detrend(conns(1:29,:));
% conns(30:65,:)=detrend(conns(30:65,:));

% correlations separately for each significant connection
for scorei=3:length(scores.Properties.VariableNames)
    if isnumeric(scores{:,scorei}) && nnz(scores{:,scorei})>0
%         [r,p]=corr(conns,scores{subs,scorei});
%         fdrP = mafdr(p,'BHFDR',true);
%         if any(fdrP<0.05)
% %             disp(['significant connection -> ' scores.Properties.VariableNames{scorei}])
%             significants=find(fdrP<0.05);
%             for n=significants'
%                 i_lab=nbs.NBS.node_label{i(n)};
%                 j_lab=nbs.NBS.node_label{j(n)};
%                 stat=nbs.NBS.test_stat(i(n),j(n));
% %                 fprintf('%s to %s. Test stat: %0.2f\n',i_lab,j_lab,stat);
%             end        
%         end
        % global correlation, i.e. with the average connectivity in the network
        if strcmp(group,'all') && ~strcmp(cond,'both')
            subsCons=1:n_PC;
            subsPats=n_PC+1:size(Mat,1);
            connsPats=Mat(subsPats,inds);
            connsCons=Mat(subsCons,inds);
                       
            if strcmp(partialCorr,'full')
                [rGpats,pGpats]=partialcorr(mean(connsPats,2),scores{subsPats,scorei},scores{subsPats,5:6});
                [rGcons,pGcons]=partialcorr(mean(connsCons,2),scores{subsCons,scorei},scores{subsCons,5:6});
                [rG,pG]=partialcorr(mean(conns,2),scores{subs,scorei},scores{subs,5:6});
            else
                [rGpats,pGpats]=corr(mean(connsPats,2),scores{subsPats,scorei});
                [rGcons,pGcons]=corr(mean(connsCons,2),scores{subsCons,scorei});
                [rG,pG]=corr(mean(conns,2),scores{subs,scorei});
            end
            [PcompR,z] = corr_rtest(rGpats,rGcons,n_FEP,n_PC);
            
            fprintf('\n')
            disp(['global correlation -> ' scores.Properties.VariableNames{scorei} ': r=' num2str(rG) ', p=' num2str(pG)])
            disp(['global correlation patients: r=' num2str(rGpats) ', p=' num2str(pGpats)])
            disp(['global correlation controls: r=' num2str(rGcons) ', p=' num2str(pGcons)])
            disp(['difference in correlation coefficients between patients and controls: z=' num2str(z) ', p=' num2str(PcompR(2))])

        elseif ~strcmp(group,'all') && strcmp(cond,'both')
            if strcmp(partialCorr,'full')
                [rG_BL,pG_BL]=partialcorr(mean(conns_BL,2),scores_BL{subs,scorei},scores_BL{subs,5:6});
                [rG_FU,pG_FU]=partialcorr(mean(conns_FU,2),scores_FU{subs,scorei},scores_FU{subs,5:6});
                [rG_diff,pG_diff]=partialcorr(mean(conns_BL,2)-mean(conns_FU,2),scores_diff{subs,scorei},scores_diff{subs,5:6});
            else
                [rG_BL,pG_BL]=corr(mean(conns_BL,2),scores_BL{subs,scorei});
                [rG_FU,pG_FU]=corr(mean(conns_FU,2),scores_FU{subs,scorei});
                [rG_diff,pG_diff]=corr(mean(conns_BL,2)-mean(conns_FU,2),scores_diff{subs,scorei});
            end
            z_BL = atanh(rG_BL);
            z_FU = atanh(rG_FU);
            z = (z_BL-z_FU)/sqrt(1/(n-3));
            PcompR = (1-normcdf(abs(z),0,1))*2;

            disp(['global correlation ' group ' -> ' scores.Properties.VariableNames{scorei}])
            disp(['global correlation in BL: r=' num2str(rG_BL) ', p=' num2str(pG_BL)])
            disp(['global correlation in FU: r=' num2str(rG_FU) ', p=' num2str(pG_FU)])
            disp(['global correlation in diff BL-FU: r=' num2str(rG_diff) ', p=' num2str(pG_diff)])
            disp(['difference in correlation coefficients between BL and FU: z=' num2str(z) ', p=' num2str(PcompR)])

        elseif ~strcmp(group,'all') && strcmp(cond,'diff')
            if strcmp(partialCorr,'full')
                [rG,pG]=partialcorr(mean(conns,2),scores{subs,scorei},scores_FU{subs,5:6});
            else
                [rG,pG]=corr(mean(conns,2),scores{subs,scorei});
            end
            
            disp(['global correlation ' group ' -> ' scores.Properties.VariableNames{scorei}])
            disp(['global correlation diff BL-FU: r=' num2str(rG) ', p=' num2str(pG)])
        end
%         DF=bramila_autocorr(scores{subs,scorei},scores{subs,scorei});
%         DF=bramila_autocorr(mean(conns,2),scores{subs,scorei});
%         pGadjusted=pvalPearson('b', rG, DF+2);     

    end
end

%% memory score correlation 

load('BL_brainnetome_EPIgroupMask_connMat.mat')
% load('BL_patsN77_brainnetome_EPIgroupMask_connMat')
Mat=permute(Mat,[3 1 2]);
load('memdata')

memScores = bl_WMSLOG_consN29_patsN36;
memScoresCons = memScores(1:29);
memScoresPats = memScores(30:end);
% memScores = bl_WMSLOG_patsN46(indsN46);

inds_memScoresSubs=find(memScores);
inds_memScoresCons=find(memScoresCons);
inds_memScoresPats=29+find(memScoresPats);

conns = Mat(inds_memScoresSubs,inds);
connsPats = Mat(inds_memScoresPats,inds);
connsCons = Mat(inds_memScoresCons,inds);
% conns=Mat(indsN77,inds);
% conns=Mat(indsN77,inds_roi);
memScores = memScores(find(memScores));
memScoresCons = memScoresCons(find(memScoresCons));
memScoresPats = memScoresPats(find(memScoresPats));

[r,p] = corr(conns,memScores);
fdrP = mafdr(p,'BHFDR',true);
if any(fdrP<0.05)
    disp('significant correlation found')
    significants=find(fdrP<0.05);
    for n=significants'
        i_lab=nbs.NBS.node_label{i(n)};
        j_lab=nbs.NBS.node_label{j(n)};
        stat=nbs.NBS.test_stat(i(n),j(n));
        fprintf('%s to %s. Test stat: %0.2f\n',i_lab,j_lab,stat);
    end        
end

% comparison between group-specific correlation
[r_pats,p_pats]=corr(connsPats,memScoresPats);
[r_cons,p_cons]=corr(connsCons,memScoresCons);
% fdr_p_pats = mafdr(p_pats,'BHFDR',true);
% fdr_p_cons = mafdr(p_cons,'BHFDR',true);
% p_compR = compare_correlation_coefficients(r_pats,r_cons,27,29);
for i=1:size(conns,2)
    [corr_rtest_results(i).p, corr_rtest_results(i).z, corr_rtest_results(i).z_pats, corr_rtest_results(i).z_cons] = corr_rtest(r_pats(i),r_cons(i),size(connsPats,1),size(connsCons,1));
end
[rG_pats,pG_pats] = corr(mean(connsPats,2),memScoresPats);
[rG_cons,pG_cons] = corr(mean(connsCons,2),memScoresCons);
% [rG_pats,pG_pats]=corr(mean(atanh(connsPats),2),memScoresPats);
% [rG_cons,pG_cons]=corr(mean(atanh(connsCons),2),memScoresCons);
% p_compRG = compare_correlation_coefficients(rG_pats,rG_cons,27,29);
[p_compRG, z_compRG, zG_pats, zG_cons] = corr_rtest(rG_pats,rG_cons,size(connsPats,1),size(connsCons,1));

% global correlation, i.e. with the average connectivity in the network
[rG,pG]=corr(mean(conns,2),memScores);
DF=bramila_autocorr(mean(conns,2),memScores);
pGadjusted=pvalPearson('b', rG, DF+2);     
disp(['global correlation -> r=' num2str(rG) ', p=' num2str(pGadjusted)])
% disp(['global correlation patients: r=' num2str(rGpats) ', p=' num2str(pGpats)])
% disp(['global correlation controls: r=' num2str(rGcons) ', p=' num2str(pGcons)])
% disp(['difference in correlation coefficients between patients and controls p=' num2str(PcompR)])
fprintf('\n')
%         if length(subs)==65
%             [rP,pP]=partialcorr(mean(conns,2),scores{subs,scorei},[ones(29,1); 2*ones(36,1)]);
%         end
%         [rP,pP]=partialcorr(mean(conns,2),scores{subs,scorei},scores{subs,18});
%         [rP,pP]=partialcorr(mean(conns,2),scores{subs,scorei},scores{subs,5});
%         pPadjusted=pvalPearson('b', rP, DF+2); 
%         disp(['partial correlation -> ' scores.Properties.VariableNames{scorei} ': r=' num2str(rP) ', p=' num2str(pPadjusted)])


%% correlation between symptom/behavioral scores

for i=3:size(scores,2)
    disp(scores.Properties.VariableNames{i})
    [r,p]=corr(table2array(scores(subs,22)),double(table2array(scores(subs,i))))
end
%% scatter plot symptom correlation

load('diffBLFU_brainnetome_EPIgroupMask_connMat.mat')
Mat=permute(Mat,[3 1 2]);
load('main_symptoms_diffBLFU.mat')
score='identification';
scorei=find(strcmp(scores.Properties.VariableNames,score));

load('inds_subgroups.mat')
subs=1:size(Mat,1);
conns=Mat(subs,inds);
scatter(mean(conns,2),scores{subs,scorei},[],'k','o')
hold on

conns=Mat(1:29,inds);
scatter(mean(conns,2),scores{1:29,scorei},[],'b','o')

conns=Mat(30:65,inds);
scatter(mean(conns,2),scores{30:65,scorei},[],'r','o')

% subs=inds_FEPnotRemis;
% conns=Mat(subs,inds);
% scatter(mean(conns,2),scores{subs,scorei},[],'r','.')
% errorbar(mean(mean(conns,2)),mean(scores{subs,scorei}),std(scores{subs,scorei})/sqrt(length(scores{subs,scorei}))/2,std(scores{subs,scorei})/sqrt(length(scores{subs,scorei}))/2,std(mean(conns,2))/sqrt(length(mean(conns,2)))/2,std(mean(conns,2))/sqrt(length(mean(conns,2)))/2,'r')
% 
% subs=inds_FEPremis;
% conns=Mat(subs,inds);
% scatter(mean(conns,2),scores{subs,scorei},[],'m','.')
% errorbar(mean(mean(conns,2)),mean(scores{subs,scorei}),std(scores{subs,scorei})/sqrt(length(scores{subs,scorei}))/2,std(scores{subs,scorei})/sqrt(length(scores{subs,scorei}))/2,std(mean(conns,2))/sqrt(length(mean(conns,2)))/2,std(mean(conns,2))/sqrt(length(mean(conns,2)))/2,'m')
% 
% subs=inds_HCs;
% conns=Mat(subs,inds);
% scatter(mean(conns,2),scores{subs,scorei},[],'g','.')
% errorbar(mean(mean(conns,2)),mean(scores{subs,scorei}),std(scores{subs,scorei})/sqrt(length(scores{subs,scorei}))/2,std(scores{subs,scorei})/sqrt(length(scores{subs,scorei}))/2,std(mean(conns,2))/sqrt(length(mean(conns,2)))/2,std(mean(conns,2))/sqrt(length(mean(conns,2)))/2,'g')
% 
% subs=[inds_depres_PCs inds_anxiety_PCs];
% conns=Mat(subs,inds);
% scatter(mean(conns,2),scores{subs,scorei},[],'b','.')
% errorbar(mean(mean(conns,2)),mean(scores{subs,scorei}),std(scores{subs,scorei})/sqrt(length(scores{subs,scorei}))/2,std(scores{subs,scorei})/sqrt(length(scores{subs,scorei}))/2,std(mean(conns,2))/sqrt(length(mean(conns,2)))/2,std(mean(conns,2))/sqrt(length(mean(conns,2)))/2,'b')


% subs=inds_depres_PCs;
% conns=Mat(subs,inds);
% scatter(mean(conns,2),scores{subs,scorei},[],'c','.')
% errorbar(mean(mean(conns,2)),mean(scores{subs,scorei}),std(scores{subs,scorei})/sqrt(length(scores{subs,scorei}))/2,std(scores{subs,scorei})/sqrt(length(scores{subs,scorei}))/2,std(mean(conns,2))/sqrt(length(mean(conns,2)))/2,std(mean(conns,2))/sqrt(length(mean(conns,2)))/2,'c')
% 
% subs=inds_anxiety_PCs;
% conns=Mat(subs,inds);
% scatter(mean(conns,2),scores{subs,scorei},[],'b','.')
% errorbar(mean(mean(conns,2)),mean(scores{subs,scorei}),std(scores{subs,scorei})/sqrt(length(scores{subs,scorei}))/2,std(scores{subs,scorei})/sqrt(length(scores{subs,scorei}))/2,std(mean(conns,2))/sqrt(length(mean(conns,2)))/2,std(mean(conns,2))/sqrt(length(mean(conns,2)))/2,'b')

%% symptom correlation twISFC

dataID='10TRwin1TRstep_refconsN15_brainnetome_EPIgroupMask_2mm';

nroi=272;
ncon=29;
npat=36;
nsub=ncon+npat;

consBL=load(['BL_consN' num2str(ncon) '_' dataID '_corrs.mat']);
patsBL=load(['BL_patsN' num2str(npat) '_' dataID '_corrs.mat']);
consFU=load(['FU_consN' num2str(ncon) '_' dataID '_corrs.mat']);
patsFU=load(['FU_patsN' num2str(npat) '_' dataID '_corrs.mat']);

ntw=length(consBL.results);
% subs=1:nsub; % all
% subs=ncon+1:nsub; % only patients

% inds=inds_roi';
inds=find(nbs.NBS.con_mat{1}); % hippocampus network
% inds=find(triu(ones(272),1)); % whole brain
%
for twi=1:ntw
    disp(['time window ' num2str(twi)])
    % create connectivity matrices
    MatBL=zeros(nsub,nroi,nroi);
    MatFU=zeros(nsub,nroi,nroi);
    triu_inds=find(triu(ones(nroi),1));

    MatBL(1:ncon,triu_inds)=consBL.results(twi).corrs';
    MatBL(ncon+1:nsub,triu_inds)=patsBL.results(twi).corrs';
    MatBL=permute(MatBL,[2 3 1]);
    MatBL=MatBL+permute(MatBL,[2 1 3])+repmat(eye(nroi),[1,1,nsub]);
    MatBL=permute(MatBL,[3 1 2]);
    
    % store group mean ISFC in the interaction network of each time window
    temp=MatBL(1:ncon,inds);
    cons_twISFCmean_bl(twi)=mean(mean(temp,2));
    temp=MatBL(ncon+1:nsub,inds);
    pats_twISFCmean_bl(twi)=mean(mean(temp,2));

    MatFU(1:ncon,triu_inds)=consFU.results(twi).corrs';
    MatFU(ncon+1:nsub,triu_inds)=patsFU.results(twi).corrs';
    MatFU=permute(MatFU,[2 3 1]);
    MatFU=MatFU+permute(MatFU,[2 1 3])+repmat(eye(nroi),[1,1,nsub]);
    MatFU=permute(MatFU,[3 1 2]);

    % store group mean ISFC in the interaction network of each time window
    temp=MatFU(1:ncon,inds);
    cons_twISFCmean_fu(twi)=mean(mean(temp,2));
    temp=MatFU(ncon+1:nsub,inds);
    pats_twISFCmean_fu(twi)=mean(mean(temp,2));
    
    MatDiff=MatBL-MatFU;    

    load('main_symptoms_diffBLFU.mat')

    % 
         %   

     subs=1:size(scores,1); % all
    % subs=30:size(scores,1); % patients
%     subs=1:29; % controls

    conns=MatDiff(subs,inds);
    
    conns_patsN26=MatBL(inds_memScoresPats,inds);
    conns_consN29=MatBL(inds_memScoresCons,inds);

    % correlation with memory scores
    [rGmem_patsN26(twi),pGmem_patsN26(twi)]=corr(mean(conns_patsN26,2),memScoresPats);
    [rGmem_consN29(twi),pGmem_consN29(twi)]=corr(mean(conns_consN29,2),memScoresCons);
    
% correlations separately for each significant connections
    for scorei=3:length(scores.Properties.VariableNames)
        if isnumeric(scores{:,scorei}) && nnz(scores{:,scorei})>0
%             [r,p]=corr(conns,scores{subs,scorei});
%             fdrP = mafdr(p,'BHFDR',true);
%             if any(fdrP<0.05)
%                 disp(['significant connection -> ' scores.Properties.VariableNames{scorei}])
%                 significants=find(fdrP<0.05);
%                 for n=significants'
%                     i_lab=nbs.NBS.node_label{i(n)};
%                     j_lab=nbs.NBS.node_label{j(n)};
%                     stat=nbs.NBS.test_stat(i(n),j(n));
%                     fprintf('%s to %s. Test stat: %0.2f\n',i_lab,j_lab,stat);
%                 end        
%             end
            % global correlation, i.e. with the average connectivity in the network
            [rG(scorei,twi),pG(scorei,twi)]=corr(mean(conns,2),scores{subs,scorei});
%             disp(['global correlation -> ' scores.Properties.VariableNames{scorei} ': r=' num2str(rG(scorei,twi)) ', p=' num2str(pG(scorei,twi))])
        end
    end

end
    

%% correlate twISFC vs symptoms correlation time series with ratings

load('twISFC_symptom_BL_correlation_patients_largerHippNetwork')
load('ratings')

for sympti=1:size(rG,1)     
    for rati=1:length(labels)
%         rating=ratings_convHRF_zscored(5:end-6,rati);
        rating=ratings_zscored(5:end-6,rati);
        rR=corr(rG(sympti,:)',rating);
        DF=bramila_autocorr(rG(sympti,:)',rating);
        pR=pvalPearson('b', rR, DF+2);

        disp(['symptom: ' scores.Properties.VariableNames{sympti} ', rating: ' labels{rati} ', r=' num2str(rR) ', p=' num2str(pR)])
    end
end

%% correlate with ratings
load('ratings.mat')

%score='';
% index=find(strcmp(scores.Properties.VariableNames,score));

ISFC=mean([pats_twISFCmean_bl; cons_twISFCmean_bl; pats_twISFCmean_fu; cons_twISFCmean_fu]);
for rati=1:length(labels)
%     rating=ratings_convHRF_zscored(5:end-6,rati);
    rating=ratings_zscored(5:end-6,rati);
%     rR=corr(rG(index,:)',rating);
    rR=corr(ISFC',rating);
    
%     DF=bramila_autocorr(rG(index,:)',rating);  
    DF=bramila_autocorr(ISFC',rating);
    pR=pvalPearson('b', rR, DF+2);
    
    disp(['rating: ' labels{rati} ', r=' num2str(rR) ', p=' num2str(pR)])
end

%% partial correlation
rati=1;
ratiControl=2;
rP=partialcorr(rG(index,:)',ratings_convHRF_zscored(5:end-6,rati),ratings_convHRF_zscored(5:end-6,ratiControl))
%%
% for i=1:length(scores_labels)
%     [r,p]=corr(conns,scores(:,i));
%     fdrP = mafdr(p,'BHFDR',true);
%     if any(fdrP<0.05)
%         disp(['significant connection -> ' scores_labels{i}])
%         significants=find(fdrP<0.05);
%         for n=significants'
%             i_lab=nbs.NBS.node_label{i(n)};
%             j_lab=nbs.NBS.node_label{j(n)};
%             stat=nbs.NBS.test_stat(i(n),j(n));
%             fprintf('%s to %s. Test stat: %0.2f\n',i_lab,j_lab,stat);
%         end        
%     end
%     [rG,pG]=corr(mean(conns,2),scores(:,i));
%     disp(['global correlation -> ' scores_labels{i} ': r=' num2str(rG) ', p=' num2str(pG)])
% end
%%  
% [r,p]=corr(conns,anxiety);
% [r,p]=corr(conns,delusions);
% [r,p]=corr(conns,bprs10);
% global correlation
% [rG,pG]=corr(mean(conns,2),anxiety);
% [rG,pG]=corr(mean(conns,2),delusions);
% [rG,pG]=corr(mean(conns,2),bprs10);


for n=1:length(inds)
    if p(n)<0.05
        i_lab=nbs.NBS.node_label{i(n)};
        j_lab=nbs.NBS.node_label{j(n)};
        stat=nbs.NBS.test_stat(i(n),j(n));
        fprintf('%s to %s; r: %0.2f\n',i_lab,j_lab,r(n));
    end
end

%% regress out nuisance covariates
file='FU_brainnetome_EPIgroupMask_connMat.mat';
load(file)
Mat=permute(Mat,[3 1 2]);
triu_inds=find(triu(ones(size(Mat,2)),1));
data=Mat(:,triu_inds);
load('designMat_2x2anova_FDcpzNuisance.mat')
if contains(file,'BL')
    model=[ones(size(Mat,1),1) design(1:65,69:70)];
elseif contains(file,'FU')
    model=[ones(size(Mat,1),1) design(66:end,69:70)];
end
b=data'/model';
dataReg=data-model*b';
Mat=zeros(size(Mat,1),size(Mat,2),size(Mat,3));
Mat(:,triu_inds)=dataReg;
Mat=permute(Mat,[2 3 1]);
Mat=Mat+permute(Mat,[2 1 3])+repmat(eye(size(Mat,1)),[1,1,size(Mat,3)]);

%% post-hoc t-test2 pats vs cons

load('FU_brainnetome_EPIgroupMask_connMat.mat')
Mat=permute(Mat,[3 1 2]);

pats=Mat(30:65,inds);
cons=Mat(1:29,inds);

% pats_globalConn=mean(pats,2);
% cons_globalConn=mean(cons,2);
pats_globalConn=tanh(mean(atanh(pats),2));
cons_globalConn=tanh(mean(atanh(cons),2));
pats_std=std(pats_globalConn);
cons_std=std(cons_globalConn);
pats_sem = pats_std/sqrt(length(pats_globalConn));
cons_sem = cons_std/sqrt(length(cons_globalConn));

% plot
xtick1=0.2;
xtick2=0.5;

figure
bar(mean(pats_globalConn),'xData',xtick1,'FaceColor',[1 0 0],'BarWidth',0.2);
hold on
errorbar(mean(pats_globalConn),pats_sem,'k','xData',xtick1,'LineWidth',1);

bar(mean(cons_globalConn),'xData',xtick2,'FaceColor',[0 0 1],'BarWidth',0.2);
errorbar(mean(cons_globalConn),cons_sem,'k','xData',xtick2,'LineWidth',1);
% set(gca,'XTick',[xtick1 xtick2],'YTick',[0 .01 .02],'XLim',[0 0.7],'YLim',[0 .02])

% two-sample t-test
[~,p,~,stats]=ttest2(pats_globalConn,cons_globalConn);
