%% A Pipeline for High-Dimensional Time Course Gene Expression Data
% M. Carey, S. Wu and H. Wu
%
% Here we present a pipeline for high-dimensional time course gene
% expression data. This pipeline utilizes functional data analysis (FDA)
% techniques in order to study the dynamics of gene networks.
%
% The pipeline includes:

%% Preprocessing
%
% Affymetrix Genechip arrays are currently among the most widely used high-throughput technologies for the genome-wide measurement
% of expression profiles. To minimize mis- and cross-hybridization problems, this technology includes both perfect match (PM) and
% mismatch (MM) probe pairs as well as multiple probes per gene (Lipshutz et al., 1999). As a result, significant preprocessing is
% required before an absolute expression level for a specific gene may be accurately assessed. In general, preprocessing probe-level
% expression data consists of three steps: background adjustment (remove local artifacts and ``noise"), normalization (remove array effects),
% and summarization at the probe set level (combine probe intensities across arrays to obtain a measure of the expression level of corresponding mRNA).
% The Pipeline allows the user to select from the following popular preprocessing techniques: Microarray suite 5 (MAS5), Robust Multi-array Average (RMA) and Guanine Cytosine Robust Multi-Array Analysis (GCRMA).


% -----------------------------------------------------------------------
%  Data Preprocessing
% -----------------------------------------------------------------
% No. of Subjects
N = 1;

% No. of Probe sets
n = size(Data,1);

%Obtain the log 2 centered gene expression data for each subject.
gexp = cell(1,N);
%Centered Gene Expression
gexp2 = cell(1,N);
Time  = cell(1,N);
for i = 1:N
    %Get Subject
    ind = 1:length(Subject);
    %Get time
    [tmp,ix] = sort(Pos(ind));
    %Check time > -1
    tmp = tmp(tmp>-1); 
    %Check for replicates
    tb       = tabulate(tmp);
    Time{i}  = tb(:,1);
    nt       = length(Time{i});
    gexp_wr{i}        = zeros(n,max(tb(:,2)),nt);
    for j = 1:nt
    %gexp2 time gene 
    idx               = find(tb(j,1)==tmp);
    gexp_wr{i}(:,1:tb(j,2),j) = log2(Data(:,ind(ix(idx))));
    gexp{i}(:,j)      = nanmean(log2(Data(:,ind(ix(idx)))),2);
    end
    gexp2{i} = gexp{i} - repmat(nanmean(gexp{i},2),1,nt);
end

%%
%
% Lets assume that the centred expression profile of the $i^{th}$ gene, belonging to subject $j$, $X_{i,j}(t)$,
% is a smooth function over time and that the time course gene expression measurements are discrete observations
% from this smooth function, which have been distorted by noise, i.e.,
% $Y_{i,j}(t_{k}) -\mu_{i,j} =X_{i,j}(t_{k})+\epsilon_{i,j}(t_{k})$, for $i=1,\ldots,n,$ $j=1,\ldots,N$ and $k=1,\ldots,K_{i,j},$
% where $n$ is the number of genes, $N$ is the number of subjects, and $K_{i,j}$ is the number of time points observed for
% the $i^{th}$ gene, belonging to subject $j$. The noise denoted by $\epsilon_{i,j}(t_{k})$ is assumed to be an independently
% identically distributed (i.i.d.) random variable with mean $0$ and
% variance $\sigma^{2}$. Figure (1) provides an illustration the time course gene expression measurements over time for each subject.

figure(1);
clear title
ind= 0;
for sub = 1:N
    surf(gexp2{i},'FaceColor','interp','EdgeColor','none')
    xlim([Time{sub}(1),length(Time{sub})])
    set(gca,'XTick',1:length(Time{sub}),'Xticklabel',Time{sub})
    ylim([1,n])
    zlim([min(min(gexp2{i})),max(max(gexp2{i}))])
    xlabel('Time')
    ylabel('Probe Set')
    zlabel('Mean Expression Level')
    title([char(Subject_name),' all genes measured'])
    hold on;
    hold off;
end


%%
% _*Figure 1* Log 2 Centered Gene Expression Level for each probe set
% across all timepoints for each subject._  


%% Spline Smoothing

%%
%
% The pipeline incorporates the option to obtain the functional entity $X_{i,j}(t)$ by spline smoothing (Ramsay and Silverman 2005). Spline smoothing involves
% minimizing the following penalized residual sum of squares
%
% $$(Y_{i,j}(t_{k}) -\mu_{i,j} - X_{i,j}(t))^{2} + \lambda \int \left[D^{2} X_{i,j}(t) \textrm{d}t\right]^{2}.$$
%
% The first term defines the discrepancy between the observed centred gene expression measurements
% $Y_{i,j}(t_{k}) -\mu_{i,j}$ and the functional entity $X_{i,j}(t),$ while the second term requires $X_{i,j}(t)$ to be sufficiently smooth.
% The parameter $\lambda$ controls the trade-off between these two competing effects and hence ensures that $X_{i,j}(t)$ has an appropriate
% amount of regularity or smoothness. As smoothing each individual gene is computationally expensive and is likely to obtain biased estimates
% as our data does not contain replicates and is expected to contain a low signal-to-noise ratio. We follow an approach similar to
% (Yao et al. 2005, Wu and Wu, 2013) and pool the genes together in order to obtain an estimate of the amount of regularity required.
%
% As the majority of the genes have a relatively flat expression profile over time (see Figure 1) estimating the smoothing parameter by pooling all
% the genes together is not ideal as it will tend to over smooth the data. Therefore, we choose a subset of the genes that exhibit time course
% patters that have relatively smooth trajectories that do not fluctuate widely. Then we rank these genes by their interquartile range and select
% the top genes for our estimation subset. Figure 2 illustartes the
% genes expresion curves in the estimation subset over time for each subject.

%  -----------------------------------------------------------------------
%      Select the estimation subset for the smoothing parameter selection
%  -----------------------------------------------------------------------

yCR = Est_Sub_Sel(Time,gexp2,N);

h=figure(2);
ind= 0 ;
for sub = 1:N
    surf(yCR{sub},'FaceColor','interp','EdgeColor','none')
    xlim([Time{sub}(1),length(Time{sub})])
    set(gca,'XTick',1:length(Time{sub}),'Xticklabel',Time{sub})
    ylim([1,size(yCR{sub},1)])
    zlim([min(min(yCR{sub})),max(max(yCR{sub}))])
    xlabel('Time')
    ylabel('Probe Set')
    zlabel('Probe Set Expression Level')
    hold on;
     title([char(Subject_name),' subset of genes with smooth',...
         ' trajectories with dof<5'])
    set(gcf,'PaperPosition',[0 0 8.5 11]);
    hold off;
end

%%
% _*Figure 2* The Log 2 Centered Gene Expression Level for each probe set
% in the estimation subset across all timepoints for each subject._

%%
% The smoothing parameter is estimated by the conventional generalized cross validation approach using the genes in our estimation subset
% (see Figure 2). Once the smoothing parameter is selected one can proceed to the estimate the functional entity $X_{i,j}(t)$ using all gene
% expression profiles. Table 1 provides the optimal smoothing parameters, the corresponding degrees of freedom,
% the generalized cross validation and the sum of squared errors for the fitted curves $X_{i,j}(t)$ for each subject in the study respectively.

%  -----------------------------------------------------------------------
%                        FDA
%  -----------------------------------------------------------------------
fdgenens    = cell(N,1);
dfgenens    = cell(N,1);
gcvgenens   = cell(N,1);
lambdagenes = cell(N,1);
yhat        = cell(N,1);
dyhat       = cell(N,1);
STDERR = cell(N,1);
SSE = cell(N,1);

options = optimset('LargeScale',  'off',  ...
    'Display',     'on',   ...
    'Diagnostics', 'off', ...
    'GradObj',     'off',     ...
    'Hessian',     'off',    ...
    'TolFun',      1e-8, ...
    'TolX',        1e-8);

for i = 1:N
    %  ---------------  set up the b-spline basis  ------------------------
    knots    = Time{i}';
    norder   = 3;
    nbasis   = length(Time{i}) + norder - 1;
    basisobj = create_bspline_basis([min(Time{i}) max(Time{i})], nbasis, norder, knots);
    
    %  -----------  Otain optimal smoothing parameter  -----------------
    B              = eval_basis(Time{i},basisobj);
    R              = eval_penalty(basisobj,2);
    lambdagenes{i} = fminbnd(@multiple_GCV_fun, 10.^-6, 10^6, options, B, yCR{i}', R);
    fdParobj       = fdPar(basisobj, 2, lambdagenes{i});
    [fdgenens{i}, dfgenens{i}, gcvgenens{i},~,SSE{i}]...  
                   = smooth_basis(Time{i}, gexp2{i}', fdParobj);
    yhat{i}        = eval_fd(Time{i}, fdgenens{i});
    dyhat{i}       = eval_fd(Time{i}, fdgenens{i},1);
    STDERR{i}      = sqrt(sum(SSE{i})/(n*(length(Time{i})-dfgenens{i})));
end

h=figure(3);
ind= 0 ;
for sub = 1:N
    surf(yhat{sub}','FaceColor','interp','EdgeColor','none')
    xlim([Time{sub}(1),length(Time{sub})])
    set(gca,'XTick',1:length(Time{sub}),'Xticklabel',Time{sub})
    ylim([1,size(yhat{sub},2)])
    zlim([min(min(yhat{sub})),max(max(yhat{sub}))])
    xlabel('Time')
    ylabel('Probe Set')
    zlabel('Probe Set Expression Level')
    hold on;
     title([char(Subject_name),' smooth gene expresion curves'])
    set(gcf,'PaperPosition',[0 0 8.5 11]);
    hold off;
end

col_hed = {'Df','GCV','log10(\lambda)','Std Error'};
row_hed = strcat(repmat({'Subject '},N,1),cellstr(arrayfun(@num2str, 1:N, 'UniformOutput', false))');
tmp = round2([cell2mat(dfgenens), mean(cell2mat(gcvgenens),2),log10(cell2mat(lambdagenes)),sum(cell2mat(STDERR),2)]);
makeHtmlTable(tmp,[],row_hed,col_hed);

% Caption = 'The degrees of freedom (Df), the generalized cross validation (GCV), the smoothing paramter $\lambda$ and the standard error of the fitted curves produced by spline smoothing'; 
% Label = 'fit1';
% filename = [con,'fitted_vals.tex'];
% Generate_Latex_Tabels(tmp,row_hed,col_hed,Caption,Label,filename)


%%
% _*Table 1* The optimal smoothing parameters, the corresponding degrees of freedom,
% the generalized cross validation and the sum of squared errors for the fitted curves $X_{i,j}(t)$ for each subject in the study._

    %Create Excel Files with the fitted curves and the derivatives of the
    %fitted curves.   
    for i = 1:N
        create_exel_file('Fitted_curves.xls',yhat{i}',i,[],path);
        create_exel_file('Derivative_Fitted_Curves.xls',dyhat{i}',i,[],path);
    end
    
disp(strcat('This is a link to the Fitted Curves <a href="',flder,'/Fitted_curves.xls">Fitted_Curves</a>.'))

disp(strcat('This is a link to the Derivatives of the Fitted Curves <a href="',flder,'/Derivative_Fitted_Curves.xls">Derivative_Fitted_Curves</a>.'))


%% Hypothesis Testing to Idenify the Dynamic Responce Genes (DRGs)
%
% In a typical gene expression experiment, tens of thousands of genes are measured simultaneously, but only a fraction
% of them are associated with the biological process of interest or a particular stimulus, such as a
% therapeutic treatment or virus infection. Since it is reasonable to
% include only these responsive
% genes in our analysis, here we identify dynamic response genes (DRGs), i.e., genes with expressions that have changed significantly over time.
%
% The pipeline uses the following hypothesis to identify the DRGs,
%
% $H_{0}$: $X_{i,j}(t)$ equal to $0$
%
% $H_{a}$: $X_{i,j}(t)$ not equal to $0$.
%
% The test statistic is the conventional F-statistic, which compares the goodness-of-fit of the null model to the alternative model:
%
% $F_{i,j}=\frac{RSS_{i,j}^0-RSS_{i,j}^1}{RSS_{i,j}^1},$
%
% where $RSS_{i,j}^0= (Y_{i,j}(t_{k}) -\mu_{i,j} )^{2}$ and $RSS_{i,j}^1=(Y_{i,j}(t_{k}) -\mu_{i,j} - X_{i,j}(t))^{2}$ are the
% residual sum of squares under the null and the alternative models for the $i$-th gene, belonging to subject $j$. The genes with large F-ratios
% can be considered as exhibiting notable changes with respect to time.

F    = cell(N,1);
INDF = cell(N,1);

for i = 1:N
    
    F{i} = Ftest(gexp2{i}, Time{i},  fdgenens{i}, dfgenens{i});
    
    [SF, INDF{i}] = sort(F{i},'descend');
    
end

%Create Excel Files with the F-values, Index for Ranking of the F-values, Probe set IDs for
%the DRGs and the time-course expresion level for the DRGs.

for i = 1:N
        cutoff =3000;
        IND_DRG{i} = INDF{i}(1:cutoff);
        GID_DRG{i} = gid(IND_DRG{i});
        DRG{i}= gexp2{i}(IND_DRG{i},:)';
end

Col = 'A':'X';
    for i = 1:N
        xlRange = [Col(i) '1'];
        create_exel_file('F_value.xls',F{i},1,xlRange,path);
        create_exel_file('Index_Ftest.xls',INDF{i},1,xlRange,path);
        create_exel_file('Index_Ftest_DRG.xls',IND_DRG{i},1,xlRange,path);
        create_exel_file('Probe_set_ID_Ftest_DRG.xls',GID_DRG{i},1,xlRange,path);
        create_exel_file('DRG.xls',DRG{i}',i,[],path);
    end
    
disp(strcat('This is a link to the F statistics <a href="',flder,'/F_value.xls">F_value</a>.'))

disp(strcat('This is a link to the Index F statistics <a href="',flder,'/Index_Ftest.xls">Index_Ftest</a>.'))

disp(strcat('This is a link to the Index of the DRGs <a href="',flder,'/Index_Ftest_DRG.xls">Index_DRGs</a>.'))

disp(strcat('This is a link to the Probe set Ids for TRGs <a href="',flder,'/Probe_set_ID_Ftest_DRG.xls">Probe_sets_DRGs</a>.'))

disp(strcat('This is a link to the DRG values <a href="',flder,'/DRG.xls">Index_DRG</a>.'))


%%
% As we wish to have an equal amount of DRGs for each subject we rank the F-ratios and select 3000 of the top ranking genes (TRG's).
% These genes are shown in Figure 3.

%  ---------------  Surfaces of Top Genes from F-test  ------------------------

h=figure(4);
ind= 0 ;
for sub = 1:N
    surf(yhat{sub}(:,IND_DRG{i}),'FaceColor','interp','EdgeColor','none')
    ylim([1,length(Time{sub})])
    set(gca,'YTick',1:length(Time{sub}),'Yticklabel',Time{sub})
    xlim([1,cutoff])
    zlim([min(min(yhat{sub}(:,IND_DRG{i}))),max(max(yhat{sub}(:,IND_DRG{i})))])
    hold on;
    title([char(Subject_name),' Dynamic Response Genes'])
    set(gcf,'PaperPosition',[0 0 8.5 11]);
    hold off;
end


%%
% _*Figure 3* The Log 2 Centered Gene Expression Level for each probe set
% contained in the TRG list across all timepoints for each subject._


%% Functional Annotation of the genes in the temporal gene response modules
%
% Using the DAVID Functional Annotation Tool devloped by (Da Wei Huang and
% Lempicki, 2008) we produce the following four annotaion reports.
%
% # Table Report
% # Chart Report
% # Cluster Report
%
% These Cluster Report provides the gene-enrichment analysis, pathway mapping, gene/term similarity search, homologue match, ID translation, etc.;
%
% _Table Report_ provides the ID,Gene Name,Species,BBID,BIOCARTA,COG_ONTOLOGY,GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART
% and SP_PIR_KEYWORDS for all the genes inputed.
%
% _Chart Report_
% is an annotation-term-focused view which lists annotation terms and their associated genes under study.
% The threshold of EASE Score, a modified Fisher Exact P-Value, for gene-enrichment analysis. It ranges from 0 to 1.
% Fisher Exact P-Value = 0 represents perfect enrichment. Usually P-Value is equal or smaller than 0.05 to be considered strongly enriched in the annotation categories. Default is 0.1.
%
% _Cluster Report_
% To reduce the redundancy, the Functional Annotation Clustering report groups/displays similar annotations together which makes the biology clearer and more
% focused to be read vs. traditional chart report. The grouping algorithm is based on the hypothesis that similar annotations should have similar gene members.
% The Functional Annotation Clustering integrates the same techniques of  Kappa statistics to measure the degree of the common genes between two annotations, and
% fuzzy heuristic clustering (used in Gene Functional Classification Tool ) to classify the groups of similar annotations according kappa values. In this sense, the more common genes
% annotations share, the higher chance they will be grouped together.
   
%   for  i = 1:N  
%     [tableReport{i},chartReport{i},ClusterReport{i}] = gene_annotation(GID_DRG{i});
%   end

% if(ispc)
% for i = 1:N
%     struc2xls('Annotation',chartReport{i},'Sheet',i)
% end
% 
% disp(strcat('This is a link to the Annotation <a href="',flder,'/Annotation.xls">Annotation</a>.'))
% 
% end
% 
% end
%  
% %%
% % *WARNING: DAVID has resource limits so if the gene lists are large it may be very slow. Try again after traffic peak time (10 am - 5 pm EST).*

% https://ashokragavendran.wordpress.com/2015/10/09/fixing-rdavidwebservice-on-yosemite/

%% Cluster the DRG's into temporal gene response modules


%
% As many of the DRG's exhibit similar expression patterns over time we
% wish to cluster them into temporal gene response modules
% (groups of genes which have similar gene expression profiles over time). In order to achieve this we adopt the Iterative Hierarchal
% Clustering (IHC) method introduced in (Carey et al. 2015). The Iterated Hierarchical clustering (IHC) algorithm can identify inhomogeneous
% clusters, capture both the large and very small clusters and provides an automated selection of the optimal number of clusters. The IHC
% algorithm is only reliant on a single parameter $\alpha,$ which controls the trade-off for the between- and within-cluster correlations and
% is outlined below:

%%
%
% #  Initialization: Cluster the gene expression curves using hierarchical clustering. Let the distance metric be the Spearman correlation with a threshold of $\alpha$.
% #  Merge: Treat each of the cluster centres (exemplars) as `new genes', use the same rule as in Step 1 to merge the exemplars into new clusters.
% #  Prune: Let $c_{i}$ be the centre of cluster $i$. If the correlation between the cluster centre and gene $j$, which will be denoted by $\rho_{i,j},$ is less than $\alpha$ remove $gene_{j}$ from the cluster $i$. Let $M$ be the number of genes removed from the existing $N$ clusters. Assign all $M$ genes into single-clusters. Hence there is now $(N+M)$ clusters in total.
% #  Repeat Steps 2 and 3 until either the index of clusters fully or partially converges or the number of iterations reaches the maximum number of iterations.
% #  Repeat Step 2 until the between-cluster correlations are less than $\alpha$.
%

%  -----------------------------------------------------------------------
%                       Cluster (IHC)
%  -----------------------------------------------------------------------

%Theshold
alpha = 0.75;

for  i = 1:N
    std_data{i}     = zscore(gexp2{i}(INDF{i}(1:cutoff),:)')';
    [fidxcluster{i},rmclusters{i},c{i},mean_clusters_mat{i},clusters{i}] = IHC(std_data{i},alpha);
    n_clusters{i}   = cellfun(@(x) size(x,1),clusters{i},'UniformOutput', false);
end

for  i = 1:N
for l = 1:length(fidxcluster{i})
        Cluster_IDX{i}(fidxcluster{i}{l}) = l;
end
end


    for i = 1:N
        create_exel_file('Cluster_IDX.xls',Cluster_IDX{i}',i,[],path);
    end
    
disp(strcat('This is a link to the Cluster Indexs <a href="',flder,'/Cluster_IDX.xls">Cluster_Index</a>.'))


for i=1:N
    [s,ind]=sort(cell2mat(n_clusters{i}),'descend');
    for b = 1:floor(size(mean_clusters_mat{i},1)./30)
        h8=figure(b);
        for gen = 1:30
            subplot(5,6,gen)
            ind2 = ind(gen+(b-1).*30);
            plot(clusters{i}{ind2}','-*b')
            xlabel('Time')
            ylabel('Probe Set Expression Level')
            hold on;
            plot(mean_clusters_mat{i}(ind2,:),'o-r','LineWidth',1.5)
            xlim([0,size(mean_clusters_mat{i}(ind2,:),2)])
            ylim([min(min(clusters{i}{ind2}))-.05,max(max(clusters{i}{ind2}))+.05])
            v = axis;
            handle=title([' No. Probe Sets ',num2str(s(gen+(b-1).*30))]);
            set(handle,'Position',[2.5 v(4)*1. 0]);
            hold off;
        end
        print(h8,'-dpsc2', '-append', 'Cluster.ps');
    end
end

% The temporal gene response modules for each subject are shown in the
% following file:
disp(strcat('This is a link to the Cluster Plots <a href="',flder,'/Cluster.ps">Cluster_Plots</a>.'))
    
%%
% Moreover, these temporal gene response modules can be classified into single-gene modules (SGM) with only one gene in each cluster,
% small-size modules (SSM) that contain between 2-10 genes in each cluster, medium-size modules (MSM) that consist of 11-99 genes in each of
% the clusters and large-size modules (LSM) which contain over 100 genes in each cluster.

for k=1:N
    sz{k}       = cell2mat(n_clusters{k});
    ind         = find(sz{k}>99);
    ind1        = find(sz{k}>9 & sz{k}<100);
    ind2        = find(sz{k}>1 & sz{k}<10);
    ind3        = find(sz{k}==1);
    lrg_id{k}   = GID_DRG{k}(vertcat(fidxcluster{k}{ind}));
    med_id{k}   = GID_DRG{k}(vertcat(fidxcluster{k}{ind1}));
    smal_id{k}  = GID_DRG{k}(vertcat(fidxcluster{k}{ind2}));
    sin_id{k}   = GID_DRG{k}(vertcat(fidxcluster{k}{ind3}));
    lrg_ts{k}   = mean_clusters_mat{k}(ind,:);
    med_ts{k}   = mean_clusters_mat{k}(ind1,:);
    smal_ts{k}  = mean_clusters_mat{k}(ind2,:);
    sin_ts{k}   = mean_clusters_mat{k}(ind3,:);
    sizes{k}    = [size(sz{k},1),length(ind),length(ind1),length(ind2),length(ind3)];
end

figure(floor(size(mean_clusters_mat{i},1)./30)+1);
subplot(2,2,1)
if(~isempty(lrg_ts{1}))
ribbon(lrg_ts{1}');
ylim([1,size(lrg_ts{1},2)])
xlim([1,size(lrg_ts{1},1)])
zlim([min(min(lrg_ts{1})),max(max(lrg_ts{1}))])
ylabel('Time (hours)')
xlabel('ith Cluster Center')
title('LSM')
end
subplot(2,2,2)
ribbon(med_ts{1}');
ylim([1,size(med_ts{1},2)])
xlim([1,size(med_ts{1},1)])
zlim([min(min(med_ts{1})),max(max(med_ts{1}))])
ylabel('Time (hours)')
xlabel('ith Cluster Center')
title('MSM')
subplot(2,2,3)
ribbon(smal_ts{1}');
ylim([1,size(smal_ts{1},2)])
xlim([1,size(smal_ts{1},1)])
zlim([min(min(smal_ts{1})),max(max(smal_ts{1}))])
ylabel('Time (hours)')
xlabel('ith Cluster Center')
title('SSM')
subplot(2,2,4)
ribbon(sin_ts{1}');
ylim([1,size(sin_ts{1},2)])
xlim([1,size(sin_ts{1},1)])
zlim([min(min(sin_ts{1})),max(max(sin_ts{1}))])
ylabel('Time (hours)')
xlabel('ith Cluster Center')
title('SGM')

% *Table 2* provides the number of clusters for each subject and the number of clusters in each category.
col_hed = {'No. of Modules','No. of LSM','No. of MSM','No. of SSM','No. of SGM'};
row_hed = strcat(repmat({'Subject '},N,1),cellstr(arrayfun(@num2str, 1:N, 'UniformOutput', false))');
tmp = round2(vertcat(sizes{:}));
makeHtmlTable(tmp,[],row_hed,col_hed);

% Caption = 'The number of clusters for each subject and the number of clusters in each category.'; 
% Label = 'clu1';
% filename = [con,'cluster_cat.tex'];
% Generate_Latex_Tabels(tmp,row_hed,col_hed,Caption,Label,filename)

%%
% %% Construct high-dimensional gene regulation networks (GRNs) using differential equation models
% %
% % Here we utizie the ordinary differential equation (ODE) modeling approach in order to reconstruct the high-dimensional gene regulation networks (GRN).
% % In ODE network models, gene regulations are modeled by rate equations, which quantify the rate of change (derivative) of the expression level
% % of one temporal gene response module in the system as a function of expression levels of all related temporal gene response modules.
% % Both up and down regulations as well as self-regulations can be appropriately captured by the ODE model.
% % The general form of the ODE model can be written as
% %
% % $\frac{d(M_{q,j}}{dt} = \alpha_{0,q,j} + \sum_{p=1}^{P} \alpha_{p,q,j}M_{p,j},$
% %
% % for i =q,..,Q where $M_{q,j}$ is the $q^{th}$ gene response module, belonging to subject $j$, $\alpha_{0,q,j}$ is the intercept for the $q^{th}$
% % gene response module, belonging to subject $j$ and the coefficients $\{\alpha_{p,q,j}\}_{i=1}^{P}$ quantify the regulation effects of other
% % gene response modules on the rate of expression change of the q-th gene response module, belonging to subject $j$.
% %
% % Although the model size Q is comparably smaller than that of the original model N,
% % simultaneous model selection and parameter optimization of ODE parameters $\{\alpha{p,q,j}\}_{q,p=1}^{Q,P}$ are still computationally very expensive,
% % because it involves costly numerical integration and complicated parameter regularization.
% % The two-stage smoothing-based estimation method (Voit and Almeida, 2004; Liang and Wu, 2008), decouples the system of
% % differential equations into a set of pseudo-regression models.  This method avoids numerically solving the differential equations directly and 
% % does not require the initial or boundary conditions of the state variables. More importantly, it allows us to perform model selection
% % and parameter estimation for one equation at a time, which significantly reduces the computational
% % cost. 
% %
% % _Step one_
% %
% % # Obtain the estimates of the temporal gene response modules $M_{q,j}$ and their derivatives $\frac{d(M_{q,j}}{dt}$ using the smoothing splines estimates obtained in step one.
% %
% % _Step two_
% %
% % We plug the estimated mean expression curves $\hat{M}_{q,j}$ and their derivatives $\frac{d(\hat{M}_{q,j}}{dt}$ into the ODE model to obtain the following set of pseudo linear regression models.
% %
% % $\frac{d(\hat{M}_{q,j}}{dt} = \beta_{0,q,j} + \sum_{p=1}^{P} \beta_{p,q,j}\hat{M}_{p,j} + \epsilon_{p,j}$
% %
% % for q=1,...,Q. For linear regression models, many penalized methods have been proposed in the regularization framework to conduct variable selection and estimation, such as the
% % least absolute shrinkage and selection operator (LASSO) (Tibshirani,
% % 1996), smoothly clipped absolute deviation (SCAD) (Fan and Li, 2001) and
% % so on. Here we use LASSO regularization to conduct variable selection and estimation using the the follwing following penalized objective function
% %
% % $\epsilon_{p,j} + \lambda \sum_{p=0}^{P} \| \beta_{p,q,j} \|$
% %
% % where $\sum_{p=0}^{P} \| \beta_{p,q,j} \|$ is the LASSO penalty which
% % shrinks the parameters $\beta_{p,q,j}$ to zero and $\lambda$ is the sparisty parameter which controls the trade-off between minimizing the error of the
% % pseudo linear regression model and requiring the parameters $\beta_{p,q,j}$ to be zero. A large $\lambda$ sets all the parameters $\beta_{p,q,j}$ to zero.
% % A low $\lambda$ minmizes the error of the pseudo linear regression model.
% %
% % _Refinement Step_
% %
% % The two-stage smoothing-based estimation method is employed to simplify the computation
% % of the ODE models and also facilitate the variable selection procedure. However, the
% % parameter estimates from the two-stage method are not efficient in terms of estimation accuracy,
% % as there can be considerable approximation error in the estimates of the
% % modules expression curves $\hat{M}_{q,j}$ and their derivatives
% % $\frac{d(\hat{M}_{q,j}}{dt}$ and decoupling the system of
% % differential equations into a set of pseudo-regression models only accounts for the direct effect, that is, how the $q^{th}$ gene response module
% % is regulated by the $p^{th}$ gene response module (q $\rightarrow$ p) and does not account for the indirect effect, that is, how the $r^{th}$ gene response module
% % is regulated by the $q^{th}$ gene response module which in turn regululates the $p^{th}$ gene response module (r $\rightarrow$ q $\rightarrow$ p) which 
% % is a fundamental property of any system of differential equations.
% 
% % To overcome the estimation deficiency and the decoupling effect of the two-stage method, we propose to refine the parameter estimates for the
% % selected ODE model using the nonlinear least squares (NLS) method. The parameter estimates
% % from the two-stage method are used to estimate the sturcture of the
% % estimates only, that is, to inform us which parameters are non-zero.
% %
% %  -----------------------------------------------------------------------
% %                      Two - Stage
% %  -----------------------------------------------------------------------
% 
% %Obtain Smoothed Estimates of the derivative and trajectory.
for i = 1:N
for j = 1:length(fidxcluster{i})
group               = IND_DRG{i}(fidxcluster{i}{j});
meanfd              = mean_grouped(fdgenens{i},group);
TimeEx{i}           = linspace(Time{i}(1),Time{i}(end),100)';
yhatEx{i}(:,j)      = eval_fd(TimeEx{i}, meanfd);
dyhatEx{i}(:,j)     = eval_fd(TimeEx{i}, meanfd,1);
end
end

%Obtain LASSO estimate of the parameters.
EAS   = cell(N);
Stats = cell(N);
for i = 1:N
for j = 1:size(dyhatEx{i},2)
 [EAS{i}(:,j),Stats{i}{j}]     = lasso(yhatEx{i},dyhatEx{i}(:,j));
end       
G{i}     = (EAS{i}~=0); 
A0{i}    = EAS{i}(G{i});
end

for i = 1:N
create_exel_file('Networks.xls',EAS{i}',i,[],path);
end

disp(strcat('This is the parameters of the ODE obtained by two-stage method <a href="',flder,'/Networks.xls">Network</a>.'))


% Obtain Refined estimates of the parameters.
%  optim_options = optimset('Display', 'iter','Algorithm','levenberg-marquardt','TolFun',1.0000e-08,'TolX', 1.0000e-08);
% % 
 A = cell(N);
 for i = 1:N 
  A{i} = EAS{i};%lsqnonlin(@rss_sp,A0{i},[],[],optim_options,TimeEx{i},yhatEx{i},G{i});
 end 

for i = 1:N
create_exel_file('Networks_Refined.xls',A{i}',i,[],path);
end

disp(strcat('This is the estimated parameters of the ODE <a href="',flder,'/Networks_Refined.xls">Parameters</a>.'))

% %% Obtain Network Analysis of the gene regulation networks (GRNs).
% 
% 
% %   Graph Statistics
% %   ----------------
% %
% %       graph_clustercoeff       - Calculate overall clustering coefficient of the graph.
% %       graph_diameter           - Calculate diameter of the graph.
% %       graph_meandist           - Calculate mean distance of the graph.
% %       graph_density            - Calculate density of the graph.
 for i = 1:N
 GS(i,:) = [graph_clustercoeff(sparse(G{i})),graph_diameter(G{i}),graph_meandist(G{i}),graph_density(G{i})];
 end
 
for i = 1:N
create_exel_file('Graph Statistics.xls',GS,i,[],path);
end
 
col_hed = {'Subject','Global clustering coefficient','Diameter of the graph','Mean distance of the graph','Density of the graph'};
row_hed = strcat(repmat({'Subject '},N,1),cellstr(arrayfun(@num2str, 1:N, 'UniformOutput', false))');
makeHtmlTable(GS,[],row_hed,col_hed);
 
% Caption = 'The graph statistics for the GRN.'; 
% Label = 'grn1';
% filename = [con,'GRN_metrics.tex'];
% Generate_Latex_Tabels(GD,row_hed,col_hed,Caption,Label,filename)
 
% %   Node Statistics
% %   ---------------
% %
% %       bridging_centrality      - Calculate Bridging Centrality.
% %       closeness_centrality     - Calculate Closeness Centrality.
% %       eccentricity_centrality  - Calculate Eccentricity Centrality.
% %       delta_centrality         - Calculate Delta Centrality.
% %       current_info_flow        - Calculate Current Information Flow.
% %       assortativitycoeff       - Calculate Assortativity Coefficient.
% %       bridgingcoeff            - Calculate Bridging Coefficient.
% %       brokeringcoeff           - Calculate Brokering Coefficient.
% %       clusteringcoeff          - Calculate Clustering Coefficient.
% %       order2cc                 - Calculate Clustering Coefficient considering 2nd level of connections.
% %       hierarchy                - Calculate Degree of Hierarchy.
% %       kclique                  - Caldulate K-clique.
% %       kcore                    - Calculate K-core of network.
% %       locavgcon                - Calculate Local Average Connectivity.
% %       neighborhood_conn        - Calculate Neighborhood Connectivity.      
% %       participationcoeff       - Calculate Participation Coefficient.
% %       richclubcoeff            - Calculate Rich Club Coefficient.
% %       smallworldindex          - Calculate Small World Index.
% %       within_module_deg        - Calculate Within Module Degree.
% %
  for i = 1:N
  NS{i} = [bridging_centrality(G{i}),closeness_centrality(sparse(G{i})),eccentricity_centrality(sparse(G{i}))'];
  SWI(i) = smallworldindex(G{i});
  end

 for i = 1:N
create_exel_file('Node Statistics.xls',GS,i,[],path);
end 
  
  
col_hed = {'Subject','Bridging Centrality','Closeness Centrality','Eccentricity Centrality'};
row_hed = strcat(repmat({'Subject '},N,1),cellstr(arrayfun(@num2str, 1:N, 'UniformOutput', false))');
makeHtmlTable(GS,[],row_hed,col_hed); 

% Caption = 'The node statistics for the GRN.'; 
% Label = 'grn2';
% filename = [con,'GRN_metrics2.tex'];
% Generate_Latex_Tabels(GS,row_hed,col_hed,Caption,Label,filename)
 
% % % %   Visualization
% % % %   -------------
% % % %
% % % %       plotnet                  - View network in a plot.
% % % %       plotnet_edgewidth        - View network in a plot with defined edge width.
% % % %       plotnet_curve            - View network in a plot with curved edges.
% % % %       plotnet_treering         - View network in a treering plot.
% % % %       viewnetsvg               - View network using SVG.
% % % %       powerlawplot             - Plot relationship between degree and number of nodes.
% % % %       sbe_layout               - Gateway function for layout methods.

%View network in a plot
% options.sym = 1;
% h=figure(13);
% ind= 0 ;
% for sub = 1:N
%     if(sub>9)
%         h1=figure(14);
%         ind=9;
%     end
%     if(sub>18)
%         h2=figure(15);
%         ind=18;
%     end
%     cytoscaperun(G{1},reshape(EAS{1},[size(EAS{1},1).*size(EAS{1},1),1]))
%     title([char(Subject_name),' GRN'])
% end
% % 

