
%% Manuscript Title
%% Co-responding Author1,*, Co-author2 and Co-Author2

%% Department of XXXXXXX, Address XXXX etc.
%% Department of XXXXXXX, Address XXXX etc.
%% Received on XXXXX; revised on XXXXX; accepted on XXXXX 

%% Associate Editor: XXXXXXX

%% ABSTRACT

% Please add abstract here .....

%% 1	INTRODUCTION 

% Please add introduction here .....



%% 2	METHODS

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
N = max(Subject);

% No. of Probe sets
n = size(Data,1);

%Obtain the log 2 centered gene expression data for each subject.
gexp2 = cell(1,N);
Time  = cell(1,N);
for i = 1:N
    ind = find(Subject==i);
    tmp = sort(Pos(ind));
    tmp      = tmp(tmp>-1);  
    [Time{i},~,ix]  = unique(tmp);
    gexp2{i} = nan(n,length(Time{i}));
    gexp2{i} = log2(Data(:,ind(unique(ix))));
    gexp2{i} = gexp2{i}-repmat(nanmean(gexp2{i},2),1,length(Time{i}));
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
ind= 0 ;
for sub = 1:N
    if(sub>9)
        figure(2);
        ind=9;
    end
    if(sub>18)
        figure(3);
        ind=18;
    end
    %subplot(3,3,sub-ind)
    surf(gexp2{sub},'FaceColor','interp','EdgeColor','none')
    xlim([Time{sub}(1),length(Time{sub})])
    set(gca,'XTick',1:length(Time{sub}),'Xticklabel',Time{sub})
    ylim([1,n])
    zlim([min(min(gexp2{sub})),max(max(gexp2{sub}))])
    xlabel('Time')
    ylabel('Probe Set')
    zlabel('Gene Expression Level')
    title(char(Subject_name))
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

h=figure(4);
ind= 0 ;
for sub = 1:N
    if(sub>9)
        h1=figure(5);
        ind=9;
    end
    if(sub>18)
        h2=figure(6);
        ind=18;
    end
 %   subplot(3,3,sub-ind)
    surf(yCR{sub},'FaceColor','interp','EdgeColor','none')
    xlim([Time{sub}(1),length(Time{sub})])
    set(gca,'XTick',1:length(Time{sub}),'Xticklabel',Time{sub})
    ylim([1,size(yCR{sub},1)])
    zlim([min(min(yCR{sub})),max(max(yCR{sub}))])
    xlabel('Time')
    ylabel('Probe Set')
    zlabel('Probe Set Expression Level')
    hold on;
    title(['Subject ',num2str(sub)])
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
    nbasis   = length(Time{i}) + 4;
    norder   = 6;
    basisobj = create_bspline_basis([min(Time{i}) max(Time{i})], nbasis, norder, knots);
    
    %  -----------  Otain optimal smoothing parameter  -----------------
    B              = eval_basis(Time{i},basisobj);
    R              = eval_penalty(basisobj,2);
    lambdagenes{i} = fminbnd(@multiple_GCV_fun, 10.^-6, 10^6, options, B, yCR{i}', R);
    fdParobj       = fdPar(basisobj, 2, lambdagenes{i});
    [fdgenens{i}, dfgenens{i}, gcvgenens{i},~,SSE{i}]  = smooth_basis(Time{i}, gexp2{i}', fdParobj);
    yhat{i}        = eval_fd(Time{i}, fdgenens{i});
    dyhat{i}       = eval_fd(Time{i}, fdgenens{i},1);
    STDERR{i}      = sqrt(sum(SSE{i})/(n*(length(Time{i})-dfgenens{i})));
end


col_hed = {'Df','GCV','log10(\lambda)','Std Error'};
row_hed = strcat(repmat({'Subject '},N,1),cellstr(arrayfun(@num2str, 1:N, 'UniformOutput', false))');
tmp = round2([cell2mat(dfgenens), mean(cell2mat(gcvgenens),2),log10(cell2mat(lambdagenes)),sum(cell2mat(STDERR),2)]);
makeHtmlTable(tmp,[],row_hed,col_hed)


%%
% _*Table 1* The optimal smoothing parameters, the corresponding degrees of freedom,
% the generalized cross validation and the sum of squared errors for the fitted curves $X_{i,j}(t)$ for each subject in the study._

if(ispc)
    %Create Excel Files with the fitted curves and the derivatives of the
    %fitted curves.
    warning('off','all')
    
    for i = 1:N
        xlswrite('Fitted_curves.xls',yhat{i}',i)
        xlswrite('Derivative_Fitted_Curves.xls',dyhat{i}',i)
    end
    
else
    %Still problems with xlwrite on mac 
    for i = 1:N
        xlwrite('Fitted_curves.xlsx',yhat{i}',i);
        %xlwrite('Derivative_Fitted_Curves.xlsx',dyhat{i}',i);
    end
    
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

F = cell(N,1);
INDF = cell(N,1);
yhat       = cell(N,1);
for i = 1:N
    
    F{i} = Ftest(gexp2{i}', Time{i},  fdgenens{i}, dfgenens{i});
    
    %  ---------------  Plot Top genes  ------------------------
    yhat{i}    = eval_fd(Time{i}, fdgenens{i});
    
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

if(ispc)
    warning('off','all')
    
    Col = 'A':'X';
    for i = 1:N
        xlRange = [Col(i) '1'];
        xlswrite('F_value.xls',F{i}',1,xlRange)
        xlswrite('Index_Ftest.xls',INDF{i}',1,xlRange)
        xlswrite('Index_Ftest_DRG.xls',IND_DRG{i}',1,xlRange)
        xlswrite('Probe_set_ID_Ftest_DRG.xls',GID_DRG{i},1,xlRange)
        xlswrite('DRG.xls',DRG{i}',i)
    end
    
else
    
    Col = 'A':'X';
    for i = 1:N
        xlRange = [Col(i) '1'];
        xlwrite('F_value.xls',F{i}',1,xlRange)
        xlwrite('Index_Ftest.xls',INDF{i}',1,xlRange)
        xlwrite('Index_Ftest_DRG.xls',IND_DRG{i}',1,xlRange)
        xlwrite('Probe_set_ID_Ftest_DRG.xls',GID_DRG{i},1,xlRange)
        xlwrite('DRG.xls',DRG{i}',i)
    end
    
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

h=figure(7);
ind= 0 ;
for sub = 1:N
    if(sub>9)
        h1=figure(8);
        ind=9;
    end
    if(sub>18)
        h2=figure(9);
        ind=18;
    end
    %subplot(3,3,sub-ind)
    cutoff =3000;
    surf(yhat{sub}(:,IND_DRG{i}),'FaceColor','interp','EdgeColor','none')
    ylim([1,length(Time{sub})])
    set(gca,'YTick',1:length(Time{sub}),'Yticklabel',Time{sub})
    xlim([1,cutoff])
    zlim([min(min(yhat{sub}(:,IND_DRG{i}))),max(max(yhat{sub}(:,IND_DRG{i})))])
    hold on;
    title(['Subject',num2str(sub)])
    set(gcf,'PaperPosition',[0 0 8.5 11]);
    hold off;
end


%%
% _*Figure 3* The Log 2 Centered Gene Expression Level for each probe set
% contained in the TRG list across all timepoints for each subject._


% The TRG's that are common across all subjects.
if(size(GID_DRG,2)>1)
com_gene = intersect2(GID_DRG);
if(~isempty(com_gene))
    makeHtmlTable(com_gene,[],[],{'Common Probe Sets'})
    %%
    % _*Table 2* The TRGs that are common across all subjects._
end
end

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
    std_data{i}=zscore(gexp2{i}(INDF{i}(1:cutoff),:)')';
    [fidxcluster{i},rmclusters{i},c{i},mean_clusters_mat{i},clusters{i}] = IHC(std_data{i},alpha,1);
    n_clusters{i}   = cellfun(@(x) size(x,1),clusters{i},'UniformOutput', false);
end

for  i = 1:N
for l = 1:length(fidxcluster{i})
        Cluster_IDX{i}(fidxcluster{i}{l}) = l;
end
end

if(ispc)
    for i = 1:N
        xlswrite('Cluster_IDX.xls',Cluster_IDX{i}',i)
    end
else
    for i = 1:N
        xlwrite('Cluster_IDX.xls',Cluster_IDX{i}',i)
    end
end
    
disp(strcat('This is a link to the Cluster Indexs <a href="',flder,'/Cluster_IDX.xls">Cluster_Index</a>.'))

if(cluster_plots)
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
       % print(h8, '-depsc', [num2str(i),num2str(b),'ClusterCenH.eps']);
    end
end

% The temporal gene response modules for each subject are shown in the
% following file:
disp(strcat('This is a link to the Cluster Plots <a href="',flder,'/Cluster.ps">Cluster_Plots</a>.'))

end
    
%%
% Moreover, these temporal gene response modules can be classified into single-gene modules (SGM) with only one gene in each cluster,
% small-size modules (SSM) that contain between 2-10 genes in each cluster, medium-size modules (MSM) that consist of 11-99 genes in each of
% the clusters and large-size modules (LSM) which contain over 100 genes in each cluster.

for k=1:N
    sz{k}       = cell2mat(n_clusters{k});
    ind         = find(sz{k}>99);
    ind1        = find(sz{k}>9 & sz{k}<99);
    ind2        = find(sz{k}>1 & sz{k}<10);
    ind3        = find(sz{k}==1);
    lrg_id{k}   = GID_DRG{k}(vertcat(fidxcluster{k}{ind}));
    med_id{k}   = GID_DRG{k}(vertcat(fidxcluster{k}{ind1}));
    smal_id{k}  = GID_DRG{k}(vertcat(fidxcluster{k}{ind2}));
    sin_id{k}   = GID_DRG{k}(vertcat(fidxcluster{k}{ind3}));
    sizes{k}    = [size(sz{k},1),length(ind),length(ind1),length(ind2),length(ind3)];
end

% *Table 2* provides the number of clusters for each subject and the number of clusters in each category.
col_hed = {'No. of Modules','No. of LSM','No. of MSM','No. of SSM','No. of SGM'};
tmp = round2(vertcat(sizes{:}));
makeHtmlTable(tmp,[],row_hed,col_hed);

% The TRG's that are common across all subjects.
if(N>1)
com_gene_LSM = intersect2(lrg_id);
com_gene_MSM = intersect2(med_id);
com_gene_SSM = intersect2(smal_id);
com_gene_SGM = intersect2(sin_id);

if(~isempty(com_gene_LSM))
    makeHtmlTable(com_gene_LSM,[],[],{'Common Probe Sets for LSM'})
    %%
    % _*Table 2* The TRGs that are common across all subjects for LSM._
end
if(~isempty(com_gene_MSM))
    makeHtmlTable(com_gene_MSM,[],[],{'Common Probe Sets for MSM'})
    %%
    % _*Table 2* The TRGs that are common across all subjects for MSM._
end
if(~isempty(com_gene_SSM))
    makeHtmlTable(com_gene_SSM,[],[],{'Common Probe Sets for SSM'})
    %%
    % _*Table 2* The TRGs that are common across all subjects for SSM._
end
if(~isempty(com_gene_SGM))
    makeHtmlTable(com_gene,[],[],{'Common Probe Sets for SGM'})
    %%
    % _*Table 2* The TRGs that are common across all subjects for SGM._
end
end

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

if(Annotation==1)
  
  for  i = 1:N  
    [tableReport{i},chartReport{i},ClusterReport{i}] = gene_annotation(GID_DRG{i});
  end
  
if(ispc)
for i = 1:N
    struc2xls('Annotation',chartReport{i},'Sheet',i)
end

disp(strcat('This is a link to the Annotation <a href="',flder,'/Annotation.xls">Annotation</a>.'))

end

end
 
% %%
% % *WARNING: DAVID has resource limits so if the gene lists are large it may be very slow. Try again after traffic peak time (10 am - 5 pm EST).*
% 

%% RESULTS

%% ACKNOWLEDGEMENTS

%% REFERENCES
